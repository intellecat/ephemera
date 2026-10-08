#!/usr/bin/env python3
"""Pack geocentric positions from a JPL DE kernel into a compact Chebyshev
table (our own format, public-domain source).

Each body is represented as piecewise Chebyshev polynomials of the
geocentric ecliptic-J2000 rectangular coordinates (x, y, z) in AU.  The
intervals are chosen adaptively: the longest interval whose fit stays
within the target error is used, so slow regions get long segments and
fast regions short ones.  Coefficients are quantized per segment and
coordinate to 24-bit integers with a float32 scale.

The result is written to data/ephem.bin; it is loaded at runtime by
src/compressed.cpp.  Development tool: needs numpy + skyfield + pyerfa.

    /tmp/eph/bin/python tools/pack_ephemeris.py [--target 0.01]
"""
from __future__ import annotations

import argparse
import math
import os
import struct
import sys

import erfa
import numpy as np
from skyfield.api import load

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = os.path.join(ROOT, "data", "ephem.bin")


# (eph_body id, DE440s key, frame, sample dt, degree, max interval days)
BODIES = [
    (2, "earth", "helio", 0.5, 12, 128),
    (10, "moon", "geo", 0.125, 16, 32),
    (0, "mercury barycenter", "helio", 0.25, 16, 64),
    (1, "venus barycenter", "helio", 0.25, 14, 64),
    (3, "mars barycenter", "helio", 0.5, 14, 64),
    (4, "jupiter barycenter", "helio", 1.0, 18, 128),
    (5, "saturn barycenter", "helio", 1.0, 16, 128),
    (6, "uranus barycenter", "helio", 1.0, 14, 256),
    (7, "neptune barycenter", "helio", 1.0, 20, 512),
    (8, "pluto barycenter", "helio", 1.0, 20, 512),
]


def sample(key, frame, jd0, jd1, dt):
    eph = load("de440s.bsp")
    ts = load.timescale()
    center = eph["earth"] if frame == "geo" else eph["sun"]
    jds = np.arange(jd0, jd1, dt)
    t = ts.tt_jd(jds)
    h = eph[key].at(t) - center.at(t)
    ra, dec, dist = h.radec()
    ra = np.asarray(ra.radians)
    dec = np.asarray(dec.radians)
    r = np.asarray(dist.au)
    eps0 = erfa.obl80(2451545.0, 0.0)
    x = r * np.cos(dec) * np.cos(ra)
    y = r * np.cos(dec) * np.sin(ra)
    z = r * np.sin(dec)
    X = np.stack([x,
                  y * math.cos(eps0) + z * math.sin(eps0),
                  -y * math.sin(eps0) + z * math.cos(eps0)], axis=1)
    return jds, X


def seg_fit(X, i, j, N):
    """Fit samples [i, j] (inclusive) with degree N.  Returns (maxerr, C)."""
    Y = X[i:j + 1]
    M = len(Y)
    if M < N + 3:
        return 1e9, None
    uu = 2.0 * np.arange(M) / (M - 1) - 1.0
    V = np.polynomial.chebyshev.chebvander(uu, N)
    C, _, _, _ = np.linalg.lstsq(V, Y, rcond=None)
    rec = V @ C
    cr = np.linalg.norm(np.cross(rec, Y), axis=1) / np.maximum(
        np.linalg.norm(rec, axis=1) * np.linalg.norm(Y, axis=1), 1e-30)
    return np.arcsin(np.clip(cr, -1, 1)).max() * 206265.0, C


def pack_body(jds, X, dt, N, max_days, target):
    n = len(X)
    maxM = max(N + 4, int(round(max_days / dt)))
    minM = N + 4
    segs = []  # (i, j, C)
    i = 0
    while i < n - 1:
        if i + minM - 1 > n - 1:
            # leftover tail: extend the previous segment to the end
            if segs:
                pi, _, _ = segs.pop()
                _, C = seg_fit(X, pi, n - 1, N)
                segs.append((pi, n - 1, C))
            break
        lo, hi = i + minM - 1, min(n - 1, i + maxM - 1)
        best, bestC = None, None
        while lo <= hi:
            mid = (lo + hi) // 2
            e, C = seg_fit(X, i, mid, N)
            if e <= target:
                best, bestC = mid, C
                lo = mid + 1
            else:
                hi = mid - 1
        if best is None:
            best, bestC = i + minM - 1, seg_fit(X, i, i + minM - 1, N)[1]
        segs.append((i, best, bestC))
        i = best
    return segs


def pack_body_q(segs, jds, X, N, bits):
    ncoef = N + 1
    nseg = len(segs)
    tstart = np.empty(nseg + 1)
    scale = np.empty((nseg, 3), np.float32)
    qs = np.empty((nseg, 3, ncoef), np.int32)
    maxerr = 0.0
    for k, (i, j, C) in enumerate(segs):
        tstart[k] = jds[i]
        sc = np.abs(C).max(axis=0) / IQ
        sc = np.maximum(sc, 1e-300)
        q = np.round(C / sc).astype(np.int32)
        scale[k] = sc
        qs[k] = q.T  # [coord][coef]
        # verify
        M = j - i + 1
        uu = 2.0 * np.arange(M) / (M - 1) - 1.0
        V = np.polynomial.chebyshev.chebvander(uu, N)
        rec = V @ (q * sc)
        Y = X[i:j + 1]
        cr = np.linalg.norm(np.cross(rec, Y), axis=1) / np.maximum(
            np.linalg.norm(rec, axis=1) * np.linalg.norm(Y, axis=1), 1e-30)
        maxerr = max(maxerr, np.arcsin(np.clip(cr, -1, 1)).max() * 206265.0)
    tstart[nseg] = jds[segs[-1][1]]
    return tstart, scale, qs, maxerr


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", type=float, default=0.01)
    ap.add_argument("--bits", type=int, default=32, choices=(24, 32))
    args = ap.parse_args()
    global IQ
    IQ = float((1 << (args.bits - 1)) - 1)
    ts = load.timescale()
    jd0 = float(ts.utc(1850, 1, 1).tt)
    jd1 = float(ts.utc(2140, 1, 1).tt)
    samp0 = jd0 - 5.0
    samp1 = jd1 + 200.0

    chunks = []
    index = []
    offset = 0
    for body, key, frame, dt, N, maxd in BODIES:
        jds, X = sample(key, frame, samp0, samp1, dt)
        segs = pack_body(jds, X, dt, N, maxd, args.target)
        tstart, scale, qs, err = pack_body_q(segs, jds, X, N, args.bits)
        ncoef = N + 1
        nseg = len(segs)
        blob = struct.pack(f"<{len(tstart)}d", *tstart)
        blob += scale.astype('<f4').tobytes()
        if args.bits == 32:
            blob += qs.astype('<i4').tobytes()
        else:
            u = (qs & 0xFFFFFF).astype(np.uint32)
            blob += u.view(np.uint8).reshape(-1, 4)[:, :3].tobytes()
        index.append((body, nseg, ncoef, offset, len(tstart)))
        chunks.append(blob)
        offset += len(blob)
        sys.stderr.write(f"  body {body:2d}: {nseg:5d} seg, N={N:2d}, "
                         f"err={err:.4f}\", {len(blob)/1e6:.3f} MB\n")

    version = 1 if args.bits == 32 else 2
    header = b"EPHP" + struct.pack("<IIddI", version, len(BODIES), jd0, jd1,
                                   len(index))
    data_base = 36 + 40 * len(index)
    idx = b""
    for (body, nseg, ncoef, off, ntstart) in index:
        off += data_base
        so = off + 8 * ntstart
        co = so + 4 * nseg * 3
        idx += struct.pack("<IIIIQQQ", body, nseg, ncoef, ntstart, off, so, co)
    header += struct.pack("<I", len(idx)) + idx
    with open(OUT, "wb") as f:
        f.write(header)
        for c in chunks:
            f.write(c)
    sys.stderr.write(f"wrote {OUT} ({os.path.getsize(OUT)/1e6:.2f} MB)\n")


if __name__ == "__main__":
    main()
