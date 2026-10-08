#!/usr/bin/env python3
"""Pack geocentric positions from a JPL DE kernel into a compact Chebyshev
table (our own format, public-domain source).

Each body is represented as piecewise Chebyshev polynomials of the
geocentric ecliptic-J2000 rectangular coordinates (x, y, z) in AU.  Each
segment stores, per coordinate, a float32 scale and int32 coefficients
quantized to that scale -- giving ~31 bits of relative precision (better
than float32 at the same 4 bytes/coefficient).

The result is written to data/ephem.bin; it is loaded at runtime by
src/compressed.cpp.  Development tool: needs numpy + skyfield.

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
INT32 = 2147483647.0

# (eph_body id, DE440s key, frame, sample dt days, interval days, degree)
# frame: "helio" = body - Sun, "geo" = body - Earth.  Planets and Earth are
# stored heliocentric (so light-time and the Sun's geocentric place work
# correctly); the Moon is stored geocentric to avoid cancellation.
BODIES = [
    (2, "earth", "helio", 0.5, 32, 12),
    (10, "moon", "geo", 0.125, 16, 16),
    (0, "mercury barycenter", "helio", 0.25, 32, 16),
    (1, "venus barycenter", "helio", 0.25, 32, 14),
    (3, "mars barycenter", "helio", 0.5, 32, 14),
    (4, "jupiter barycenter", "helio", 1.0, 64, 18),
    (5, "saturn barycenter", "helio", 1.0, 64, 16),
    (6, "uranus barycenter", "helio", 1.0, 64, 14),
    (7, "neptune barycenter", "helio", 1.0, 128, 20),
    (8, "pluto barycenter", "helio", 1.0, 128, 20),
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
    # ICRS -> mean ecliptic J2000 (the frame precess_ecliptic_j2000_to_date
    # expects), using the IAU-1980 J2000 obliquity.
    eps0 = erfa.obl80(2451545.0, 0.0)
    x = r * np.cos(dec) * np.cos(ra)
    y = r * np.cos(dec) * np.sin(ra)
    z = r * np.sin(dec)
    X = np.stack([x,
                  y * math.cos(eps0) + z * math.sin(eps0),
                  -y * math.sin(eps0) + z * math.cos(eps0)], axis=1)
    return jds, X


def fit_body(jds, X, dt, L, N):
    M = int(round(L / dt))          # number of steps in the interval
    nseg = (len(X) - 1) // M
    # M+1 samples per segment (segments share endpoints), so the domain
    # [tstart[i], tstart[i+1]] is exactly the fitted interval.
    Y = np.stack([X[i * M:i * M + M + 1] for i in range(nseg)])
    uu = 2.0 * np.arange(M + 1) / M - 1.0
    V = np.polynomial.chebyshev.chebvander(uu, N)
    P = np.linalg.pinv(V)
    C = np.einsum('km,imc->ikc', P, Y)  # (nseg, N+1, 3) float64

    # Quantize: per segment, per coordinate scale (float32) + int32.
    scale = np.abs(C).max(axis=1) / INT32            # (nseg, 3)
    scale = np.maximum(scale, 1e-300)
    q = np.round(C / scale[:, None, :]).astype(np.int32)
    Cq = q.astype(np.float64) * scale[:, None, :]

    rec = np.einsum('mk,ikc->imc', V, Cq)
    cr = np.linalg.norm(np.cross(rec, Y), axis=2) / np.maximum(
        np.linalg.norm(rec, axis=2) * np.linalg.norm(Y, axis=2), 1e-30)
    err = np.arcsin(np.clip(cr, -1, 1)).max() * 206265.0

    tstart = jds[::M][:nseg + 1].astype(np.float64)
    return tstart, scale.astype(np.float32), q, err, nseg


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", type=float, default=0.01)
    args = ap.parse_args()

    ts = load.timescale()
    # Advertised range; sample well beyond it so (a) light-time retardation
    # stays inside and (b) every body's last full segment reaches jd1.
    jd0 = float(ts.utc(1850, 1, 1).tt)
    jd1 = float(ts.utc(2140, 1, 1).tt)
    samp0 = jd0 - 5.0
    samp1 = jd1 + 200.0

    chunks = []
    index = []
    offset = 0
    total_err = 0.0
    for body, key, frame, dt, L, N in BODIES:
        jds, X = sample(key, frame, samp0, samp1, dt)
        tstart, scale, q, err, nint = fit_body(jds, X, dt, L, N)
        ncoef = N + 1
        # layout: tstart (double), scale (float32, nint*3), q (int32, nint*3*ncoef)
        blob = struct.pack(f"<{len(tstart)}d", *tstart)
        blob += scale.astype('<f4').tobytes()
        # layout [segment][coord][coef] to match the C loader
        blob += np.transpose(q, (0, 2, 1)).astype('<i4').tobytes()
        index.append((body, nint, ncoef, offset, len(tstart)))
        chunks.append(blob)
        offset += len(blob)
        total_err = max(total_err, err)
        sys.stderr.write(f"  body {body:2d}: {nint:5d} seg, N={N:2d}, "
                         f"err={err:.4f}\", {len(blob)/1e6:.3f} MB\n")

    header = b"EPHP" + struct.pack("<IIddI", 1, len(BODIES), jd0, jd1,
                                   len(index))
    # index entries: body, nseg, ncoef, ntstart, tstart_off, scale_off, coef_off
    # Data offsets are absolute (header is 36 bytes + 40 per index entry).
    data_base = 36 + 40 * len(index)
    idx = b""
    for (body, nint, ncoef, off, ntstart) in index:
        off += data_base
        so = off + 8 * ntstart
        co = so + 4 * nint * 3
        idx += struct.pack("<IIIIQQQ", body, nint, ncoef, ntstart, off, so, co)
    header += struct.pack("<I", len(idx)) + idx
    with open(OUT, "wb") as f:
        f.write(header)
        for c in chunks:
            f.write(c)
    sys.stderr.write(f"wrote {OUT} ({os.path.getsize(OUT)/1e6:.2f} MB), "
                     f"max err {total_err:.4f}\"\n")


if __name__ == "__main__":
    main()
