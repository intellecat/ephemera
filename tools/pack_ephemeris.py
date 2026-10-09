#!/usr/bin/env python3
"""Pack data/ephem.bin: DE440 minus the library's own series.

For each body the reference ephemeris (JPL DE440s via Skyfield) is compared
with a truncated copy of the library's analytical series (VSOP87D,
ELP2000-82B), or with nothing for bodies where the series would cost more
time than it saves space.  The residual is small and
smooth, so it is stored as Chebyshev polynomials over fixed-length
segments with one absolute quantization step per body and a fixed bit
width per coefficient.  For every body the segment length, degree and step
giving the smallest file within --target (arcsec, as seen from the Earth at
the body's closest approach) are chosen.  The format is described in
src/compressed.cpp.

The residuals come from tools/residual.cpp, built here against
build/libephemeris.a (run `make` first), so they match the runtime exactly.
A table must be regenerated whenever the series change.

Development tool; needs numpy, pyerfa and skyfield.

    /tmp/eph/bin/python tools/pack_ephemeris.py [--target 0.002]
"""
from __future__ import annotations

import argparse
import os
import struct
import subprocess
import sys
import tempfile

import erfa
import numpy as np
from skyfield.api import load

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
ARCSEC = np.pi / 648000.0

# (eph_body id, DE key, centre, closest geocentric distance (AU), sample
#  step (days), base series threshold in radians or None for no base).
# The thresholds trade file size for speed: the residual of the dropped
# terms costs table space.  For the slow outer planets a base saves little.
# The Earth's entry is the Earth-Moon barycentre (smooth, unlike the Earth);
# its error shows in every body, so it is held to Venus' closest distance.
BODIES = [
    (10, "moon", "earth", 0.00236, 0.25, 1e-7),
    (2, "earth barycenter", "sun", 0.264, 1.0, None),
    (0, "mercury barycenter", "sun", 0.549, 1.0, 1e-7),
    (1, "venus barycenter", "sun", 0.264, 1.0, 1e-7),
    (3, "mars barycenter", "sun", 0.372, 1.0, 1e-7),
    (4, "jupiter barycenter", "sun", 3.9, 4.0, None),
    (5, "saturn barycenter", "sun", 7.9, 4.0, None),
    (6, "uranus barycenter", "sun", 17.2, 8.0, None),
    (7, "neptune barycenter", "sun", 28.7, 8.0, None),
    (8, "pluto barycenter", "sun", 28.5, 8.0, None),
]
LENGTHS = [2.0 ** k for k in range(2, 12)]  # DE440s ends in 2150
MAX_COEF = 32


def build_helper(tmp):
    exe = os.path.join(tmp, "residual")
    subprocess.run(["c++", "-O2", "-std=c++17", "-I" + os.path.join(ROOT, "src"),
                    "-I" + os.path.join(ROOT, "include"),
                    os.path.join(HERE, "residual.cpp"),
                    os.path.join(ROOT, "build", "libephemeris.a"), "-o", exe],
                   check=True)
    return exe


def reference(eph, ts, key, center, jds):
    """Mean ecliptic J2000 (IAU-1976 obliquity) vectors, AU."""
    t = ts.tt_jd(jds)
    x = (eph[key].at(t) - eph[center].at(t)).position.au
    e = erfa.obl80(2451545.0, 0.0)
    c, s = np.cos(e), np.sin(e)
    return np.stack([x[0], c * x[1] + s * x[2], -s * x[1] + c * x[2]], axis=1)


def residual(exe, body, threshold, jd0, jd1, jds, X):
    base = threshold is not None
    args = [exe, str(body), "1" if base else "0", repr(threshold or 0.0),
            repr(jd0), repr(jd1)]
    inp = "\n".join("%.10f %.17g %.17g %.17g" % (j, *x) for j, x in zip(jds, X))
    out = subprocess.run(args, input=inp, capture_output=True, text=True,
                         check=True).stdout
    return np.array(out.split(), float).reshape(-1, 3)


def fit(R, dt, span, quantum_au, speed_au):
    """Best (bytes, L, N, quantum, Q, err_au) over segment lengths, degrees and
    quantization steps for residual samples R (step dt) covering `span`
    days, or None.  The position error must stay within quantum_au and the
    rounding error of the derivative (velocities come from it; it grows
    as N^2) within speed_au per day."""
    best = None
    for L in LENGTHS:
        m = int(round(L / dt))
        if m < 8 or abs(m * dt - L) > 1e-9:
            continue
        nseg = int(np.ceil(span / L))
        Y = R[np.arange(nseg)[:, None] * m + np.arange(m + 1)[None, :]]
        u = np.linspace(-1.0, 1.0, m + 1)
        for N in range(1, min(MAX_COEF, m - 2)):
            V = np.polynomial.chebyshev.chebvander(u, N)
            D = V @ np.array([np.pad(np.polynomial.chebyshev.chebder(e),
                                     (0, 1)) for e in np.eye(N + 1)]).T
            C = np.linalg.lstsq(V, Y.transpose(1, 0, 2).reshape(m + 1, -1),
                                rcond=None)[0].reshape(N + 1, nseg, 3)
            fit_err = np.linalg.norm(np.einsum("ij,jsc->sic", V, C) - Y,
                                     axis=-1).max()
            if fit_err > quantum_au:
                continue
            for f in (2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64):
                q = quantum_au / f
                Q = np.round(C / q)
                err = np.linalg.norm(np.einsum("ij,jsc->sic", V, Q * q) - Y,
                                     axis=-1).max()
                derr = np.linalg.norm(np.einsum("ij,jsc->sic", D, Q * q - C),
                                      axis=-1).max() * 2.0 / L
                if err > quantum_au or derr > speed_au:
                    continue
                bits = widths(Q)
                nbytes = nseg * bits.sum() / 8.0
                if best is None or nbytes < best[0]:
                    best = (nbytes, L, N, q, Q, err)
                break
    return best


def widths(Q):
    """Two's complement bit width per (coordinate, coefficient)."""
    a = np.abs(Q).max(axis=1).T  # (3, ncoef)
    b = np.zeros(a.shape, int)
    nz = a > 0
    b[nz] = np.floor(np.log2(a[nz])).astype(int) + 2
    return b


def pack_bits(Q, bits):
    """Segments of x, y, z coefficients as an LSB-first bit stream."""
    nseg = Q.shape[1]
    out = bytearray()
    acc, nacc = 0, 0
    for s in range(nseg):
        for c in range(3):
            for k in range(Q.shape[0]):
                w = int(bits[c][k])
                if w == 0:
                    continue
                acc |= (int(Q[k, s, c]) & ((1 << w) - 1)) << nacc
                nacc += w
                while nacc >= 8:
                    out.append(acc & 0xFF)
                    acc >>= 8
                    nacc -= 8
    if nacc:
        out.append(acc & 0xFF)
    return bytes(out) + b"\0" * 8


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", type=float, default=0.002,
                    help="max error per body, arcsec (default 0.002)")
    ap.add_argument("--speed-target", type=float, default=0.01,
                    help="max rounding error of the speeds, arcsec/day")
    ap.add_argument("--kernel", default=os.path.join(ROOT, "de440s.bsp"))
    ap.add_argument("--start", type=int, default=1850)
    ap.add_argument("--end", type=int, default=2140)
    ap.add_argument("--out", default=os.path.join(ROOT, "data", "ephem.bin"))
    args = ap.parse_args()

    eph = load(args.kernel)
    ts = load.timescale()
    jd0 = float(ts.utc(args.start, 1, 1).tt)
    jd1 = float(ts.utc(args.end, 1, 1).tt)

    records, blobs = [], []
    with tempfile.TemporaryDirectory() as tmp:
        exe = build_helper(tmp)
        for body, key, center, dmin, dt, thr in BODIES:
            # Two days of margin; sample far enough for the longest segment.
            start = jd0 - 2.0
            n = int(np.ceil((jd1 + 2.0 - start + max(LENGTHS)) / dt)) + 1
            jds = start + dt * np.arange(n)
            R = residual(exe, body, thr, jd0, jd1, jds,
                         reference(eph, ts, key, center, jds))
            quantum = args.target * ARCSEC * dmin
            nbytes, L, N, q, Q, err = fit(R, dt, jd1 + 2.0 - start, quantum,
                                          args.speed_target * ARCSEC * dmin)
            nseg = Q.shape[1]
            bits = widths(Q)
            blob = pack_bits(Q, bits)
            records.append((body, thr, nseg, N + 1, start, L, q, int(bits.sum()),
                            bits))
            blobs.append(blob)
            sys.stderr.write(
                "  body %2d: L=%5gd N=%2d %5d segs %7.1f KB  max %.4f\"\n"
                % (body, L, N, nseg, len(blob) / 1e3, err / dmin / ARCSEC))

    header = b"EPHP" + struct.pack("<IIdd", 3, len(records), jd0, jd1)
    offset = len(header) + 160 * len(records)
    index = b""
    for (body, thr, nseg, ncoef, start, L, q, segbits, bits), blob in zip(
            records, blobs):
        b = np.zeros((3, MAX_COEF), np.uint8)
        b[:, :ncoef] = bits
        index += struct.pack("<IIIIddddQII", body, thr is not None, nseg,
                             ncoef, thr or 0.0, start, L, q, offset, segbits, 0)
        index += b.tobytes()
        offset += len(blob)
    with open(args.out, "wb") as f:
        f.write(header + index)
        for blob in blobs:
            f.write(blob)
    sys.stderr.write("wrote %s (%.0f KB)\n" % (args.out,
                                               os.path.getsize(args.out) / 1e3))


if __name__ == "__main__":
    main()
