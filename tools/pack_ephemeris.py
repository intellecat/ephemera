#!/usr/bin/env python3
"""Pack the DE440 tables (data/de440/): DE440 minus the library's series.

For each body the reference ephemeris (JPL DE440 via Skyfield) is compared
with a truncated copy of the library's analytical series (VSOP87D,
ELP2000-82B), or with nothing for bodies where the series would cost more
time than it saves space.  The truncated series are stored in the same
file (written by build/pack_series), so the table needs no other data.  The residual is small and
smooth, so it is stored as Chebyshev polynomials over fixed-length
segments with one absolute quantization step per body and a fixed bit
width per coefficient.  For every body the segment length, degree and step
giving the smallest file within --target (arcsec, as seen from the Earth at
the body's closest approach) are chosen.  The format is described in
src/data.h and src/compressed.cpp.

The residuals come from tools/residual.cpp, built here against
build/libephemeris.a (run `make` and `make tools` first), so they match the
runtime exactly.
A table must be regenerated whenever the series change.

Development tool; needs numpy, pyerfa and skyfield.

    python3 tools/pack_ephemeris.py                 # data/de440/de440-YEAR.bin
    python3 tools/pack_ephemeris.py --chunk 0 --start 1850 --end 2140 \
        --out data/ephem.bin                        # one file
"""
from __future__ import annotations

import argparse
import multiprocessing
import os
import struct
import subprocess
import sys
import tempfile
import zlib

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
LENGTHS = [2.0 ** k for k in range(2, 12)]
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


def residual(exe, body, base_file, jds, X):
    args = [exe, str(body)] + ([base_file] if base_file else [])
    inp = "\n".join("%.10f %.17g %.17g %.17g" % (j, *x) for j, x in zip(jds, X))
    out = subprocess.run(args, input=inp, capture_output=True, text=True,
                         check=True).stdout
    return np.array(out.split(), float).reshape(-1, 3)


def fit(R, t0, dt, lo, hi, quantum_au, speed_au):
    """Best (bytes, L, N, quantum, Q, err_au, start) over segment lengths,
    degrees and quantization steps for residual samples R (at t0 + i dt)
    covering lo..hi (on the sample grid), or None.  The segments start at
    lo, or end at the last sample if that is too close.  The position error
    must stay within quantum_au and the rounding error of the derivative
    (velocities come from it; it grows as N^2) within speed_au per day."""
    best = None
    for L in LENGTHS:
        m = int(round(L / dt))
        if m < 8 or abs(m * dt - L) > 1e-9:
            continue
        nseg = int(np.ceil((hi - lo) / L))
        i0 = int(round((lo - t0) / dt))
        if i0 + nseg * m > len(R) - 1:
            i0 = len(R) - 1 - nseg * m
        if i0 < 0:
            continue
        start = t0 + i0 * dt
        Y = R[i0 + np.arange(nseg)[:, None] * m + np.arange(m + 1)[None, :]]
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
                    best = (nbytes, L, N, q, Q, err, start)
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


def pack_range(args, year0, year1, exe, base_file, base_sections):
    """One data file: the base series and the residual table for
    year0..year1 (1 January, TT)."""
    eph = load(args.kernel)
    ts = load.timescale()
    k0 = max(sg.spk_segment.start_jd for sg in eph.segments)
    k1 = min(sg.spk_segment.end_jd for sg in eph.segments)
    jd0 = float(ts.tt(year0, 1, 1).tt)
    jd1 = float(ts.tt(year1, 1, 1).tt)
    records, blobs, log = [], [], []
    for body, key, center, dmin, dt, thr in BODIES:
        # Segments cover jd0 - 0.5 .. jd1 + 0.5; sample as far as the
        # longest segment could reach, within the kernel.
        lo, hi = jd0 - 0.5, jd1 + 0.5
        t0 = lo - dt * np.floor((lo - max(lo - max(LENGTHS), k0 + 0.5)) / dt)
        n = int(np.floor((min(hi + max(LENGTHS), k1 - 0.5) - t0) / dt)) + 1
        jds = t0 + dt * np.arange(n)
        R = residual(exe, body, base_file if thr is not None else None,
                     jds, reference(eph, ts, key, center, jds))
        quantum = args.target * ARCSEC * dmin
        nbytes, L, N, q, Q, err, start = fit(R, t0, dt, lo, hi, quantum,
                                             args.speed_target * ARCSEC * dmin)
        bits = widths(Q)
        records.append((body, thr, Q.shape[1], N + 1, start, L, q,
                        int(bits.sum()), bits))
        blobs.append(pack_bits(Q, bits))
        log.append("  %d body %2d: L=%5gd N=%2d %5d segs %7.1f KB  max %.4f\""
                   % (year0, body, L, N, Q.shape[1], len(blobs[-1]) / 1e3,
                      err / dmin / ARCSEC))

    base_crc = 0
    for t, _, payload in base_sections:
        if t in (1, 2):  # base VSOP and ELP sections, in file order
            base_crc = zlib.crc32(payload, base_crc)
    header = b"EPHP" + struct.pack("<IIddII", 4, len(records), jd0, jd1,
                                   base_crc, 0)
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
    table = header + index + b"".join(blobs)
    info = ("DE440 residual table %d-%d, %g\" (base: %s)"
            % (year0, year1, args.target,
               next(p for t, _, p in base_sections if t == 7).decode())).encode()
    sections = [(7, 0, info)] + [s for s in base_sections if s[0] != 7]
    sections.append((6, 0, table))
    return write_sections(sections), "\n".join(log)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", type=float, default=0.002,
                    help="max error per body, arcsec (default 0.002)")
    ap.add_argument("--speed-target", type=float, default=0.01,
                    help="max rounding error of the speeds, arcsec/day")
    ap.add_argument("--kernel", default=os.path.join(ROOT, "de440.bsp"))
    ap.add_argument("--start", type=int, default=1550)
    ap.add_argument("--end", type=int, default=2650)
    ap.add_argument("--chunk", type=int, default=100,
                    help="years per file (0: one file)")
    ap.add_argument("--out", default=os.path.join(ROOT, "data", "de440"),
                    help="directory for the de440-YEAR.bin chunks, or the "
                         "file with --chunk 0")
    ap.add_argument("--jobs", type=int, default=4)
    args = ap.parse_args()

    # Chunks on multiples of --chunk years (1900-2000, 2000-2100, ...).
    if args.chunk:
        cuts = [y for y in range(args.start + 1, args.end)
                if y % args.chunk == 0]
        bounds = [args.start] + cuts + [args.end]
    else:
        bounds = [args.start, args.end]
    ranges = list(zip(bounds[:-1], bounds[1:]))
    with tempfile.TemporaryDirectory() as tmp:
        exe = build_helper(tmp)
        # The base series: Mercury..Neptune thresholds (-1: none), ELP.
        thresholds = {b: t for b, _, _, _, _, t in BODIES}
        vsop = ",".join(repr(thresholds[p] if thresholds[p] is not None and
                             p != 2 else -1.0) for p in range(8))
        base_file = os.path.join(tmp, "base.bin")
        subprocess.run([os.path.join(ROOT, "build", "pack_series"), "base",
                        base_file, vsop,
                        repr(thresholds[10] if thresholds[10] is not None
                             else -1.0)], check=True)
        with open(base_file, "rb") as f:
            base_sections = read_sections(f.read())
        if args.chunk:
            os.makedirs(args.out, exist_ok=True)
        jobs = [(args, y0, y1, exe, base_file, base_sections) for y0, y1 in ranges]
        with multiprocessing.Pool(min(args.jobs, len(jobs))) as pool:
            for (y0, _), (data, log) in zip(ranges, pool.starmap(pack_range, jobs)):
                path = (os.path.join(args.out, "de440-%d.bin" % y0)
                        if args.chunk else args.out)
                with open(path, "wb") as f:
                    f.write(data)
                sys.stderr.write("%s\nwrote %s (%.0f KB)\n"
                                 % (log, path, len(data) / 1e3))


def read_sections(buf):
    """[(type, id, payload)] of a data file (format in src/data.h)."""
    assert buf[:4] == b"EPHD"
    _, n, _ = struct.unpack_from("<III", buf, 4)
    out = []
    for i in range(n):
        t, sid, off, size = struct.unpack_from("<IIQQ", buf, 16 + 24 * i)
        out.append((t, sid, buf[off:off + size]))
    return out


def write_sections(sections):
    """A data file; the header ends with the CRC-32 of the rest."""
    off = 16 + 24 * len(sections)
    index, body = b"", b""
    for t, sid, payload in sections:
        off_aligned = (off + len(body) + 7) // 8 * 8
        body += b"\0" * (off_aligned - off - len(body))
        index += struct.pack("<IIQQ", t, sid, off_aligned, len(payload))
        body += payload
    rest = index + body
    return (b"EPHD" + struct.pack("<III", 2, len(sections), zlib.crc32(rest))
            + rest)


if __name__ == "__main__":
    main()
