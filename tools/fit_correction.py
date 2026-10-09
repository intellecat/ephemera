#!/usr/bin/env python3
"""Fit the DE440 corrections for VSOP87D and ELP2000-82B.

VSOP87D and ELP2000-82B are fitted to DE200/LE200 and drift from DE440:
slowly (Uranus/Neptune by ~10" over 1550-2650) and, because their orbital
elements' rates differ slightly, periodically with growing amplitude (the
Moon by ~45", Mars by ~2" at the ends of 1550-2650).  This tool takes the
uncorrected series from the built library, samples their residual against
a JPL kernel in heliocentric (geocentric for the Moon) spherical
coordinates of the IAU-1976 ecliptic of date, and fits each coordinate as

    P0(u) + sum_k [ Ck(u) cos(theta_k) + Sk(u) sin(theta_k) ]

with Chebyshev polynomials in u (T, centuries from J2000, scaled to the
window) and theta_k the body's fundamental arguments: harmonics of a
planet's mean longitude, or combinations of the Moon's Delaunay arguments.
The correction is clamped to the window at runtime (src/correction.cpp).

Emits data/corr_data.{h,cpp}.  Needs numpy, skyfield and pyerfa, a kernel
covering the window (de440.bsp: 1550-2650), and `make && make tools`
(build/libephemeris.a and build/pack_series).  Run `make data` afterwards.

    python3 tools/fit_correction.py --kernel de440.bsp --start 1550 --end 2650
"""
from __future__ import annotations

import argparse
import os
import subprocess
import sys
import tempfile

import erfa
import numpy as np
from skyfield.api import load

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
J2000 = 2451545.0
AS = np.pi / 648000.0
DEG = np.pi / 180.0

# Mean longitudes, degrees and degrees/century (Standish 1992, J2000).
MEAN_LONGITUDE = {
    "MERCURY": (252.25032350, 149472.67411175),
    "VENUS": (181.97909950, 58517.81538729),
    "EARTH": (100.46457166, 35999.37244981),
    "MARS": (-4.55343205, 19140.30268499),
    "JUPITER": (34.39644051, 3034.74612775),
    "SATURN": (49.95424423, 1222.49362201),
    "URANUS": (313.23810451, 428.48202785),
    "NEPTUNE": (-55.12002969, 218.45945325),
}
# Delaunay arguments l, l', F, D (arcsec, arcsec/century; Simon et al. 1994).
DELAUNAY = [(485868.249036, 1717915923.2178), (1287104.79305, 129596581.0481),
            (335779.526232, 1739527262.8478), (1072260.70369, 1602961601.2090)]
# The Moon's arguments as multipliers of (l, l', F, D): the largest terms of
# the lunar longitude and latitude.
MOON_ARGS = [(1, 0, 0, 0), (0, 0, 1, 0), (0, 0, 0, 2), (-1, 0, 0, 2),
             (2, 0, 0, 0), (0, 1, 0, 0), (0, 0, 0, 1), (0, 0, 2, 0),
             (0, 0, -1, 2), (0, -1, 0, 2), (1, 0, 1, 0), (1, 0, -1, 0),
             (0, 0, 1, 2), (-1, 0, 1, 2), (-1, 0, -1, 2), (1, 0, 0, 2),
             (-2, 0, 0, 2), (1, -1, 0, 0), (1, 1, 0, 0), (-1, -1, 0, 2)]

# (C++ tag, body index in the probe output, kernel key, centre, harmonics)
BODIES = [
    ("EARTH", 2, "earth", "sun", 2),
    ("MERCURY", 0, "mercury barycenter", "sun", 3),
    ("VENUS", 1, "venus barycenter", "sun", 2),
    ("MARS", 3, "mars barycenter", "sun", 2),
    ("JUPITER", 4, "jupiter barycenter", "sun", 2),
    ("SATURN", 5, "saturn barycenter", "sun", 2),
    ("URANUS", 6, "uranus barycenter", "sun", 2),
    ("NEPTUNE", 7, "neptune barycenter", "sun", 2),
    ("MOON", 10, "moon", "earth", len(MOON_ARGS)),
]
SECULAR = [20, 40, 60, 80, 100, 120]  # candidate secular degrees (+1 terms)
NH = 8                                # terms per harmonic amplitude

# Uncorrected series, as src/vsop87.cpp and src/moon.cpp evaluate them.
PROBE = r"""
#include <cstdio>
#include "ephemeris/ephemeris.h"
#include "internal.h"
using namespace eph;
int main(int argc, char** argv) {
    if (argc != 2 || eph_load_file(argv[1]) != 0) return 1;
    double jd;
    while (std::scanf("%lf", &jd) == 1) {
        for (int p = 0; p < 8; ++p) {
            double l, b, r;
            vsop87_base(p, jd, l, b, r);
            std::printf("%.17g %.17g %.17g ", l, b, r);
        }
        double x, y, z;
        elp_base_j2000(jd, x, y, z);
        double l, b, r;
        cartesian_to_spherical({x, y, z}, l, b, r);
        precess_ecliptic_j2000_to_date(jd, l, b);
        std::printf("%.17g %.17g %.17g\n", l, b, r / AU_KM);
    }
    return 0;
}
"""


def library(jds, tmp):
    """Uncorrected series: the full series packed as a table base (the base
    role carries no correction)."""
    base = os.path.join(tmp, "base.bin")
    subprocess.run([os.path.join(ROOT, "build", "pack_series"), "base", base,
                    "0", "0"], check=True, stderr=subprocess.DEVNULL)
    src, exe = os.path.join(tmp, "probe.cpp"), os.path.join(tmp, "probe")
    with open(src, "w") as f:
        f.write(PROBE)
    subprocess.run(["c++", "-O2", "-std=c++17",
                    "-I" + os.path.join(ROOT, "include"),
                    "-I" + os.path.join(ROOT, "src"), src,
                    os.path.join(ROOT, "build", "libephemeris.a"), "-o", exe],
                   check=True)
    out = subprocess.run([exe, base], input="\n".join("%.9f" % x for x in jds),
                         capture_output=True, text=True, check=True).stdout
    a = np.array(out.split(), float).reshape(len(jds), 9, 3)
    return {i if i < 8 else 10: a[:, i] for i in range(9)}


def kernel_spherical(eph, ts, key, center, jds):
    """Mean ecliptic/equinox of date (IAU 1976) spherical coordinates."""
    t = ts.tt_jd(jds)
    v = (eph[key].at(t) - eph[center].at(t)).position.au.T  # ICRS
    dt = jds - J2000
    v = np.einsum("nij,nj->ni", erfa.pmat76(J2000, dt), v)
    eps = erfa.obl80(J2000, dt)
    c, s = np.cos(eps), np.sin(eps)
    y = v[:, 1] * c + v[:, 2] * s
    z = -v[:, 1] * s + v[:, 2] * c
    r = np.linalg.norm(v, axis=1)
    return np.stack([np.arctan2(y, v[:, 0]), np.arcsin(z / r), r], axis=1)


def arguments(tag, k):
    """[(phase, rate)] in radians and radians/century."""
    if tag == "MOON":
        return [(sum(m * a for m, (a, _) in zip(mult, DELAUNAY)) * AS,
                 sum(m * b for m, (_, b) in zip(mult, DELAUNAY)) * AS)
                for mult in MOON_ARGS[:k]]
    a, b = MEAN_LONGITUDE[tag]
    return [(n * a * DEG, n * b * DEG) for n in range(1, k + 1)]


def design(u, T, n0, harm):
    cols = [np.polynomial.chebyshev.chebvander(u, n0 - 1)]
    vh = np.polynomial.chebyshev.chebvander(u, NH - 1)
    for phase, rate in harm:
        th = phase + rate * T
        cols += [vh * np.cos(th)[:, None], vh * np.sin(th)[:, None]]
    return np.hstack(cols)


def solve(X, y):
    """Least squares by normal equations in chunks (the columns are nearly
    orthogonal, so this is well conditioned)."""
    A = np.zeros((X.shape[1], X.shape[1]))
    b = np.zeros(X.shape[1])
    for i in range(0, len(y), 50000):
        A += X[i:i + 50000].T @ X[i:i + 50000]
        b += X[i:i + 50000].T @ y[i:i + 50000]
    return np.linalg.solve(A, b)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--kernel", default=os.path.join(ROOT, "de440.bsp"))
    ap.add_argument("--start", type=int, default=1550)
    ap.add_argument("--end", type=int, default=2650)
    args = ap.parse_args()

    eph = load(args.kernel)
    ts = load.timescale()
    jds = np.arange(float(ts.tt(args.start, 1, 1).tt),
                    float(ts.tt(args.end, 1, 1).tt), 1.0)
    T = (jds - J2000) / 36525.0
    tmin, tmax = float(T[0]), float(T[-1])
    u = (2.0 * T - (tmin + tmax)) / (tmax - tmin)
    with tempfile.TemporaryDirectory() as tmp:
        lib = library(jds, tmp)

    models = {}
    for tag, idx, key, center, nharm in BODIES:
        step = 1 if tag == "MOON" else 4  # samples (days)
        tr = np.arange(0, len(T), 2 * step)
        te = np.arange(step, len(T), 2 * step)
        d = kernel_spherical(eph, ts, key, center, jds) - lib[idx]
        d[:, 0] = (d[:, 0] + np.pi) % (2 * np.pi) - np.pi
        harm = arguments(tag, nharm)
        best = None
        for n0 in SECULAR:
            X, Xt = design(u[tr], T[tr], n0, harm), design(u[te], T[te], n0, harm)
            coefs = [solve(X, d[tr, k]) for k in range(3)]
            err = [np.abs(Xt @ c - d[te, k]).max() for k, c in enumerate(coefs)]
            score = max(err[0], err[1])
            if best is None or score < 0.9 * best[0] - 0.001 * AS:
                best = (score, n0, coefs, err)
        _, n0, coefs, err = best
        models[tag] = (n0, harm, coefs)
        raw = np.abs(d).max(axis=0)
        print("%-8s n0=%3d harmonics=%2d  lon %7.3f -> %6.4f\"  lat %7.3f -> "
              "%6.4f\"  rad %.1e -> %.1e AU"
              % (tag, n0, len(harm), raw[0] / AS, err[0] / AS, raw[1] / AS,
                 err[1] / AS, raw[2], err[2]), flush=True)
    emit(models, tmin, tmax, args)


def emit(models, tmin, tmax, args):
    note = ("/* Generated by tools/fit_correction.py (%s, %d-%d) -- do not "
            "edit. */\n" % (os.path.basename(args.kernel), args.start, args.end))
    h = os.path.join(ROOT, "data", "corr_data.h")
    cpp = os.path.join(ROOT, "data", "corr_data.cpp")
    with open(h, "w") as f:
        f.write(note + "#pragma once\n\nnamespace eph {\nnamespace corr {\n\n")
        f.write("/* theta = phase + rate * T (radians; T in centuries from J2000). */\n"
                "struct Harmonic {\n    double phase, rate;\n};\n\n"
                "/* Per coordinate (longitude and latitude in radians, distance in\n"
                " * AU): n0 Chebyshev coefficients of the secular part, then for each\n"
                " * harmonic nh of its cosine and nh of its sine amplitude, all in\n"
                " * u = (2T - tmin - tmax) / (tmax - tmin). */\n"
                "struct Series {\n    double tmin, tmax;\n    int n0, nh, nharm;\n"
                "    const Harmonic* harm;\n    const float* coef; "
                "/* 3 * (n0 + 2 * nharm * nh) */\n};\n\n")
        f.write("extern const Series earth, mercury, venus, mars, jupiter,\n"
                "                    saturn, uranus, neptune, moon;\n\n")
        f.write("}  // namespace corr\n}  // namespace eph\n")
    with open(cpp, "w") as f:
        f.write(note + '#include "corr_data.h"\n\nnamespace eph {\nnamespace corr {\n\n')
        for tag, (n0, harm, coefs) in models.items():
            f.write("static const Harmonic %s_harm[%d] = {\n" % (tag, len(harm)))
            for phase, rate in harm:
                f.write("    {%.17g, %.17g},\n" % (phase, rate))
            f.write("};\n\n")
            flat = np.concatenate(coefs).astype(np.float32)
            f.write("static const float %s_coef[%d] = {\n" % (tag, len(flat)))
            for i in range(0, len(flat), 6):
                f.write("    " + ", ".join("%.9gf" % x for x in flat[i:i + 6]) + ",\n")
            f.write("};\n\n")
        for tag, (n0, harm, coefs) in models.items():
            f.write("const Series %s = {%.17g, %.17g, %d, %d, %d, %s_harm, %s_coef};\n"
                    % (tag.lower(), tmin, tmax, n0, NH, len(harm), tag, tag))
        f.write("\n}  // namespace corr\n}  // namespace eph\n")
    sys.stderr.write("wrote %s\nwrote %s\n" % (h, cpp))


if __name__ == "__main__":
    main()
