#!/usr/bin/env python3
"""Generate embedded C++ tables from the ELP2000-82B lunar theory.

Downloads the 36 IMCCE/CDS data files and folds the theory constants,
main-problem corrections and argument combinations into per-coordinate
Poisson series:

    C(t) = sum_k amp_k * t^tpow_k * sin( ph0 + ph1 t + ph2 t^2
                                          + ph3 t^3 + ph4 t^4 )

with t in Julian centuries from J2000.  The result is the Moon's
geocentric position in the mean dynamical ecliptic and inertial equinox of
J2000 (the coordinate change to rectangular form is done at runtime).

Development tool; requires network access on first run.  Emits
data/elp_data.{h,cpp}.

    python3 tools/gen_elp.py
"""
from __future__ import annotations

import argparse
import math
import os
import sys
import urllib.request

BASE = "https://ftp.imcce.fr/pub/ephem/moon/elp82b"
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
RAW = os.path.join(ROOT, "data", "elp")
OUT = os.path.join(ROOT, "data")

DEG = math.pi / 180.0
RAD = 648000.0 / math.pi  # arcseconds per radian


def download(name):
    os.makedirs(RAW, exist_ok=True)
    path = os.path.join(RAW, name)
    if not os.path.exists(path):
        sys.stderr.write(f"downloading {BASE}/{name}\n")
        with urllib.request.urlopen(f"{BASE}/{name}", timeout=120) as r, \
                open(path, "wb") as f:
            f.write(r.read())
    return path


def constants():
    c = {}
    c["am"] = 0.074801329518
    c["alpha"] = 0.002571881335
    c["dtasm"] = 2.0 * c["alpha"] / (3.0 * c["am"])

    w = [[0.0] * 5 for _ in range(3)]
    w[0][0] = (218 + 18 / 60 + 59.95571 / 3600) * DEG
    w[1][0] = (83 + 21 / 60 + 11.67475 / 3600) * DEG
    w[2][0] = (125 + 2 / 60 + 40.39816 / 3600) * DEG
    w[0][1] = 1732559343.73604 / RAD
    w[1][1] = 14643420.2632 / RAD
    w[2][1] = -6967919.3622 / RAD
    w[0][2] = -5.8883 / RAD
    w[1][2] = -38.2776 / RAD
    w[2][2] = 6.3622 / RAD
    w[0][3] = 0.6604e-2 / RAD
    w[1][3] = -0.45047e-1 / RAD
    w[2][3] = 0.7625e-2 / RAD
    w[0][4] = -0.3169e-4 / RAD
    w[1][4] = 0.21301e-3 / RAD
    w[2][4] = -0.3586e-4 / RAD
    c["w"] = w

    eart = [0.0] * 5
    eart[0] = (100 + 27 / 60 + 59.22059 / 3600) * DEG
    eart[1] = 129597742.2758 / RAD
    eart[2] = -0.0202 / RAD
    eart[3] = 0.9e-5 / RAD
    eart[4] = 0.15e-6 / RAD
    peri = [0.0] * 5
    peri[0] = (102 + 56 / 60 + 14.42753 / 3600) * DEG
    peri[1] = 1161.2283 / RAD
    peri[2] = 0.5327 / RAD
    peri[3] = -0.138e-3 / RAD
    peri[4] = 0.0
    c["eart"], c["peri"] = eart, peri

    c["precess"] = 5029.0966 / RAD

    p = [[0.0, 0.0] for _ in range(8)]
    p[0][0] = (252 + 15 / 60 + 3.25986 / 3600) * DEG
    p[1][0] = (181 + 58 / 60 + 47.28305 / 3600) * DEG
    p[2][0] = eart[0]
    p[3][0] = (355 + 25 / 60 + 59.78866 / 3600) * DEG
    p[4][0] = (34 + 21 / 60 + 5.34212 / 3600) * DEG
    p[5][0] = (50 + 4 / 60 + 38.89694 / 3600) * DEG
    p[6][0] = (314 + 3 / 60 + 18.01841 / 3600) * DEG
    p[7][0] = (304 + 20 / 60 + 55.19575 / 3600) * DEG
    p[0][1] = 538101628.68898 / RAD
    p[1][1] = 210664136.43355 / RAD
    p[2][1] = eart[1]
    p[3][1] = 68905077.59284 / RAD
    p[4][1] = 10925660.42861 / RAD
    p[5][1] = 4399609.65932 / RAD
    p[6][1] = 1542481.19393 / RAD
    p[7][1] = 786550.32074 / RAD
    c["p"] = p

    c["delnu"] = +0.55604 / RAD / w[0][1]
    c["dele"] = +0.01789 / RAD
    c["delg"] = -0.08066 / RAD
    c["delnp"] = -0.06424 / RAD / w[0][1]
    c["delep"] = -0.12879 / RAD

    del_ = [[0.0] * 5 for _ in range(4)]  # del_[i-1][k]
    for i in range(5):
        del_[0][i] = w[0][i] - eart[i]  # D
        del_[3][i] = w[0][i] - w[2][i]  # F
        del_[2][i] = w[0][i] - w[1][i]  # l'
        del_[1][i] = eart[i] - peri[i]  # l (Sun M)
    del_[0][0] += math.pi
    c["del"] = del_
    c["zeta"] = [w[0][0], w[0][1] + c["precess"]]
    return c


def parse_main(path):
    out = []
    with open(path) as f:
        f.readline()
        for line in f:
            if len(line) < 99:
                continue
            try:
                ilu = [int(line[3 * i:3 * i + 3]) for i in range(4)]
                coef = [float(line[14:27])]
                for j in range(6):
                    a = 29 + 12 * j
                    coef.append(float(line[a:a + 10]))
            except ValueError:
                continue
            out.append((ilu, coef))
    return out


def parse_pert(path, nint):
    out = []
    with open(path) as f:
        f.readline()
        for line in f:
            if len(line) < nint * 3 + 20:
                continue
            try:
                ints = [int(line[3 * i:3 * i + 3]) for i in range(nint)]
                pha = float(line[nint * 3 + 1:nint * 3 + 10])
                amp = float(line[nint * 3 + 11:nint * 3 + 20])
            except ValueError:
                continue
            out.append((ints, pha, amp))
    return out


def build(c, min_amp):
    """Return (main, pert): per coordinate, main-problem terms (5 phase
    coefficients) and perturbation terms (2 phase coefficients + t^tpow)."""
    main = [[], [], []]
    pert = [[], [], []]
    for ific in range(1, 37):
        itab = (ific + 2) // 3
        iv = (ific - 1) % 3 + 1  # 1=Lon,2=Lat,3=Dist
        name = f"ELP{ific}"
        path = download(name)

        if ific <= 3:
            rows = parse_main(path)
            for ilu, coef in rows:
                tgv = coef[1] + c["dtasm"] * coef[5]
                a1 = coef[0]
                if ific == 3:
                    a1 = a1 - 2.0 * a1 * c["delnu"] / 3.0
                xx = (a1 + tgv * (c["delnp"] - c["am"] * c["delnu"])
                      + coef[2] * c["delg"] + coef[3] * c["dele"]
                      + coef[4] * c["delep"])
                if abs(xx) < min_amp:
                    continue
                ph = [0.0] * 5
                for k in range(5):
                    ph[k] = sum(ilu[i] * c["del"][i][k] for i in range(4))
                if iv == 3:
                    ph[0] += math.pi / 2.0
                main[iv - 1].append((xx, ph))
        else:
            nint = 11 if 10 <= ific <= 21 else 5
            rows = parse_pert(path, nint)
            for ints, pha, amp in rows:
                if abs(amp) < min_amp:
                    continue
                ph = [0.0] * 5
                for k in range(2):
                    if nint == 5:
                        iz = ints[0]
                        ilu = ints[1:5]
                        y = iz * c["zeta"][k] + sum(
                            ilu[i] * c["del"][i][k] for i in range(4))
                    elif ific < 16:
                        ipla = ints
                        y = (ipla[8] * c["del"][0][k]
                             + ipla[9] * c["del"][2][k]
                             + ipla[10] * c["del"][3][k]
                             + sum(ipla[i] * c["p"][i][k] for i in range(8)))
                    else:
                        ipla = ints
                        y = (sum(ipla[i + 7] * c["del"][i][k] for i in range(4))
                             + sum(ipla[i] * c["p"][i][k] for i in range(7)))
                    if k == 0:
                        y += pha * DEG
                    ph[k] = y
                tpow = 1 if itab in (3, 5, 7, 9) else (2 if itab == 12 else 0)
                pert[iv - 1].append((amp, ph[0], ph[1], tpow))
    return main, pert


TWO_PI = 2.0 * math.pi
PHASE_BITS = 30
PHASE_SCALE = 1 << PHASE_BITS


def qph(x):
    """Quantize a phase (radians, reduced mod 2pi) to PHASE_BITS bits."""
    v = x % TWO_PI
    if v < 0.0:
        v += TWO_PI
    return int(round(v * PHASE_SCALE / TWO_PI)) & (PHASE_SCALE - 1)


def c_float(x):
    return f"{x:.17g}"


def emit(main, pert):
    h = os.path.join(OUT, "elp_data.h")
    cpp = os.path.join(OUT, "elp_data.cpp")
    with open(h, "w") as f:
        f.write("/* Generated by tools/gen_elp.py -- do not edit. */\n")
        f.write("#pragma once\n\nnamespace eph {\nnamespace elp {\n\n")
        f.write("struct MainTerm { double ph1, ph2, ph3, ph4; float amp;\n"
                "                  unsigned int packed; };\n")
        f.write("struct PertTerm { double ph1; float amp; unsigned int packed; };\n")
        f.write("struct Coord { const MainTerm* main; int nmain;\n"
                "               const PertTerm* pert; int npert; };\n\n")
        f.write("extern const Coord coords[3];\n\n")
        f.write("}  // namespace elp\n}  // namespace eph\n")
    with open(cpp, "w") as f:
        f.write("/* Generated by tools/gen_elp.py -- do not edit. */\n")
        f.write('#include "elp_data.h"\n\nnamespace eph {\nnamespace elp {\n\n')
        for i in range(3):
            f.write(f"static const MainTerm m{i}[] = {{\n")
            for amp, ph in main[i]:
                f.write("    {%s, %s, %s, %s, %s, %du},\n"
                        % (c_float(ph[1]), c_float(ph[2]), c_float(ph[3]),
                           c_float(ph[4]), c_float(amp), qph(ph[0])))
            f.write("};\n\n")
            f.write(f"static const PertTerm p{i}[] = {{\n")
            for amp, ph0, ph1, tpow in pert[i]:
                f.write("    {%s, %s, %du},\n"
                        % (c_float(ph1), c_float(amp),
                           (tpow << PHASE_BITS) | qph(ph0)))
            f.write("};\n\n")
        f.write("const Coord coords[3] = {\n")
        for i in range(3):
            f.write(f"    {{m{i}, {len(main[i])}, p{i}, {len(pert[i])}}},\n")
        f.write("};\n\n}  // namespace elp\n}  // namespace eph\n")
    return h, cpp


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--min-amp", type=float, default=1e-4,
                    help="Drop terms with |amplitude| below this (arcsec/km).")
    args = ap.parse_args()
    c = constants()
    main_, pert_ = build(c, args.min_amp)
    for i in range(3):
        sys.stderr.write(f"coord{i}: main={len(main_[i])} pert={len(pert_[i])}\n")
    sys.stderr.write("total %d terms\n" % sum(
        len(main_[i]) + len(pert_[i]) for i in range(3)))
    emit(main_, pert_)


if __name__ == "__main__":
    main()
