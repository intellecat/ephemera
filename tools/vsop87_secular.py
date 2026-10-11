#!/usr/bin/env python3
"""Check the mean orbital elements in src/orbits.cpp (development only):
the secular (non-periodic) terms of VSOP87 in elliptic elements (IMCCE
files VSOP87.mer ... VSOP87.nep, J2000 frame) give e, i, node and
perihelion; rotated to the IAU 2006 mean ecliptic of date (ERFA ecm06) and
fitted by cubics in Julian centuries from J2000 they reproduce Meeus'
Table 31.A (Simon et al. 1994) to its long-period terms and its IAU 1976
precession (~1e-4 deg/century in the node and perihelion rates).

    python3 tools/vsop87_secular.py DIR mer ven mar jup sat ura nep

DIR holds the files from https://ftp.imcce.fr/pub/ephem/planets/vsop87/.
"""
import sys, math, re
import numpy as np, erfa
D = sys.argv[1]
def secular(p):
    poly = {v: [0.0] * 6 for v in range(1, 7)}
    var = power = None
    for ln in open(f'{D}/VSOP87.{p}'):
        if 'VARIABLE' in ln:
            var = int(re.search(r'VARIABLE (\d)', ln).group(1)); power = int(re.search(r'\*T\*\*(\d)', ln).group(1)); continue
        if all(int(ln[10+3*j:13+3*j]) == 0 for j in range(12)):
            A = float(ln[79:97]); B = float(ln[97:111]); poly[var][power] += A * math.cos(B)
    return poly
eps0 = math.radians(84381.406 / 3600)
Rx = np.array([[1, 0, 0], [0, math.cos(eps0), -math.sin(eps0)], [0, math.sin(eps0), math.cos(eps0)]])  # ecl -> eq
for p in sys.argv[2:]:
    poly = secular(p)
    Ts = np.linspace(-10, 10, 81)  # centuries
    rows = []
    for T in Ts:
        tm = T / 10.0  # VSOP87 time: millennia
        val = {v: sum(c * tm**k for k, c in enumerate(poly[v])) for v in range(3, 7)}
        k, h, q, pp = val[3], val[4], val[5], val[6]
        e = math.hypot(k, h)
        # orbit normal and perihelion direction in the J2000 ecliptic
        sin_half = math.hypot(q, pp); i = 2 * math.asin(sin_half); node = math.atan2(pp, q)
        w = math.atan2(h, k) - node
        n = np.array([math.cos(node), math.sin(node), 0.0])
        hu = np.array([math.sin(i) * math.sin(node), -math.sin(i) * math.cos(node), math.cos(i)])
        qv = np.cross(hu, n); peri_dir = math.cos(w) * n + math.sin(w) * qv
        jd = 2451545.0 + T * 36525
        R = erfa.ecm06(2400000.5, jd - 2400000.5) @ Rx
        hu2, pd2 = R @ hu, R @ peri_dir
        i2 = math.acos(hu2[2]); node2 = math.atan2(hu2[0], -hu2[1])
        n2 = np.array([math.cos(node2), math.sin(node2), 0.0])
        w2 = math.atan2(np.dot(np.cross(n2, pd2), hu2), np.dot(n2, pd2))
        rows.append((e, math.degrees(i2), math.degrees(node2), math.degrees(node2 + w2)))
    rows = np.array(rows)
    for col, name in ((0, 'e'), (1, 'i'), (2, 'node'), (3, 'peri')):
        y = np.unwrap(np.radians(rows[:, col])) if col >= 2 else rows[:, col]
        if col >= 2: y = np.degrees(y)
        c = np.polynomial.polynomial.polyfit(Ts, y, 3)
        print(p, name, ' '.join('%.10g' % x for x in c))
