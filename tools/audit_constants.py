#!/usr/bin/env python3
"""Audit the hand-entered constants in src/ against their authoritative
sources (development only; needs network access, numpy, pyerfa, skyfield):

  * leap seconds (src/util.cpp)        IERS Leap_Second.dat
  * minor bodies' H, G, radius          JPL Small-Body Database API
    (src/phenomena.cpp)
  * planet / Sun mass ratios            DE440 constants (NAIF gm_de440.tpc)
    (src/orbits.cpp, src/stars.cpp)
  * Sun's barycentric velocity model    DE440 (Skyfield)
    (src/stars.cpp)
  * mean orbital elements at J2000      VSOP87 secular terms (IMCCE)
    (src/orbits.cpp)

Each line prints OK or DIFF with the two values.

    DE440=/path/to/de440.bsp python3 tools/audit_constants.py
"""
import json
import math
import os
import re
import urllib.request

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
problems = 0


def src(name):
    return open(os.path.join(ROOT, "src", name)).read()


def get(url):
    with urllib.request.urlopen(url, timeout=120) as r:
        return r.read().decode("utf-8", "replace")


def report(ok, what):
    global problems
    print(("OK    " if ok else "DIFF  ") + what)
    problems += not ok


def leap_seconds():
    text = get("https://hpiers.obspm.fr/iers/bul/bulc/Leap_Second.dat")
    iers = [(float(l.split()[0]) + 2400000.5, int(l.split()[4]))
            for l in text.splitlines() if re.match(r"\s+\d{5}\.0", l)]
    s = src("util.cpp")
    body = s[s.index("kLeapSeconds[] = {"):s.index("};", s.index("kLeapSeconds[] = {"))]
    ours = [(float(a), int(b)) for a, b in re.findall(r"\{\s*([0-9.]+)\s*,\s*(\d+)\s*\}", body)]
    report(ours == iers, "leap seconds: %d entries, last %s" % (len(ours), ours[-1]))


def minor_bodies():
    s = src("phenomena.cpp")
    for name, sstr in (("CERES", "1"), ("PALLAS", "2"), ("JUNO", "3"), ("VESTA", "4"),
                       ("CHIRON", "2060")):
        d = json.loads(get("https://ssd-api.jpl.nasa.gov/sbdb.api?sstr=%s&phys-par=true" % sstr))
        p = {x["name"]: x["value"] for x in d.get("phys_par", [])}
        h, g = re.search(r"case EPH_%s: return hg_magnitude\(([\d.]+), ([\d.]+)" % name, s).groups()
        r = re.search(r"case EPH_%s: return ([\d.]+);" % name, s).group(1)
        jh, jg = float(p["H"]), float(p["G"]) if p.get("G") else 0.15
        jr = float(p["diameter"]) / 2
        report(abs(float(h) - jh) <= 0.005 and abs(float(g) - jg) <= 0.005,
               "%s H, G: %s %s vs SBDB %s %s" % (name.title(), h, g, jh, jg))
        report(abs(float(r) / jr - 1) < 0.005, "%s radius: %s vs SBDB %.1f km" % (name.title(), r, jr))


def masses():
    text = get("https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/gm_de440.tpc")
    gm = {int(a): float(b.replace("D", "E"))
          for a, b in re.findall(r"BODY(\d+)_GM\s*=\s*\(\s*([0-9.D+E-]+)", text)}
    naif = {"MERCURY": 1, "VENUS": 2, "MARS": 4, "JUPITER": 5, "SATURN": 6, "URANUS": 7,
            "NEPTUNE": 8, "PLUTO": 9}
    s = src("orbits.cpp")
    for name, n in naif.items():
        ours = float(re.search(r"case EPH_%s: return 1\.0 / ([\d.]+);" % name, s).group(1))
        ref = gm[10] / gm[n]
        report(abs(ours / ref - 1) < 1e-6, "%s mass ratio: %s vs DE440 %.6f" % (name.title(), ours, ref))


def sun_velocity():
    import erfa
    from skyfield.api import load
    s = src("stars.cpp")
    body = s[s.index("kOrbits[4] = {"):s.index("};", s.index("kOrbits[4] = {"))]
    nums = []
    for tok in re.findall(r"1\.0 / [\d.]+|-?\d+\.\d+", body):
        nums.append(1.0 / float(tok.split("/")[1]) if "/" in tok else float(tok))
    orbits = [nums[i * 13:(i + 1) * 13] for i in range(4)]
    k, deg = 0.01720209895, math.pi / 180

    def vsun(jd):
        T = (jd - 2451545.0) / 36525
        v, total = np.zeros(3), 1.0
        for a, e, i, l, w, o, da, de, di, dl, dw, do, m in orbits:
            a, e, inc = a + da * T, e + de * T, (i + di * T) * deg
            w, node = (w + dw * T) * deg, (o + do * T) * deg
            M = math.fmod((l + dl * T) * deg - w, 2 * math.pi)
            E = M
            for _ in range(6):
                E = M + e * math.sin(E)
            ed = k / (a * math.sqrt(a)) / (1 - e * math.cos(E))
            xd, yd = -a * math.sin(E) * ed, a * math.sqrt(1 - e * e) * math.cos(E) * ed
            cw, sw, cn, sn = math.cos(w - node), math.sin(w - node), math.cos(node), math.sin(node)
            ci, si = math.cos(inc), math.sin(inc)
            v -= m * np.array([(cw * cn - sw * sn * ci) * xd - (sw * cn + cw * sn * ci) * yd,
                               (cw * sn + sw * cn * ci) * xd - (sw * sn - cw * cn * ci) * yd,
                               sw * si * xd + cw * si * yd])
            total += m
        return v / total

    eph = load(os.environ.get("DE440", os.path.join(ROOT, "de440.bsp")))
    ts = load.timescale()
    rot = erfa.ecm06(2451544.5, 0.5)
    worst = max(np.linalg.norm(vsun(jd) - rot @ eph["sun"].at(ts.tt_jd(jd)).velocity.au_per_d)
                for jd in np.linspace(2396759, 2506332, 300)) * 149597870700 / 86400
    report(worst < 0.5, "Sun barycentric velocity 1850-2150: worst %.2f m/s" % worst)


def mean_elements():
    s = src("orbits.cpp")
    rows = re.findall(r"/\* (\w+)(?:(?!\*/).)*\*/\s*\{\{([^}]*)\},\s*\{([^}]*)\},\s*"
                      r"\{([^}]*)\},\s*\{([^}]*)\},\s*\{([^}]*)\}\}", s, re.S)
    files = {"Mercury": "mer", "Venus": "ven", "Mars": "mar", "Jupiter": "jup", "Saturn": "sat",
             "Uranus": "ura", "Neptune": "nep"}
    found = [r[0] for r in rows]
    report(all(n in found for n in files), "mean element rows found: " + " ".join(found))
    for name, a, e, i, node, peri in rows:
        if name not in files:
            continue
        text = get("https://ftp.imcce.fr/pub/ephem/planets/vsop87/VSOP87." + files[name])
        const, var, power = {}, None, None
        for ln in text.splitlines():
            if "VARIABLE" in ln:
                var = int(re.search(r"VARIABLE (\d)", ln).group(1))
                power = int(re.search(r"\*T\*\*(\d)", ln).group(1))
                continue
            if power == 0 and all(int(ln[10 + 3 * j:13 + 3 * j]) == 0 for j in range(12)):
                const[var] = const.get(var, 0) + float(ln[79:97]) * math.cos(float(ln[97:111]))
        kk, hh, q, p = const[3], const[4], const[5], const[6]
        ref = (math.hypot(kk, hh), math.degrees(2 * math.asin(math.hypot(q, p))),
               math.degrees(math.atan2(p, q)) % 360, math.degrees(math.atan2(hh, kk)) % 360)
        ours = [float(x.split(",")[0]) for x in (e, i, node, peri)]
        d = [ours[0] - ref[0], ours[1] - ref[1], (ours[2] - ref[2] + 180) % 360 - 180,
             (ours[3] - ref[3] + 180) % 360 - 180]
        # Meeus' constants include long-period terms: up to 1e-4 in e,
        # 5e-4 deg in the perihelion.
        ok = abs(d[0]) < 1e-4 and abs(d[1]) < 1e-5 and abs(d[2]) < 1e-4 and abs(d[3]) < 6e-4
        report(ok, "%s mean elements at J2000 vs VSOP87: de %+.1e di %+.1e dnode %+.1e dperi %+.1e"
               % (name, *d))


def main():
    for check in (leap_seconds, minor_bodies, masses, sun_velocity, mean_elements):
        try:
            check()
        except Exception as exc:  # network or parse failure: report, continue
            report(False, "%s: could not check (%s)" % (check.__name__, exc))
    print("%d problem(s)" % problems)


if __name__ == "__main__":
    main()
