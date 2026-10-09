#!/usr/bin/env python3
"""Build the fixed-star catalogue data/stars_data.{h,cpp} (development only).

Sources (public data):
  * names: the IAU Working Group on Star Names list (IAU-CSN, E. Mamajek)
  * astrometry: Hipparcos, the new reduction (van Leeuwen 2007; VizieR
    I/311), epoch J1991.25
  * radial velocities: SIMBAD
for the IAU-named stars with a Hipparcos number, plus the Galactic Center
(Sgr A*).  Positions are carried to epoch J2000.0 with ERFA's rigorous
space motion (pmsafe); src/stars.cpp continues from there.  `make data`
packs the result into data/eph-stars.bin.  Needs network access (or the
cached files in --cache) and pyerfa.

    python3 tools/make_stars.py [--cache DIR]
"""
from __future__ import annotations

import argparse
import csv
import io
import math
import os
import re
import urllib.parse
import urllib.request

import erfa

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CSN_URL = "https://www.pas.rochester.edu/~emamajek/WGSN/IAU-CSN.txt"
VIZIER_TAP = "https://tapvizier.cds.unistra.fr/TAPVizieR/tap/sync"
SIMBAD_TAP = "https://simbad.cds.unistra.fr/simbad/sim-tap/sync"
GREEK = {"alf": "alpha", "bet": "beta", "gam": "gamma", "del": "delta",
         "eps": "epsilon", "zet": "zeta", "eta": "eta", "tet": "theta",
         "iot": "iota", "kap": "kappa", "lam": "lambda", "mu": "mu",
         "nu": "nu", "ksi": "xi", "omi": "omicron", "pi": "pi", "rho": "rho",
         "sig": "sigma", "tau": "tau", "ups": "upsilon", "phi": "phi",
         "chi": "chi", "psi": "psi", "ome": "omega"}


def fetch(url, cache, name, data=None):
    path = os.path.join(cache, name)
    if not os.path.exists(path):
        with urllib.request.urlopen(url, data, timeout=180) as r:
            text = r.read().decode("utf-8", "replace")
        with open(path, "w", encoding="utf-8") as f:
            f.write(text)
    with open(path, encoding="utf-8") as f:
        return f.read()


def tap(url, query, cache, name):
    data = urllib.parse.urlencode({"REQUEST": "doQuery", "LANG": "ADQL",
                                   "FORMAT": "csv", "QUERY": query}).encode()
    return list(csv.DictReader(io.StringIO(fetch(url, cache, name, data))))


def names(text):
    """(name, designation, HIP, V) from the IAU-CSN list: fixed-width name
    columns, then fields read from the right (the designation column can
    hold one or two words, a notes column may follow the date)."""
    out = []
    for ln in text.splitlines():
        if not ln[:1].isalpha():
            continue
        name = ln[:18].strip()
        tok = ln[36:].split()
        while tok and not re.match(r"\d{4}-\d\d-\d\d$", tok[-1]):
            tok.pop()
        if len(tok) < 13 or tok[-5] == "_":
            continue
        hip, vmag, con, bayer = int(tok[-5]), float(tok[-7]), tok[-10], tok[-12]
        designation = "" if bayer == "_" else "%s %s" % (bayer, con)
        out.append((name, designation, hip, vmag))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cache", default=os.path.join(ROOT, "build", "stars-cache"))
    args = ap.parse_args()
    os.makedirs(args.cache, exist_ok=True)

    stars = names(fetch(CSN_URL, args.cache, "IAU-CSN.txt"))
    hips = ",".join(str(h) for _, _, h, _ in stars)
    astro = {int(r["HIP"]): r for r in tap(
        VIZIER_TAP, 'SELECT HIP, RArad, DErad, Plx, pmRA, pmDE FROM "I/311/hip2" '
        "WHERE HIP IN (%s)" % hips, args.cache, "hip2.csv")}
    rv = {int(r["id"].split()[1]): float(r["rvz_radvel"] or 0.0) for r in tap(
        SIMBAD_TAP, "SELECT ident.id, basic.rvz_radvel FROM basic JOIN ident ON "
        "ident.oidref = basic.oid WHERE ident.id IN (%s)"
        % ",".join("'HIP %d'" % h for _, _, h, _ in stars), args.cache, "rv.csv")}

    rows = []
    for name, designation, hip, vmag in sorted(stars):
        a = astro.get(hip)
        if a is None:
            continue
        ra, dec = math.radians(float(a["RArad"])), math.radians(float(a["DErad"]))
        plx = max(float(a["Plx"]), 0.0)
        pm_ra_cos = float(a["pmRA"])
        pm_dec = float(a["pmDE"])
        mas = math.radians(1.0 / 3.6e6)
        # J1991.25 (TT) -> J2000.0 (pmsafe takes dRA/dt, not dRA/dt cos dec).
        ra2, dec2, pmr2, pmd2, px2, rv2 = erfa.pmsafe(
            ra, dec, pm_ra_cos * mas / math.cos(dec), pm_dec * mas, plx / 1000.0,
            rv.get(hip, 0.0), 2448349.0625, 0.0, 2451545.0, 0.0)
        rows.append((name, designation, hip, vmag, float(ra2), float(dec2),
                     float(pmr2) * math.cos(dec2) / mas, float(pmd2) / mas,
                     float(px2) * 1000.0, float(rv2)))
    # Sgr A*, the Galactic Center (Reid & Brunthaler 2004/2020), for the
    # galactic ayanamsas, with its apparent motion (the Sun's orbit around
    # the Galaxy, -6.4 mas/yr along the galactic plane).
    rows.append(("Galactic Center", "Sgr A*", 0, 0.0,
                 math.radians(266.4168371), math.radians(-29.0078106),
                 -3.156, -5.585, 0.0, 0.0))

    note = "/* Generated by tools/make_stars.py -- do not edit. */\n"
    with open(os.path.join(ROOT, "data", "stars_data.h"), "w") as f:
        f.write(note + "#pragma once\n\nnamespace eph {\nnamespace stars {\n\n"
                "/* ICRS position at epoch J2000.0 (radians), proper motion\n"
                " * (mas/yr; RA as dRA/dt cos dec), parallax (mas), radial\n"
                " * velocity (km/s), V magnitude. */\n"
                "struct Star {\n    const char* name;\n    const char* designation;\n"
                "    unsigned hip;\n    float vmag;\n    double ra, dec;\n"
                "    float pm_ra, pm_dec, parallax, rv;\n};\n\n"
                "extern const Star catalogue[];\nextern const int count;\n\n"
                "}  // namespace stars\n}  // namespace eph\n")
    with open(os.path.join(ROOT, "data", "stars_data.cpp"), "w") as f:
        f.write(note + '#include "stars_data.h"\n\nnamespace eph {\nnamespace stars {\n\n'
                "const Star catalogue[] = {\n")
        for r in rows:
            f.write('    {"%s", "%s", %d, %.2ff, %.12f, %.12f, %.3ff, %.3ff, %.3ff, %.2ff},\n' % r)
        f.write("};\n\nconst int count = sizeof catalogue / sizeof catalogue[0];\n\n"
                "}  // namespace stars\n}  // namespace eph\n")
    print("wrote %d stars" % len(rows))


if __name__ == "__main__":
    main()
