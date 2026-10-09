#!/usr/bin/env python3
"""Compare data files against a JPL kernel on random dates.

Apparent geocentric ecliptic positions (and longitude speeds) of the
library, with each given set of data files, against the kernel via
Skyfield (IAU-2006/2000A ecliptic of date).  Times are TT, so Delta-T does
not enter.  Development tool; needs numpy + skyfield and a built
build/libephemeris.a.

    python3 tools/check_accuracy.py --kernel de440.bsp --start 1550 --end 2650 \\
        data/eph-1as.bin data/eph-full.bin "data/de440/de440-1900.bin data/de440/de440-2000.bin"
"""
from __future__ import annotations

import argparse
import os
import subprocess
import tempfile

import numpy as np
from skyfield.api import load
from skyfield.framelib import ecliptic_frame

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
BODIES = {0: "mercury barycenter", 1: "venus barycenter", 3: "mars barycenter",
          4: "jupiter barycenter", 5: "saturn barycenter",
          6: "uranus barycenter", 7: "neptune barycenter",
          8: "pluto barycenter", 9: "sun", 10: "moon"}
NAMES = {0: "Mercury", 1: "Venus", 3: "Mars", 4: "Jupiter", 5: "Saturn",
         6: "Uranus", 7: "Neptune", 8: "Pluto", 9: "Sun", 10: "Moon"}

HELPER = r"""
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ephemeris/ephemeris.h"
int main(int argc, char** argv) {
    for (int i = 1; i < argc; ++i)
        if (eph_load_file(argv[i])) { fprintf(stderr, "cannot load %s\n", argv[i]); return 1; }
    const double d = 57.29577951308232;
    double jd;
    while (scanf("%lf", &jd) == 1)
        for (int b = 0; b <= 10; ++b) {
            eph_position p;
            if (b == 2) continue;
            if (eph_position_apparent_tt((eph_body)b, jd, &p)) printf("%d nan nan nan\n", b);
            else printf("%d %.12f %.12f %.12f\n", b, p.longitude * d, p.latitude * d,
                        p.longitude_speed * d);
        }
    return 0;
}
"""


def reference(kernel, jd):
    eph = load(kernel)
    ts = load.timescale()
    earth = eph["earth"]
    t = ts.tt_jd(jd)
    ref = {}
    for b, key in BODIES.items():
        lat, lon, _ = earth.at(t).observe(eph[key]).apparent().frame_latlon(
            ecliptic_frame)
        # Speed from +-h; bodies behind the Sun's disk are excluded below.
        h = 0.001 if b == 10 else 0.01
        l2 = [earth.at(ts.tt_jd(jd + s)).observe(eph[key]).apparent()
              .frame_latlon(ecliptic_frame)[1].degrees for s in (-h, h)]
        speed = ((l2[1] - l2[0] + 180) % 360 - 180) / (2 * h)
        ref[b] = (lon.degrees, lat.degrees, speed)
    sun = ref[9][0]
    for b in BODIES:
        elong = np.abs((ref[b][0] - sun + 180) % 360 - 180)
        ref[b] = ref[b] + (elong,)
    return ref


def run(exe, files, jd):
    out = subprocess.run([exe] + files, input="\n".join("%.10f" % x for x in jd),
                         capture_output=True, text=True, check=True).stdout
    r = {}
    for ln in out.split("\n"):
        if ln:
            f = ln.split()
            r.setdefault(int(f[0]), []).append([float(x) for x in f[1:]])
    return {b: np.array(v) for b, v in r.items()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--kernel", default=os.path.join(ROOT, "de440.bsp"))
    ap.add_argument("--start", type=float, default=1550)
    ap.add_argument("--end", type=float, default=2650)
    ap.add_argument("-n", type=int, default=3000, help="number of dates")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("sets", nargs="+",
                    help="data files to load; quote several for one set")
    args = ap.parse_args()

    ts = load.timescale()
    jd0 = float(ts.tt(args.start, 1, 1).tt) + 1
    jd1 = float(ts.tt(args.end, 1, 1).tt) - 1
    jd = np.sort(np.random.default_rng(args.seed).uniform(jd0, jd1, args.n))
    ref = reference(args.kernel, jd)

    with tempfile.TemporaryDirectory() as tmp:
        src, exe = os.path.join(tmp, "h.c"), os.path.join(tmp, "h")
        with open(src, "w") as f:
            f.write(HELPER)
        subprocess.run(["cc", "-O2", "-I" + os.path.join(ROOT, "include"), src,
                        os.path.join(ROOT, "build", "libephemeris.a"), "-lm",
                        "-o", exe], check=True)
        results = [(s, run(exe, s.split(), jd)) for s in args.sets]

    print("Apparent position vs %s, %d dates %g-%g (TT): max / median arcsec;"
          " speed: max arcsec/day (outside the Sun's disk)"
          % (os.path.basename(args.kernel), args.n, args.start, args.end))
    for name, r in results:
        print("\n" + name)
        for b in BODIES:
            a = r[b]
            ok = ~np.isnan(a[:, 0])
            lon, lat, spd, elong = (x[ok] for x in ref[b])
            if not ok.any():
                print("  %-8s no data" % NAMES[b])
                continue
            dl = ((a[ok, 0] - lon + 180) % 360 - 180) * np.cos(np.radians(lat))
            e = np.hypot(dl, a[ok, 1] - lat) * 3600
            vis = (elong > 0.3) | (b == 9)
            se = np.abs(a[ok, 2] - spd)[vis].max() * 3600
            worst = jd[ok][np.argmax(e)]
            print("  %-8s %8.3f / %6.3f   speed %7.4f   worst %s%s"
                  % (NAMES[b], e.max(), np.median(e), se,
                     ts.tt_jd(worst).utc_strftime("%Y-%m"),
                     "" if ok.all() else "   (%d dates without data)"
                     % (~ok).sum()))


if __name__ == "__main__":
    main()
