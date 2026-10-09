#!/usr/bin/env python3
"""Print reference values for tests/test_search.cpp and tests/test_time.cpp
(development only) from JPL DE440 with Skyfield:

* the first time (TT) after a date when the Sun's and the Moon's apparent
  ecliptic longitude of date (true equinox) reaches a given value;
* the equation of time at a UT date: apparent minus mean solar time, i.e.
  GAST - apparent RA of the Sun + 12 h - UT of day, wrapped to +-12 h.

    python3 tools/make_search_reference.py [path/to/de440.bsp]
"""
import math
import os
import sys

from skyfield.api import load

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
# Start (TT) and target longitude (degrees) for the crossing searches.
CROSSINGS = [
    (2461756.230099, 103.254195), (2455078.703751, 279.912270),
    (2457901.865740, 329.536843), (2461509.841230, 330.565546),
    (2440664.693362, 157.409280), (2452123.608174, 23.455511),
    (2440140.651267, 299.023717), (2464582.556107, 282.459547),
]
EOT_DATES = [2448908.50, 2451545.00, 2460400.25]  # UT


def wrap(d):
    return (d + 180.0) % 360.0 - 180.0


def crossing(eph, ts, body, start, target, step):
    earth = eph["earth"]

    def f(jd):
        lon = earth.at(ts.tt_jd(jd)).observe(eph[body]).apparent() \
            .ecliptic_latlon(epoch="date")[1].degrees
        return wrap(lon - target)

    a, fa = start, f(start)
    while True:
        b = a + step
        fb = f(b)
        if fa < 0 <= fb and fb - fa < 180:
            break
        a, fa = b, fb
    for _ in range(60):  # bisection to ~1e-10 day
        m = (a + b) / 2
        if f(m) < 0:
            a = m
        else:
            b = m
    return (a + b) / 2


def equation_of_time(eph, ts, jd_ut):
    t = ts.ut1_jd(jd_ut)
    ra = eph["earth"].at(t).observe(eph["sun"]).apparent().radec(epoch="date")[0].hours
    e = t.gast - ra + 12.0 - ((jd_ut + 0.5) % 1.0) * 24.0
    return ((e + 12.0) % 24.0 - 12.0) / 24.0


def main():
    eph = load(sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, "de440.bsp"))
    ts = load.timescale()
    print("crossings (start, longitude, Sun, Moon):")
    for jd, lon in CROSSINGS:
        print("    {%.6f, %.6f, %.9f, %.9f}," % (
            jd, lon, crossing(eph, ts, "sun", jd, lon, 5.0),
            crossing(eph, ts, "moon", jd, lon, 0.5)))
    print("equation of time (UT, days):")
    for jd in EOT_DATES:
        print("    {%.2f, %.12f}," % (jd, equation_of_time(eph, ts, jd)))


if __name__ == "__main__":
    main()
