#!/usr/bin/env python3
"""Regenerate the test reference data (development only).

This uses Skyfield with a JPL DE440s kernel to produce geometric and
apparent geocentric positions for the Sun, Moon, Pluto and the eight
planets.  The DE kernel is NOT part of this library and is never shipped;
it is only used to produce tests/reference_data.h.

Skyfield expresses vectors in the ICRS/IAU-2000 frame, but VSOP87D (and
this library) use the mean ecliptic and equinox of date under the IAU-1976
(FK5) reduction.  The references are therefore converted from the ICRS
vector with ERFA's IAU-1976 precession matrix and IAU-1980 obliquity, so
the test measures the dynamics in the library's own frame rather than a
frame-model difference.

    python3 -m venv /tmp/eph && /tmp/eph/bin/pip install skyfield pyerfa
    /tmp/eph/bin/python tools/make_reference.py /tmp/ref.json
    python3 tools/gen_reference.py /tmp/ref.json

DE440s spans 1849-2150, hence the date range below.
"""
import json
import math
import sys

import erfa
import numpy as np
from skyfield.api import load

J2000 = 2451545.0


def icrs_to_ecliptic_date(ra, dec, jd_tt, apparent=False):
    """ICRS unit direction -> ecliptic/equinox of date (IAU 1976/1980).

    apparent=False: mean equinox of date (matches VSOP87D geometric).
    apparent=True : true equinox of date, i.e. includes nutation.
    """
    v = np.array([math.cos(dec) * math.cos(ra),
                  math.cos(dec) * math.sin(ra),
                  math.sin(dec)])
    dt = jd_tt - J2000
    if apparent:
        v = erfa.pnm80(J2000, dt) @ v
        _, deps = erfa.nut80(J2000, dt)
        eps = erfa.obl80(J2000, dt) + deps
    else:
        v = erfa.pmat76(J2000, dt) @ v
        eps = erfa.obl80(J2000, dt)
    y = v[1] * math.cos(eps) + v[2] * math.sin(eps)
    z = -v[1] * math.sin(eps) + v[2] * math.cos(eps)
    lon = math.degrees(math.atan2(y, v[0])) % 360.0
    lat = math.degrees(math.asin(z))
    return lon, lat


def main():
    out = sys.argv[1] if len(sys.argv) > 1 else "/tmp/ref.json"
    eph = load("de440s.bsp")
    ts = load.timescale()
    earth = eph["earth"]

    targets = {
        "Mercury": "mercury barycenter",
        "Venus": "venus barycenter",
        "Mars": "mars barycenter",
        "Jupiter": "jupiter barycenter",
        "Saturn": "saturn barycenter",
        "Uranus": "uranus barycenter",
        "Neptune": "neptune barycenter",
        "Pluto": "pluto barycenter",
        "Sun": "sun",
        "Moon": "moon",
    }
    dates = [
        (2000, 1, 1, 12.0), (1992, 4, 12, 0.0), (1900, 1, 1, 0.0),
        (2026, 1, 1, 0.0), (1950, 6, 15, 6.0), (2000, 6, 1, 0.0),
        (2050, 1, 1, 0.0), (1987, 4, 10, 0.0), (2020, 12, 21, 18.0),
        (1930, 3, 3, 3.0), (1850, 1, 1, 0.0), (1875, 1, 1, 0.0),
        (2100, 1, 1, 0.0), (2125, 1, 1, 0.0), (2140, 1, 1, 0.0),
    ]

    ref = {
        "dates": [{"utc": list(d), "tt_jd": ts.utc(*d).tt} for d in dates],
        "bodies": {},
    }
    for name, key in targets.items():
        body = eph[key]
        rows = []
        for d in dates:
            t = ts.utc(*d)
            jd_tt = float(t.tt)
            geom = body.at(t) - earth.at(t)
            gra, gdec, gdist = geom.radec()          # ICRS, same-time
            glon, glat = icrs_to_ecliptic_date(gra.radians, gdec.radians, jd_tt)

            app = earth.at(t).observe(body).apparent()  # light-time+aberration
            ara, adec, adist = app.radec()
            alon, alat = icrs_to_ecliptic_date(ara.radians, adec.radians, jd_tt,
                                               apparent=True)

            rows.append({
                "geom_lon": glon, "geom_lat": glat, "geom_dist": gdist.au,
                "app_lon": alon, "app_lat": alat, "app_dist": adist.au,
            })
        ref["bodies"][name] = rows

    with open(out, "w") as f:
        json.dump(ref, f, indent=1)
    print("wrote", out)


if __name__ == "__main__":
    main()
