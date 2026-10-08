#!/usr/bin/env python3
"""Regenerate the test reference data (development only).

This uses Skyfield with a JPL DE421 kernel to produce geometric and
apparent geocentric ecliptic-of-date positions for the Sun, Moon, Pluto
and the eight planets.  The DE kernel is NOT part of this library and is
never shipped; it is only used to produce tests/reference_data.h.

    python3 -m venv /tmp/eph && /tmp/eph/bin/pip install skyfield
    /tmp/eph/bin/python tools/make_reference.py /tmp/ref.json
    python3 tools/gen_reference.py /tmp/ref.json

DE421 spans roughly 1899-2053, hence the date range below.
"""
import json
import sys

from skyfield.api import load


def main():
    out = sys.argv[1] if len(sys.argv) > 1 else "/tmp/ref.json"
    eph = load("de421.bsp")
    ts = load.timescale()
    earth = eph["earth"]

    targets = {
        "Mercury": "mercury barycenter",
        "Venus": "venus barycenter",
        "Mars": "mars",
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
        (1930, 3, 3, 3.0),
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
            geom = body.at(t) - earth.at(t)
            lat, lon, dist = geom.ecliptic_latlon(epoch="date")
            app = earth.at(t).observe(body).apparent()
            la, lo, di = app.ecliptic_latlon(epoch="date")
            rows.append({
                "geom_lon": lon.degrees % 360, "geom_lat": lat.degrees,
                "geom_dist": dist.au,
                "app_lon": lo.degrees % 360, "app_lat": la.degrees,
                "app_dist": di.au,
            })
        ref["bodies"][name] = rows

    with open(out, "w") as f:
        json.dump(ref, f, indent=1)
    print("wrote", out)


if __name__ == "__main__":
    main()
