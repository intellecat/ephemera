#!/usr/bin/env python3
"""Regenerate tests/houses_ref.h (development only): house cusps, angles
and ayanamsa precession from an independent Python implementation of the
published definitions, for random dates in 1900-2030 and latitudes within
+-65 degrees, plus polar cases for the systems defined there.

Everything is done with vectors on the celestial sphere: the cusps are
intersections of the ecliptic with great circles (the horizon, house
circles through the north and south points, hour circles), Placidus by
root finding on the semi-arc condition.  Sidereal time, obliquity,
nutation and precession come from ERFA (IAU 2006/2000A).  The conventions
for the polar regions follow src/houses.cpp.

    python3 tools/make_houses_reference.py
"""
import math
import os

import erfa
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
POLAR_OK = (1, 2, 5, 6, 7, 9, 10, 11, 12, 13)  # defined in the polar circles
SYSTEMS = range(1, 14)
AYANAMSA_DATES = [2341972.5, 2415020.5, 2451545.0, 2469807.5, 2597641.5]
DELTA_T = 69.0 / 86400.0  # days; only moves the obliquity by ~1e-6"
TAU = 2 * math.pi


def norm(a):
    return a % TAU


def unit(v):
    return v / np.linalg.norm(v)


def equator_point(ra, dec=0.0):
    return np.array([math.cos(dec) * math.cos(ra), math.cos(dec) * math.sin(ra),
                     math.sin(dec)])


class Sky:
    """The sphere at one instant and place (equatorial frame of date)."""

    def __init__(self, jd_ut, lat_deg, lon_deg):
        tt = jd_ut + DELTA_T
        dpsi, deps = erfa.nut06a(tt, 0.0)
        self.eps = erfa.obl06(tt, 0.0) + deps
        self.ramc = norm(erfa.gst06a(jd_ut, 0.0, tt, 0.0) + math.radians(lon_deg))
        self.phi = math.radians(lat_deg)
        se, ce = math.sin(self.eps), math.cos(self.eps)
        self.pole = np.array([0.0, -se, ce])          # ecliptic north pole

    def longitude(self, v):
        """Ecliptic longitude of a direction."""
        se, ce = math.sin(self.eps), math.cos(self.eps)
        return norm(math.atan2(v[1] * ce + v[2] * se, v[0]))

    def ecliptic(self, lon):
        se, ce = math.sin(self.eps), math.cos(self.eps)
        return np.array([math.cos(lon), math.sin(lon) * ce, math.sin(lon) * se])

    def meet(self, normal):
        """The two points where the ecliptic meets the great circle with
        this pole."""
        d = unit(np.cross(self.pole, normal))
        return d, -d

    def north_point(self, phi, ramc):
        """North point of the horizon at latitude phi, sidereal time ramc."""
        return equator_point(ramc + math.pi, math.pi / 2 - phi)

    def zenith(self, phi, ramc):
        return equator_point(ramc, phi)

    def horizon_points(self, phi, ramc):
        """(east-side, west-side) intersections of the ecliptic with the
        horizon at latitude phi."""
        a, b = self.meet(self.zenith(phi, ramc))
        east = equator_point(ramc + math.pi / 2)
        return (a, b) if np.dot(a, east) > 0 else (b, a)

    def ascendant(self, ramc=None, phi=None):
        ramc = self.ramc if ramc is None else ramc
        phi = self.phi if phi is None else phi
        return self.longitude(self.horizon_points(phi, ramc)[0])

    def mc(self):
        """The ecliptic on the meridian, on the hour circle of the RAMC."""
        normal = np.cross(np.array([0.0, 0.0, 1.0]), equator_point(self.ramc))
        a, b = self.meet(normal)
        return self.longitude(a if np.dot(a, equator_point(self.ramc)) > 0 else b)

    def by_right_ascension(self, ra):
        """Ecliptic point on the hour circle of right ascension ra."""
        x = equator_point(ra)
        a, b = self.meet(np.cross(np.array([0.0, 0.0, 1.0]), x))
        return self.longitude(a if np.dot(a, x) > 0 else b)

    def on_circle(self, n, p):
        """Ecliptic point on the great circle through n and p, on the half
        from n through p (p perpendicular to n here)."""
        a, b = self.meet(np.cross(n, p))
        return a if np.dot(a, p) > 0 else b

    def declination(self, lon):
        return math.asin(math.sin(self.eps) * math.sin(lon))

    def right_ascension(self, lon):
        v = self.ecliptic(lon)
        return norm(math.atan2(v[1], v[0]))

    def semi_arc(self, lon):
        x = -math.tan(self.phi) * math.tan(self.declination(lon))
        return None if abs(x) >= 1 else math.acos(x)


def east_of(sky, p, ramc):
    """Of the two opposite points of a house circle, the convention of
    src/houses.cpp: the one with right ascension 0..180 degrees east of
    the RAMC."""
    lon = sky.longitude(p)
    ra = sky.right_ascension(lon)
    return lon if math.sin(ra - ramc) >= 0 else norm(lon + math.pi)


def ordered(seq, sign):
    return all(norm(sign * (b - a)) < math.pi for a, b in zip(seq, seq[1:]))


def consistent_mc(mc, asc, c):
    """Polar convention of src/houses.cpp: the MC (or its opposite) that
    keeps houses 10..4 in order around the ecliptic."""
    if ordered([mc, c[0], c[1], asc, c[2], c[3], mc + math.pi], +1):
        return mc
    if ordered([mc + math.pi, c[0], c[1], asc, c[2], c[3], mc], -1):
        return norm(mc + math.pi)
    return mc


def angles(sky):
    """Asc, MC with the convention that the Asc is 0..180 degrees after the
    MC along the ecliptic."""
    mc = sky.mc()
    east, west = sky.horizon_points(sky.phi, sky.ramc)
    asc = sky.longitude(east)
    if norm(asc - mc) >= math.pi:
        asc = sky.longitude(west)
    return asc, mc


def complete(asc, mc, c11, c12, c2, c3):
    cusps = [0.0] * 12
    for i, v in zip((0, 1, 2, 9, 10, 11), (asc, c2, c3, mc, c11, c12)):
        cusps[i] = norm(v)
        cusps[(i + 6) % 12] = norm(v + math.pi)
    return cusps


def bisect(f, a, b):
    fa = f(a)
    for _ in range(200):
        m = (a + b) / 2
        fm = f(m)
        if (fm < 0) == (fa < 0):
            a, fa = m, fm
        else:
            b = m
    return (a + b) / 2


def placidus(sky, asc, mc):
    """Cusps 11, 12 (2, 3): ecliptic points whose angle east of the upper
    meridian is 1/3, 2/3 of their diurnal semi-arc (the semi-arc plus 1/3,
    2/3 of the nocturnal one), found by root finding along the ecliptic."""
    def solve(target, start, span):
        def f(s):
            lon = start + s
            dsa = sky.semi_arc(lon)
            east = norm(sky.right_ascension(lon) - sky.ramc)
            return east - target(dsa)
        xs = np.linspace(1e-6, span - 1e-6, 721)
        vs = [f(x) for x in xs]
        for i in range(720):
            if (vs[i] < 0) != (vs[i + 1] < 0) and abs(vs[i] - vs[i + 1]) < 1:
                return norm(start + bisect(f, xs[i], xs[i + 1]))
        raise ValueError("no Placidus cusp")
    up = norm(asc - mc)
    down = norm(mc + math.pi - asc)
    return complete(
        asc, mc,
        solve(lambda d: d / 3, mc, up), solve(lambda d: 2 * d / 3, mc, up),
        solve(lambda d: math.pi - 2 * (math.pi - d) / 3, asc, down),
        solve(lambda d: math.pi - (math.pi - d) / 3, asc, down))


def cusps_for(sky, system):
    asc, mc = angles(sky)
    ramc, phi, third = sky.ramc, sky.phi, math.pi / 6
    if system == 1:                                   # Equal
        return [norm(asc + i * third) for i in range(12)]
    if system == 2:                                   # Porphyry
        d1, d2 = norm(asc - mc) / 3, norm(mc + math.pi - asc) / 3
        return complete(asc, mc, mc + d1, mc + 2 * d1, asc + d2, asc + 2 * d2)
    if system == 3:
        return placidus(sky, asc, mc)
    if system == 4:                                   # Koch
        dsa = sky.semi_arc(mc)
        return complete(asc, mc, sky.ascendant(ramc - 2 * dsa / 3),
                        sky.ascendant(ramc - dsa / 3), sky.ascendant(ramc + dsa / 3),
                        sky.ascendant(ramc + 2 * dsa / 3))
    if system in (5, 6, 8):
        # House circles through the north and south points of a horizon:
        # Regiomontanus through equator points 30k degrees east of the
        # meridian, Campanus through prime-vertical points 30k degrees
        # from the zenith, Topocentric (Polich-Page) through the equator
        # points with the north point of a horizon of pole height
        # tan P = tan(phi) k / 3.
        n = sky.north_point(phi, ramc)
        zen = sky.zenith(phi, ramc)
        east = equator_point(ramc + math.pi / 2)
        c = []
        for k in (1, 2, 4, 5):
            h = k * third
            if system == 5:
                p = sky.on_circle(n, equator_point(ramc + h))
            elif system == 6:
                p = sky.on_circle(n, math.cos(h) * zen + math.sin(h) * east)
            else:
                pole = math.atan(math.tan(phi) * min(k, 6 - k) / 3)
                x = equator_point(ramc + h)
                p = sky.on_circle(equator_point(ramc + h + math.pi / 2,
                                                math.pi / 2 - pole), x)
            c.append(east_of(sky, p, ramc))
        return complete(asc, consistent_mc(mc, asc, c), *c)
    if system == 7:                                   # Alcabitius
        dsa = sky.semi_arc(asc)
        nsa = math.pi - dsa
        r = sky.by_right_ascension
        return complete(asc, mc, r(ramc + dsa / 3), r(ramc + 2 * dsa / 3),
                        r(ramc + dsa + nsa / 3), r(ramc + dsa + 2 * nsa / 3))
    if system == 9:                                   # Vehlow
        return [norm(asc + (i - 0.5) * third) for i in range(12)]
    if system == 10:                                  # Equal from the MC
        return [norm(mc + (i + 3) * third) for i in range(12)]
    if system == 11:                                  # Equal from 0 Aries
        return [i * third for i in range(12)]
    if system == 12:                                  # Meridian
        return [sky.by_right_ascension(ramc + (i + 3) * third) for i in range(12)]
    if system == 13:                                  # Morinus
        return [sky.longitude(equator_point(ramc + (i + 3) * third)) for i in range(12)]
    raise ValueError(system)


def all_angles(sky):
    """Asc, MC, ARMC, Vertex, equatorial Asc, co-Asc (Koch), co-Asc
    (Munkasey), polar Asc."""
    asc, mc = angles(sky)
    # Vertex: the ecliptic on the prime vertical (pole: the north point),
    # in the west.
    a, b = sky.meet(sky.north_point(sky.phi, sky.ramc))
    west = equator_point(sky.ramc - math.pi / 2)
    vertex = sky.longitude(a if np.dot(a, west) > 0 else b)
    equatorial = sky.ascendant(sky.ramc, 0.0)
    co_koch = sky.ascendant(sky.ramc, -sky.phi)
    # Munkasey: the horizon of the co-latitude; of its two points the one
    # whose longitude has the sign of sin(RAMC + 90) in its sine (the
    # convention of the closed formula).
    p, q = sky.meet(sky.zenith(math.pi / 2 - sky.phi, sky.ramc))
    lp = sky.longitude(p)
    co_munkasey = lp if math.sin(lp) * math.cos(sky.ramc) >= 0 else sky.longitude(q)
    return [asc, mc, sky.ramc, vertex, equatorial, co_koch, co_munkasey,
            norm(co_koch + math.pi)]


def main():
    rng = np.random.default_rng(11)
    rows, angle_rows = [], []
    for k in range(48):
        polar = k >= 40
        jd = rng.uniform(2415020.5, 2462502.5)
        lat = (rng.choice([-1, 1]) * rng.uniform(67, 80) if polar
               else rng.uniform(-65, 65))
        lon = rng.uniform(-180, 180)
        sky = Sky(jd, lat, lon)
        if not polar:
            angle_rows.append((jd, lat, lon) + tuple(math.degrees(a) for a in all_angles(sky)))
        for system in SYSTEMS:
            if polar and system not in POLAR_OK:
                continue
            rows.append((jd, lat, lon, system) +
                        tuple(math.degrees(c) for c in cusps_for(sky, system)))
    # General precession in longitude (IAU 2006, ERFA p06e) since J2000.
    precession = [math.degrees(erfa.p06e(d, 0.0)[12]) for d in AYANAMSA_DATES]
    nutation_1956 = math.degrees(erfa.nut06a(2435553.5, 0.0)[0])

    out = os.path.join(ROOT, "tests", "houses_ref.h")
    with open(out, "w") as f:
        f.write("/* Generated by tools/make_houses_reference.py: house cusps "
                "(degrees) from an\n * independent implementation of the "
                "definitions with ERFA (IAU 2006/2000A),\n * for JD UT, "
                "latitude, longitude, system. */\n")
        f.write("#pragma once\n")
        f.write("struct HouseRef { double jd_ut, lat, lon; int system; "
                "double cusp[12]; };\n\n")
        f.write("static const HouseRef kHouseRef[%d] = {\n" % len(rows))
        for r in rows:
            f.write("    {%.10f, %.10f, %.10f, %d, {%s}},\n"
                    % (r[0], r[1], r[2], r[3],
                       ", ".join("%.10f" % c for c in r[4:])))
        f.write("};\n\n")
        f.write("/* Asc, MC, ARMC, Vertex, equatorial Asc, co-Asc Koch, co-Asc\n"
                " * Munkasey, polar Asc (degrees). */\n")
        f.write("struct AngleRef { double jd_ut, lat, lon, value[8]; };\n\n")
        f.write("static const AngleRef kAngleRef[%d] = {\n" % len(angle_rows))
        for r in angle_rows:
            f.write("    {%.10f, %.10f, %.10f, {%s}},\n"
                    % (r[0], r[1], r[2], ", ".join("%.10f" % c for c in r[3:])))
        f.write("};\n\n")
        f.write("/* General precession in longitude since J2000 (degrees; ERFA "
                "p06e) at 1700,\n * 1900, 2000, 2050 and 2400 (TT). */\n")
        f.write("static const double kAyanamsaDates[5] = {%s};\n"
                % ", ".join("%.1f" % d for d in AYANAMSA_DATES))
        f.write("static const double kPrecession[5] = {%s};\n"
                % ", ".join("%.12f" % p for p in precession))
        f.write("/* Nutation in longitude (degrees; ERFA nut06a) on 1956-03-21 "
                "0h TT. */\n")
        f.write("static const double kNutation1956 = %.12f;\n" % nutation_1956)
    print("wrote", out, len(rows), "rows")


if __name__ == "__main__":
    main()
