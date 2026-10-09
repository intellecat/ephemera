#!/usr/bin/env python3
"""Regenerate tests/eclipse_ref.h (development only): eclipse circumstances
from JPL DE440 through Skyfield, computed directly from the apparent Sun
and Moon rather than from Besselian elements:

  * local circumstances: topocentric positions for an observer on the
    WGS84 ellipsoid; the contacts are where the centres' separation equals
    the sum (1st/4th) or difference (2nd/3rd) of the semidiameters;
  * global contacts: the first and last instants when any point of the
    ellipsoid is inside the penumbral (P1/P4) or umbral/antumbral (U1/U4)
    cone, found by minimising over the whole surface, and when the shadow
    axis meets the ellipsoid (centre line);
  * lunar contacts: the Moon's centre against the Earth's shadow at its
    distance, enlarged by Danjon's rule (NASA's canon).

Same constants as the library: Sun radius 696000 km, Moon radius k =
0.272488 (penumbra, 1st/4th contacts) and 0.272281 (umbra, 2nd/3rd) Earth
equatorial radii, WGS84.  UT for the topocentric places uses the library's
Delta-T table (data/delta_t_data.cpp).  The eclipses are taken from
tests/eclipse_nasa_ref.h.  Needs skyfield, numpy and a DE440 kernel:

    python3 tools/make_eclipse_reference.py [path/to/de440.bsp]
"""
import math
import os
import re
import sys

import numpy as np
from skyfield.api import Loader, wgs84

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
A_KM = 6378.137
F = 1 / 298.257223563
B_KM = A_KM * (1 - F)
SUN_KM = 696000.0
K1, K2 = 0.272488, 0.2722810
DAY = 86400.0

# Places: (east longitude, latitude, metres, first year).  The last five
# lie on central paths: 1999-08-11 (total), 2009-07-22 (total), 2012-05-20
# (annular), 2023-10-14 (annular), 2024-04-08 (total).
PLACES = [(-104.15, 25.29, 0.0, 1950), (13.40, 52.52, 34.0, 1980),
          (116.40, 39.90, 44.0, 2010), (-70.65, -33.45, 570.0, 2040),
          (151.21, -33.87, 20.0, 2070), (11.58, 48.14, 520.0, 1992),
          (121.47, 31.23, 4.0, 2003), (139.69, 35.69, 40.0, 2006),
          (-106.65, 35.08, 1500.0, 2016), (-106.42, 23.22, 10.0, 2018)]
PER_PLACE = 6


def nasa():
    text = open(os.path.join(ROOT, "tests", "eclipse_nasa_ref.h")).read()
    solar_part, lunar_part = text.split("kNasaLunar")
    num = r"\{(2\d+\.\d+), (\d+), (\d+), (-?[\d.]+)"
    return ([(float(m[0]), int(m[2]), float(m[3])) for m in re.findall(num, solar_part)],
            [(float(m[0]), int(m[2]), float(m[3])) for m in re.findall(num, lunar_part)])


def delta_t_table():
    text = open(os.path.join(ROOT, "data", "delta_t_data.cpp")).read()
    vals = [float(v) for v in re.findall(r"(-?\d+\.\d+)f", text)]
    return np.float32(vals).astype(float)


DT = delta_t_table()


def delta_t(jd):
    """The library's Delta-T (src/time.cpp), 1620-2150."""
    y = 2000.0 + (jd - 2451545.0) / 365.25
    i = int(math.floor(y)) - 1620
    return DT[i] + (y - math.floor(y)) * (DT[i + 1] - DT[i])


def root(f, a, b, tol=1e-8):
    """Zero of f in [a, b] (sign change required), Illinois."""
    fa, fb = f(a), f(b)
    if (fa < 0) == (fb < 0):
        return 0.0
    side = 0
    while abs(b - a) > tol:
        c = (fa * b - fb * a) / (fa - fb)
        fc = f(c)
        if fc == 0:
            return c
        if (fc < 0) == (fb < 0):
            b, fb = c, fc
            if side == -1:
                fa /= 2
            side = -1
        else:
            a, fa = c, fc
            if side == 1:
                fb /= 2
            side = 1
    return (a + b) / 2


def golden(f, a, b, tol=1e-8):
    g = (math.sqrt(5) - 1) / 2
    c, d = b - g * (b - a), a + g * (b - a)
    fc, fd = f(c), f(d)
    while b - a > tol:
        if fc < fd:
            b, d, fd = d, c, fc
            c = b - g * (b - a)
            fc = f(c)
        else:
            a, c, fc = c, d, fd
            d = a + g * (b - a)
            fd = f(d)
    return (a + b) / 2


class Sky:
    def __init__(self, kernel):
        load = Loader(os.path.dirname(kernel) or ".")
        self.eph = load(os.path.basename(kernel))
        self.earth, self.sun, self.moon = self.eph["earth"], self.eph["sun"], self.eph["moon"]
        self.load = load
        self.ts = load.timescale()

    def geocentric(self, jd):
        """Apparent geocentric Sun and Moon (km), true equator of date."""
        t = self.ts.tt_jd(jd)
        e = self.earth.at(t)
        s = e.observe(self.sun).apparent().position.km
        m = e.observe(self.moon).apparent().position.km
        return t.M @ s, t.M @ m

    # ---- local circumstances ------------------------------------------

    def place(self, lon, lat, h, jd):
        """Topocentric state function for a place, UT from the library's
        Delta-T at jd."""
        ts = self.load.timescale(delta_t=delta_t(jd))
        obs = self.earth + wgs84.latlon(lat, lon, elevation_m=h)

        def state(x):
            o = obs.at(ts.tt_jd(x))
            s = o.observe(self.sun).apparent()
            m = o.observe(self.moon).apparent()
            ds, dm = s.distance().km, m.distance().km
            sep = s.separation_from(m).radians
            rs = math.asin(SUN_KM / ds)
            alt = s.altaz()[0].radians
            return sep, rs, math.asin(K1 * A_KM / dm), math.asin(K2 * A_KM / dm), alt
        return state

    def visible(self, lon, lat, h, jd):
        """Whether the eclipse at jd is seen from the place with the Sun up
        at maximum (coarse)."""
        ts = self.load.timescale(delta_t=delta_t(jd))
        t = ts.tt_jd(jd + np.linspace(-0.3, 0.3, 61))
        o = (self.earth + wgs84.latlon(lat, lon, elevation_m=h)).at(t)
        s = o.observe(self.sun).apparent()
        m = o.observe(self.moon).apparent()
        g = (s.separation_from(m).radians - np.arcsin(SUN_KM / s.distance().km) -
             np.arcsin(K1 * A_KM / m.distance().km))
        i = int(np.argmin(g))
        return g[i] < -1e-4 and s.altaz()[0].degrees[i] > 2.0

    def local(self, lon, lat, h, jd):
        st = self.place(lon, lat, h, jd)
        outer = lambda x: (lambda v: v[0] - v[1] - v[2])(st(x))
        inner = lambda x: (lambda v: v[0] - abs(v[1] - v[3]))(st(x))
        tmax = golden(outer, jd - 0.25, jd + 0.25)
        c1, c4 = root(outer, tmax - 0.25, tmax), root(outer, tmax, tmax + 0.25)
        c2 = c3 = 0.0
        if inner(tmax) < 0:
            c2, c3 = root(inner, tmax - 0.05, tmax), root(inner, tmax, tmax + 0.05)
        sep, rs, rm1, rm2, alt = st(tmax)
        if sep <= abs(rs - rm2):
            mag = rm2 / rs
        else:
            mag = (rs + rm1 - sep) / (2 * rs)
        return tmax, c1, c2, c3, c4, mag, overlap(rs, rm1, sep), math.degrees(alt)

    # ---- global contacts ------------------------------------------------

    def surface_min(self, s, m, kind):
        """Minimum over the ellipsoid of separation - (rs + rm) (kind 0)
        or separation - |rs - rm| (kind 1), radians."""
        def g(beta, lam):
            p = np.stack([A_KM * np.cos(beta) * np.cos(lam),
                          A_KM * np.cos(beta) * np.sin(lam),
                          B_KM * np.sin(beta)], axis=-1)
            vs, vm = s - p, m - p
            ns, nm = np.linalg.norm(vs, axis=-1), np.linalg.norm(vm, axis=-1)
            sep = np.arctan2(np.linalg.norm(np.cross(vs, vm), axis=-1),
                             np.einsum("...i,...i", vs, vm))
            rs = np.arcsin(SUN_KM / ns)
            if kind == 0:
                return sep - rs - np.arcsin(K1 * A_KM / nm)
            return sep - np.abs(rs - np.arcsin(K2 * A_KM / nm))
        b, l = np.meshgrid(np.radians(np.arange(-89, 90, 2.0)),
                           np.radians(np.arange(0, 360, 2.0)), indexing="ij")
        v = g(b, l)
        i = np.unravel_index(np.argmin(v), v.shape)
        cb, cl, step = b[i], l[i], math.radians(2.0)
        offs = np.linspace(-2, 2, 5)
        best = v[i]
        while step > 1e-10:
            bb, ll = np.meshgrid(cb + offs * step, cl + offs * step, indexing="ij")
            vv = g(bb, ll)
            j = np.unravel_index(np.argmin(vv), vv.shape)
            if j == (2, 2):
                step /= 2
            cb, cl, best = bb[j], ll[j], vv[j]
        return best

    def axis_miss(self, s, m):
        """Distance of the shadow axis from the Earth's centre (km) in
        coordinates where the ellipsoid is a sphere, minus the radius."""
        st = np.diag([1.0, 1.0, A_KM / B_KM])
        s, m = st @ s, st @ m
        u = (m - s) / np.linalg.norm(m - s)
        return np.linalg.norm(np.cross(m, u)) - A_KM

    def global_contacts(self, jd):
        def fn(kind):
            def f(x):
                s, m = self.geocentric(x)
                return self.axis_miss(s, m) if kind == 2 else self.surface_min(s, m, kind)
            return f
        out = []
        for kind in range(3):
            f = fn(kind)
            if f(jd) < 0:
                out += [root(f, jd - 0.3, jd), root(f, jd, jd + 0.3)]
            else:
                out += [0.0, 0.0]
        return out

    # ---- lunar ------------------------------------------------------------

    def lunar(self, jd):
        def shadow(x):
            s, m = self.geocentric(x)
            ds, dm = np.linalg.norm(s), np.linalg.norm(m)
            sep = math.atan2(np.linalg.norm(np.cross(-s, m)), np.dot(-s, m))
            par, ps = math.asin(A_KM / dm), math.asin(A_KM / ds)
            rs = math.asin(SUN_KM / ds)
            moon = math.asin(K1 * A_KM / dm)
            pm = par * (1 + 1 / 85) * (1 - F / 2)
            return sep, pm + ps + rs, pm + ps - rs, moon
        tmax = golden(lambda x: shadow(x)[0], jd - 0.1, jd + 0.1)
        out = [tmax]
        for which in range(3):
            def f(x, which=which):
                sep, pen, umb, moon = shadow(x)
                return sep - (pen + moon if which == 0 else
                              umb + moon if which == 1 else umb - moon)
            if f(tmax) < 0:
                out += [root(f, tmax - 0.3, tmax), root(f, tmax, tmax + 0.3)]
            else:
                out += [0.0, 0.0]
        return out


def overlap(r1, r2, c):
    """Area of disc 2 over disc 1, as a fraction of disc 1."""
    if c >= r1 + r2:
        return 0.0
    if c <= abs(r1 - r2):
        return 1.0 if r2 >= r1 else r2 * r2 / (r1 * r1)
    a1 = math.acos((c * c + r1 * r1 - r2 * r2) / (2 * c * r1))
    a2 = math.acos((c * c + r2 * r2 - r1 * r1) / (2 * c * r2))
    return (r1 * r1 * (a1 - math.sin(2 * a1) / 2) +
            r2 * r2 * (a2 - math.sin(2 * a2) / 2)) / (math.pi * r1 * r1)


def main():
    sky = Sky(sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, "de440.bsp"))
    solar, lunar = nasa()

    local = []
    for k, (lon, lat, h, year) in enumerate(PLACES):
        start = 2451545.0 + (year - 2000) * 365.25
        n = 0
        for jd, _, _ in solar:
            if jd < start or n >= PER_PLACE:
                continue
            if sky.visible(lon, lat, h, jd):
                local.append((k, jd) + sky.local(lon, lat, h, jd))
                n += 1
        print("place", k, n, "eclipses", file=sys.stderr)

    # Global contacts: every 12th eclipse from 1950, the non-central ones
    # and the most grazing partial ones.
    recent = [e for e in solar if e[0] > 2433282.5]
    pick = set(range(0, len(recent), 12))
    pick |= {i for i, e in enumerate(recent) if e[1] & 32}
    pick |= set(sorted(range(len(recent)), key=lambda i: -abs(recent[i][2]))[:4])
    glob = [(recent[i][0],) + tuple(sky.global_contacts(recent[i][0])) for i in sorted(pick)]
    print(len(glob), "global", file=sys.stderr)

    lun = [e for e in lunar if e[0] > 2433282.5][::10]
    lun = [(e[0],) + tuple(sky.lunar(e[0])) for e in lun]
    print(len(lun), "lunar", file=sys.stderr)

    with open(os.path.join(ROOT, "tests", "eclipse_ref.h"), "w") as f:
        f.write("/* Generated by tools/make_eclipse_reference.py from JPL DE440\n"
                " * (Skyfield): eclipse circumstances computed directly from the\n"
                " * apparent Sun and Moon, TT; 0: the phase does not occur. */\n"
                "#pragma once\n\n")
        f.write("static const double kPlaces[][3] = {\n")
        for p in PLACES:
            f.write("    {%.4f, %.4f, %.1f},\n" % p[:3])
        f.write("};\n\n/* Local circumstances: place, the eclipse's greatest eclipse\n"
                " * (NASA), maximum, contacts, magnitude, obscuration, Sun's\n"
                " * altitude (degrees). */\n"
                "struct LocalRef { int place; double global_max, max, c1, c2, c3, c4, "
                "mag, obscuration, altitude; };\nstatic const LocalRef kLocal[] = {\n")
        for r in local:
            f.write("    {%d, %.6f, %.8f, %.8f, %.8f, %.8f, %.8f, %.6f, %.6f, %.5f},\n" % r)
        f.write("};\n\n/* Global contacts: penumbra, umbra/antumbra and axis on the\n"
                " * Earth. */\nstruct GlobalRef { double global_max, begin, end, "
                "central_begin, central_end, center_begin, center_end; };\n"
                "static const GlobalRef kGlobal[] = {\n")
        for r in glob:
            f.write("    {%.6f, %.8f, %.8f, %.8f, %.8f, %.8f, %.8f},\n" % r)
        f.write("};\n\n/* Lunar: NASA greatest eclipse, maximum, penumbral, partial and\n"
                " * total contacts. */\nstruct LunarRef { double nasa_max, max, pen_begin, "
                "pen_end, par_begin, par_end, tot_begin, tot_end; };\n"
                "static const LunarRef kLunar[] = {\n")
        for r in lun:
            f.write("    {%.6f, %.8f, %.8f, %.8f, %.8f, %.8f, %.8f, %.8f},\n" % r)
        f.write("};\n")
    print(len(local), "local,", len(glob), "global,", len(lun), "lunar")


if __name__ == "__main__":
    main()
