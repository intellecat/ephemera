"""The native library and its C types (ctypes).  Internal: use the functions
in the ``ephemera`` package, which take and return degrees."""

import ctypes
import os
import sys
from ctypes import POINTER, Structure, c_char_p, c_double, c_int, c_size_t, c_void_p

_HERE = os.path.dirname(os.path.abspath(__file__))


def _find_library() -> str:
    if sys.platform == "win32":
        names = ["ephemera_native.dll", "libephemera_native.dll"]
    elif sys.platform == "darwin":
        names = ["libephemera_native.dylib", "libephemera_native.so"]
    else:
        names = ["libephemera_native.so"]
    for name in names:
        path = os.path.join(_HERE, name)
        if os.path.exists(path):
            return path
    raise ImportError(f"ephemera: native library not found in {_HERE} (looked for {', '.join(names)})")


lib = ctypes.CDLL(_find_library())


class Position(Structure):
    _fields_ = [(n, c_double) for n in ("longitude", "latitude", "distance",
                                         "longitude_speed", "latitude_speed", "distance_speed")]


class Observer(Structure):
    _fields_ = [("longitude", c_double), ("latitude", c_double), ("elevation", c_double)]


class Angles(Structure):
    _fields_ = [(n, c_double) for n in ("ascendant", "mc", "armc", "vertex", "equatorial_ascendant",
                                         "co_ascendant_koch", "co_ascendant_munkasey",
                                         "polar_ascendant")]


class Phenomena(Structure):
    _fields_ = [(n, c_double) for n in ("phase_angle", "phase", "elongation", "diameter", "magnitude")]


class SolarEclipse(Structure):
    _fields_ = [("type", c_int), ("maximum", c_double), ("gamma", c_double), ("magnitude", c_double),
                ("latitude", c_double), ("longitude", c_double), ("saros", c_int),
                ("partial_begin", c_double), ("partial_end", c_double),
                ("central_begin", c_double), ("central_end", c_double),
                ("center_line_begin", c_double), ("center_line_end", c_double)]


class LocalEclipseState(Structure):
    _fields_ = [("type", c_int), ("magnitude", c_double), ("obscuration", c_double),
                ("altitude", c_double)]


class LocalSolarEclipse(Structure):
    _fields_ = [("type", c_int), ("maximum", c_double),
                ("first_contact", c_double), ("second_contact", c_double),
                ("third_contact", c_double), ("fourth_contact", c_double),
                ("magnitude", c_double), ("obscuration", c_double), ("altitude", c_double),
                ("saros", c_int)]


class LunarEclipse(Structure):
    _fields_ = [("type", c_int), ("maximum", c_double), ("gamma", c_double),
                ("umbral_magnitude", c_double), ("penumbral_magnitude", c_double), ("saros", c_int),
                ("penumbral_begin", c_double), ("penumbral_end", c_double),
                ("partial_begin", c_double), ("partial_end", c_double),
                ("total_begin", c_double), ("total_end", c_double)]


class Orbit(Structure):
    _fields_ = [("semi_major_axis", c_double), ("eccentricity", c_double), ("inclination", c_double),
                ("node", c_double), ("perihelion", c_double),
                ("heliocentric", Position * 4), ("geocentric", Position * 4)]


D = c_double
I = c_int
PD = POINTER(c_double)
PI = POINTER(c_int)
PPOS = POINTER(Position)
POBS = POINTER(Observer)

_SIGNATURES = {
    # time
    "eph_julian_date": (D, [I, I, D]),
    "eph_calendar_date": (None, [D, PI, PI, PD]),
    "eph_delta_t": (D, [D]),
    "eph_tdb_minus_tt": (D, [D]),
    "eph_ut_to_tt": (D, [D, D]),
    "eph_tt_to_ut": (D, [D, D]),
    "eph_utc_to_jd": (I, [I, I, I, I, I, D, PD, PD]),
    "eph_sidereal_time": (D, [D]),
    "eph_equation_of_time": (I, [D, PD]),
    # frames
    "eph_mean_obliquity": (D, [D]),
    "eph_nutation": (None, [D, PD, PD]),
    "eph_ecliptic_to_equatorial": (None, [D, D, D, PD, PD]),
    # positions
    "eph_position_tt": (I, [I, D, PPOS]),
    "eph_position_apparent_tt": (I, [I, D, PPOS]),
    "eph_position_ut": (I, [I, D, D, PPOS]),
    "eph_position_apparent_ut": (I, [I, D, D, PPOS]),
    "eph_heliocentric_tt": (I, [I, D, PD, PD, PD]),
    "eph_position_topocentric_ut": (I, [I, D, D, POBS, PPOS]),
    "eph_body_name": (c_char_p, [I]),
    # sidereal zodiac and houses
    "eph_ayanamsa": (D, [I, D]),
    "eph_ayanamsa_name": (c_char_p, [I]),
    "eph_ayanamsa_custom": (D, [D, D, D]),
    "eph_ayanamsa_lahiri": (D, [D]),
    "eph_asc_mc": (I, [D, POBS, PD, PD]),
    "eph_houses": (I, [D, POBS, I, PD]),
    "eph_house_angles": (I, [D, POBS, POINTER(Angles)]),
    "eph_house_position": (D, [PD, D]),
    # phenomena, eclipses
    "eph_phenomena_tt": (I, [I, D, POINTER(Phenomena)]),
    "eph_solar_eclipse_next": (I, [D, I, POINTER(SolarEclipse)]),
    "eph_solar_eclipse_where": (I, [D, PD, PD]),
    "eph_solar_eclipse_how": (I, [D, POBS, POINTER(LocalEclipseState)]),
    "eph_solar_eclipse_local": (I, [D, POBS, I, POINTER(LocalSolarEclipse)]),
    "eph_lunar_eclipse_next": (I, [D, I, POINTER(LunarEclipse)]),
    # searches
    "eph_longitude_crossing": (I, [I, D, D, I, D, PD]),
    "eph_longitude_crossing_ut": (I, [I, D, D, I, D, PD]),
    "eph_heliocentric_crossing": (I, [I, D, D, I, D, PD]),
    "eph_longitude_crossing_sidereal": (I, [I, D, I, D, I, D, PD]),
    "eph_longitude_crossing_sidereal_custom": (I, [I, D, D, D, D, I, D, PD]),
    "eph_heliocentric_crossing_sidereal": (I, [I, D, I, D, I, D, PD]),
    "eph_heliocentric_crossing_sidereal_custom": (I, [I, D, D, D, D, I, D, PD]),
    "eph_node_crossing": (I, [I, D, I, D, PD, PD]),
    "eph_station": (I, [I, D, I, D, PD]),
    # minor planets, orbits
    "eph_minor_count": (I, []),
    "eph_minor_find": (I, [c_char_p]),
    "eph_minor_number": (I, [I]),
    "eph_minor_name": (c_char_p, [I]),
    "eph_orbit_tt": (I, [I, D, I, POINTER(Orbit)]),
    # utilities
    "eph_midpoint": (D, [D, D]),
    "eph_split_longitude": (None, [D, PI, PI, PI, PD]),
    # fixed stars
    "eph_star_count": (I, []),
    "eph_star_find": (I, [c_char_p]),
    "eph_star_name": (c_char_p, [I]),
    "eph_star_designation": (c_char_p, [I]),
    "eph_star_magnitude": (D, [I]),
    "eph_star_position_apparent_ut": (I, [I, D, D, PPOS]),
    # data files
    "eph_load_data": (I, [c_void_p, c_size_t]),
    "eph_load_file": (I, [c_char_p]),
    "eph_load_embedded": (I, []),
    "eph_unload_data": (None, []),
}

for _name, (_restype, _argtypes) in _SIGNATURES.items():
    _fn = getattr(lib, _name)
    _fn.restype = _restype
    _fn.argtypes = _argtypes
