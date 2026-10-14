"""Ephemera: a pocket-sized, MIT-licensed ephemeris for astrology and astronomy.

    import ephemera as eph

    jd = eph.julian_date(2026, 1, 1.5)               # 1 Jan 2026, 12:00 UT
    moon = eph.position(eph.Body.MOON, jd)           # degrees, AU, per day
    cusps = eph.houses(jd, -0.1278, 51.5074, eph.HouseSystem.PLACIDUS)
    lahiri = eph.ayanamsa(eph.Ayanamsa.LAHIRI, eph.ut_to_tt(jd, 2026))

The 1" planetary series with Chiron (1550-2650) is built in.  More data
files come with the package (``load_package_data('stars')``) or from the
repository (``load_data(path)``).

Angles are degrees, distances AU, times Julian days (UT or TT as named);
calls that find nothing return None.  The data loaded is shared by the
whole process: load files before computing from several threads.
"""

from __future__ import annotations

import datetime as _dt
import math
import os
from ctypes import byref, c_double, c_int
from dataclasses import dataclass
from enum import IntEnum, IntFlag
from typing import List, NamedTuple, Optional, Union

from . import _native as _n
from ._native import lib as _lib

__version__ = "0.4.1"

DEG = 180.0 / math.pi
_RAD = math.pi / 180.0
MINOR_BASE = 1000
"""Body id of the first loaded minor planet (``minor_body`` gives the id)."""


class Body(IntEnum):
    MERCURY = 0
    VENUS = 1
    EARTH = 2
    MARS = 3
    JUPITER = 4
    SATURN = 5
    URANUS = 6
    NEPTUNE = 7
    PLUTO = 8
    SUN = 9
    MOON = 10
    CERES = 11
    PALLAS = 12
    JUNO = 13
    VESTA = 14
    CHIRON = 15
    MEAN_NODE = 16
    TRUE_NODE = 17
    MEAN_LILITH = 18
    OSCU_LILITH = 19
    # The Hamburg school's Uranian points (Witte, Sieggrün).
    CUPIDO = 20
    HADES = 21
    ZEUS = 22
    KRONOS = 23
    APOLLON = 24
    ADMETOS = 25
    VULKANUS = 26
    POSEIDON = 27


class Ayanamsa(IntEnum):
    LAHIRI = 0
    FAGAN_BRADLEY = 1
    RAMAN = 2
    KRISHNAMURTI = 3
    YUKTESHWAR = 4
    DE_LUCE = 5
    USHA_SHASHI = 6
    DJWHAL_KHUL = 7
    JN_BHASIN = 8
    BABYLONIAN_KUGLER1 = 9
    BABYLONIAN_KUGLER2 = 10
    BABYLONIAN_KUGLER3 = 11
    BABYLONIAN_HUBER = 12
    BABYLONIAN_ETA_PISCIUM = 13
    BABYLONIAN_ALDEBARAN = 14
    HIPPARCHOS = 15
    SASSANIAN = 16
    J2000 = 17
    J1900 = 18
    B1950 = 19
    SURYASIDDHANTA = 20
    SURYASIDDHANTA_MEAN_SUN = 21
    ARYABHATA = 22
    ARYABHATA_MEAN_SUN = 23
    SS_REVATI = 24
    SS_CITRA = 25
    SKYDRAM = 26
    ARYABHATA_522 = 27
    BABYLONIAN_BRITTON = 28
    GALACTIC_EQUATOR_FIORENZA = 29
    VETTIUS_VALENS = 30
    LAHIRI_1940 = 31
    LAHIRI_VP285 = 32
    KRISHNAMURTI_SENTHILATHIBAN = 33
    LAHIRI_ICRC = 34
    # Star-based (need the stars: load_package_data('stars')).
    TRUE_CITRA = 35
    TRUE_REVATI = 36
    TRUE_PUSHYA = 37
    TRUE_MULA = 38
    GALACTIC_CENTER = 39
    GALACTIC_CENTER_BRAND = 40
    GALACTIC_CENTER_COCHRANE = 41


class HouseSystem(IntEnum):
    WHOLE_SIGN = 0
    EQUAL = 1
    PORPHYRY = 2
    PLACIDUS = 3
    KOCH = 4
    REGIOMONTANUS = 5
    CAMPANUS = 6
    ALCABITIUS = 7
    TOPOCENTRIC = 8
    VEHLOW = 9
    EQUAL_MC = 10
    EQUAL_ARIES = 11
    MERIDIAN = 12
    MORINUS = 13


class OrbitMethod(IntEnum):
    MEAN = 0
    """Secular mean elements (Simon et al. 1994 / Meeus ch. 31): Mercury..Neptune."""
    OSCULATING = 1
    """From the position and velocity: the planets but the Earth, Pluto, minor bodies."""


class EclipseType(IntFlag):
    TOTAL = 1
    ANNULAR = 2
    HYBRID = 4
    PARTIAL = 8
    PENUMBRAL = 16
    NON_CENTRAL = 32


# ---- results ---------------------------------------------------------------


@dataclass(frozen=True)
class Position:
    longitude: float
    latitude: float
    distance: float
    longitude_speed: float
    latitude_speed: float
    distance_speed: float

    @property
    def retrograde(self) -> bool:
        return self.longitude_speed < 0


class CalendarDate(NamedTuple):
    year: int
    month: int
    day: float  # with the fraction of the day


class JulianDates(NamedTuple):
    jd_tt: float
    jd_ut: float


class Nutation(NamedTuple):
    longitude: float
    obliquity: float


class Equatorial(NamedTuple):
    right_ascension: float
    declination: float


class Heliocentric(NamedTuple):
    longitude: float
    latitude: float
    distance: float


class AscMc(NamedTuple):
    asc: float
    mc: float


class HousesInfo(NamedTuple):
    cusps: List[float]
    fallback: bool  # the system is undefined here (polar): Porphyry was used


@dataclass(frozen=True)
class HouseAngles:
    ascendant: float
    mc: float
    armc: float
    vertex: float
    equatorial_ascendant: float
    co_ascendant_koch: float
    co_ascendant_munkasey: float
    polar_ascendant: float


@dataclass(frozen=True)
class Phenomena:
    phase_angle: float
    phase: float  # illuminated fraction, 0..1
    elongation: float
    diameter: float
    magnitude: float


class CustomAyanamsa(NamedTuple):
    """A user-defined ayanamsa: its value (degrees) at a TT Julian day."""
    jd_t0: float
    value_at_t0: float


class NodeCrossing(NamedTuple):
    jd: float
    node: str  # 'ascending' or 'descending'
    longitude: float


class Station(NamedTuple):
    jd: float
    turns: str  # 'retrograde' or 'direct'


@dataclass(frozen=True)
class OrbitPoints:
    ascending_node: Position
    descending_node: Position
    perihelion: Position
    aphelion: Position


@dataclass(frozen=True)
class Orbit:
    """A heliocentric orbit, mean ecliptic and equinox of date."""
    semi_major_axis: float
    eccentricity: float
    inclination: float
    node: float
    perihelion: float
    heliocentric: OrbitPoints
    geocentric: OrbitPoints


@dataclass(frozen=True)
class SolarEclipse:
    """Times are TT Julian days; 0 marks a phase that does not occur."""
    type: EclipseType
    maximum: float
    gamma: float
    magnitude: float
    latitude: float
    longitude: float
    saros: int
    partial_begin: float
    partial_end: float
    central_begin: float
    central_end: float
    center_line_begin: float
    center_line_end: float


@dataclass(frozen=True)
class LocalSolarEclipse:
    type: EclipseType
    maximum: float
    first_contact: float
    second_contact: float
    third_contact: float
    fourth_contact: float
    magnitude: float
    obscuration: float
    altitude: float  # of the Sun at maximum; below 0: not visible
    saros: int


@dataclass(frozen=True)
class LocalEclipseState:
    type: EclipseType  # 0: none
    magnitude: float
    obscuration: float
    altitude: float


@dataclass(frozen=True)
class LunarEclipse:
    type: EclipseType
    maximum: float
    gamma: float
    umbral_magnitude: float
    penumbral_magnitude: float
    saros: int
    penumbral_begin: float
    penumbral_end: float
    partial_begin: float
    partial_end: float
    total_begin: float
    total_end: float


class GeoPoint(NamedTuple):
    latitude: float
    longitude: float


class MinorPlanet(NamedTuple):
    body: int
    number: int
    name: str


class SplitLongitude(NamedTuple):
    sign: int  # 0 = Aries .. 11 = Pisces
    degrees: int
    minutes: int
    seconds: float


# ---- helpers -----------------------------------------------------------------


def _year(jd: float) -> float:
    return 2000.0 + (jd - 2451545.0) / 365.25


def _position(p: _n.Position) -> Position:
    return Position(p.longitude * DEG, p.latitude * DEG, p.distance,
                    p.longitude_speed * DEG, p.latitude_speed * DEG, p.distance_speed)


def _observer(lon: float, lat: float, elevation: float = 0.0) -> _n.Observer:
    return _n.Observer(lon * _RAD, lat * _RAD, elevation)


def _str(b: Optional[bytes]) -> str:
    return b.decode("utf-8") if b else ""


# ---- data files -------------------------------------------------------------

PACKAGE_DATA = ("1as", "1as-chiron", "01as", "01as-chiron", "chiron", "asteroids", "stars")
"""Data files shipped in the package (data/eph-<name>.bin)."""


def data_path(name: str) -> str:
    """Path of a packaged data file, e.g. ``data_path('stars')``."""
    if name not in PACKAGE_DATA:
        raise ValueError(f"unknown data file {name!r}; packaged: {', '.join(PACKAGE_DATA)}")
    return os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", f"eph-{name}.bin")


def load_data(source: Union[bytes, bytearray, memoryview, str, os.PathLike]) -> bool:
    """Load a data file from its bytes or a path; False if it is rejected.
    Later files override what earlier ones provide."""
    if isinstance(source, (bytes, bytearray, memoryview)):
        data = bytes(source)
        return _lib.eph_load_data(data, len(data)) == 0  # the library copies it
    return _lib.eph_load_file(os.fsencode(source)) == 0


def load_package_data(*names: str) -> None:
    """Load packaged data files: '01as', 'chiron', 'asteroids', 'stars', ..."""
    for name in names:
        if not load_data(data_path(name)):
            raise ValueError(f"cannot load data file {name!r}")


def unload_data() -> None:
    """Unload every data file (the built-in data too; see load_embedded)."""
    _lib.eph_unload_data()


def load_embedded() -> bool:
    """Reload the built-in data (the 1" series with Chiron)."""
    return _lib.eph_load_embedded() == 0


# ---- time --------------------------------------------------------------------


def julian_date(year: int, month: int, day: float) -> float:
    """Julian day of a calendar date (Gregorian; Julian before 1582-10-15);
    ``day`` carries the fraction of the day (1.5 = noon)."""
    return _lib.eph_julian_date(year, month, day)


def calendar_date(jd: float) -> CalendarDate:
    y, m, d = c_int(), c_int(), c_double()
    _lib.eph_calendar_date(jd, byref(y), byref(m), byref(d))
    return CalendarDate(y.value, m.value, d.value)


def delta_t(year: float) -> float:
    """TT - UT, seconds."""
    return _lib.eph_delta_t(year)


def ut_to_tt(jd_ut: float, year: Optional[float] = None) -> float:
    return _lib.eph_ut_to_tt(jd_ut, _year(jd_ut) if year is None else year)


def tt_to_ut(jd_tt: float, year: Optional[float] = None) -> float:
    return _lib.eph_tt_to_ut(jd_tt, _year(jd_tt) if year is None else year)


def utc_to_jd(year: Union[int, _dt.datetime], month: int = 1, day: int = 1, hour: int = 0,
              minute: int = 0, second: float = 0.0) -> Optional[JulianDates]:
    """Civil UTC (leap seconds from 1972) to TT and UT Julian days; None for an
    invalid date.  Also takes a datetime: an aware one is converted to UTC, a
    naive one is taken as UTC."""
    if isinstance(year, _dt.datetime):
        t = year
        if t.tzinfo is not None:
            t = t.astimezone(_dt.timezone.utc)
        year, month, day, hour, minute = t.year, t.month, t.day, t.hour, t.minute
        second = t.second + t.microsecond / 1e6
    tt, ut = c_double(), c_double()
    if _lib.eph_utc_to_jd(year, month, day, hour, minute, second, byref(tt), byref(ut)) != 0:
        return None
    return JulianDates(tt.value, ut.value)


def tdb_minus_tt(jd_tt: float) -> float:
    """TDB - TT, seconds."""
    return _lib.eph_tdb_minus_tt(jd_tt)


def sidereal_time(jd_ut: float) -> float:
    """Greenwich apparent sidereal time (IAU 2006), degrees."""
    return _lib.eph_sidereal_time(jd_ut) * DEG


def equation_of_time(jd_ut: float) -> Optional[float]:
    """Apparent minus mean solar time, minutes; None without the Sun's data."""
    e = c_double()
    return e.value * 1440.0 if _lib.eph_equation_of_time(jd_ut, byref(e)) == 0 else None


# ---- frames ----------------------------------------------------------------


def mean_obliquity(jd_tt: float) -> float:
    """Mean obliquity of the ecliptic (IAU 2006), degrees."""
    return _lib.eph_mean_obliquity(jd_tt) * DEG


def nutation(jd_tt: float) -> Nutation:
    """Nutation in longitude and obliquity (IAU 2000B), degrees."""
    dpsi, deps = c_double(), c_double()
    _lib.eph_nutation(jd_tt, byref(dpsi), byref(deps))
    return Nutation(dpsi.value * DEG, deps.value * DEG)


def ecliptic_to_equatorial(lon: float, lat: float, jd_tt: float) -> Equatorial:
    """Apparent ecliptic to equatorial coordinates (true obliquity), degrees."""
    ra, dec = c_double(), c_double()
    _lib.eph_ecliptic_to_equatorial(lon * _RAD, lat * _RAD, jd_tt, byref(ra), byref(dec))
    return Equatorial(ra.value * DEG, dec.value * DEG)


# ---- positions ---------------------------------------------------------------


def position(body: int, jd_ut: float, apparent: bool = True,
             year: Optional[float] = None) -> Optional[Position]:
    """Geocentric ecliptic position of date for a UT Julian day, apparent by
    default; None without data.  ``year`` (fractional) selects ΔT."""
    out = _n.Position()
    fn = _lib.eph_position_apparent_ut if apparent else _lib.eph_position_ut
    rc = fn(body, jd_ut, _year(jd_ut) if year is None else year, byref(out))
    return _position(out) if rc == 0 else None


def position_tt(body: int, jd_tt: float, apparent: bool = True) -> Optional[Position]:
    out = _n.Position()
    fn = _lib.eph_position_apparent_tt if apparent else _lib.eph_position_tt
    return _position(out) if fn(body, jd_tt, byref(out)) == 0 else None


def heliocentric(body: int, jd_tt: float) -> Optional[Heliocentric]:
    """Heliocentric ecliptic of date; None for the Sun, the Moon or without data."""
    lon, lat, dist = c_double(), c_double(), c_double()
    if _lib.eph_heliocentric_tt(body, jd_tt, byref(lon), byref(lat), byref(dist)) != 0:
        return None
    return Heliocentric(lon.value * DEG, lat.value * DEG, dist.value)


def body_name(body: int) -> str:
    """'Sun', 'Moon', ... (minor planets: their names)."""
    if body >= MINOR_BASE:
        return minor_name(body)
    return _str(_lib.eph_body_name(body))


def topocentric(body: int, jd_ut: float, lon: float, lat: float, elevation: float = 0.0,
                year: Optional[float] = None) -> Optional[Position]:
    """Apparent position seen from a place (east longitude, latitude in degrees,
    elevation in metres)."""
    out = _n.Position()
    obs = _observer(lon, lat, elevation)
    rc = _lib.eph_position_topocentric_ut(body, jd_ut, _year(jd_ut) if year is None else year,
                                          byref(obs), byref(out))
    return _position(out) if rc == 0 else None


# ---- houses ------------------------------------------------------------------


def asc_mc(jd_ut: float, lon: float, lat: float) -> AscMc:
    asc, mc = c_double(), c_double()
    obs = _observer(lon, lat)
    _lib.eph_asc_mc(jd_ut, byref(obs), byref(asc), byref(mc))
    return AscMc(asc.value * DEG, mc.value * DEG)


def houses_info(jd_ut: float, lon: float, lat: float,
                system: int = HouseSystem.PORPHYRY) -> Optional[HousesInfo]:
    """The twelve cusps, and whether the system is undefined at this latitude
    (Placidus, Koch, Alcabitius inside the polar circles: Porphyry is used)."""
    cusps = (c_double * 12)()
    obs = _observer(lon, lat)
    rc = _lib.eph_houses(jd_ut, byref(obs), system, cusps)
    if rc < 0:
        return None
    return HousesInfo([c * DEG for c in cusps], rc == 1)


def houses(jd_ut: float, lon: float, lat: float,
           system: int = HouseSystem.PORPHYRY) -> Optional[List[float]]:
    """Twelve cusps (degrees); None on bad arguments."""
    info = houses_info(jd_ut, lon, lat, system)
    return info.cusps if info else None


def house_angles(jd_ut: float, lon: float, lat: float) -> HouseAngles:
    """Ascendant, MC, ARMC, Vertex, equatorial Ascendant, co-Ascendants and
    polar Ascendant, degrees."""
    out = _n.Angles()
    obs = _observer(lon, lat)
    _lib.eph_house_angles(jd_ut, byref(obs), byref(out))
    return HouseAngles(*(getattr(out, f) * DEG for f, _ in _n.Angles._fields_))


def house_position(cusps: List[float], lon: float) -> float:
    """House 1.0 .. 12.999 of a longitude, given cusps from houses()."""
    arr = (c_double * 12)(*(c * _RAD for c in cusps))
    return _lib.eph_house_position(arr, lon * _RAD)


# ---- phenomena -----------------------------------------------------------------


def phenomena(body: int, jd_tt: float) -> Optional[Phenomena]:
    """Phase angle, illuminated fraction, elongation, apparent diameter
    (degrees) and visual magnitude; None without data."""
    out = _n.Phenomena()
    if _lib.eph_phenomena_tt(body, jd_tt, byref(out)) != 0:
        return None
    return Phenomena(out.phase_angle * DEG, out.phase, out.elongation * DEG,
                     out.diameter * DEG, out.magnitude)


# ---- searches --------------------------------------------------------------------

AyanamsaOption = Union[int, CustomAyanamsa]


def longitude_crossing(body: int, lon: float, jd_ut: float, direction: int = 1,
                       max_days: float = 400, ayanamsa: Optional[AyanamsaOption] = None
                       ) -> Optional[float]:
    """When the apparent longitude reaches ``lon`` (returns, ingresses); UT in
    and out, None if not within max_days.  With ``ayanamsa`` (an Ayanamsa or a
    CustomAyanamsa) the longitude is sidereal, the ayanamsa taken at the
    crossing itself (sidereal returns and ingresses)."""
    out = c_double()
    if ayanamsa is None:
        rc = _lib.eph_longitude_crossing_ut(body, lon * _RAD, jd_ut, direction, max_days, byref(out))
        return out.value if rc == 0 else None
    jd_tt = ut_to_tt(jd_ut)
    if isinstance(ayanamsa, CustomAyanamsa):
        rc = _lib.eph_longitude_crossing_sidereal_custom(
            body, lon * _RAD, ayanamsa.jd_t0, ayanamsa.value_at_t0 * _RAD, jd_tt, direction,
            max_days, byref(out))
    else:
        rc = _lib.eph_longitude_crossing_sidereal(body, lon * _RAD, int(ayanamsa), jd_tt, direction,
                                                  max_days, byref(out))
    return tt_to_ut(out.value) if rc == 0 else None


def heliocentric_crossing(body: int, lon: float, jd_tt: float, direction: int = 1,
                          max_days: float = 40000, ayanamsa: Optional[AyanamsaOption] = None
                          ) -> Optional[float]:
    """When the heliocentric longitude reaches ``lon``; TT in and out."""
    out = c_double()
    if ayanamsa is None:
        rc = _lib.eph_heliocentric_crossing(body, lon * _RAD, jd_tt, direction, max_days, byref(out))
    elif isinstance(ayanamsa, CustomAyanamsa):
        rc = _lib.eph_heliocentric_crossing_sidereal_custom(
            body, lon * _RAD, ayanamsa.jd_t0, ayanamsa.value_at_t0 * _RAD, jd_tt, direction,
            max_days, byref(out))
    else:
        rc = _lib.eph_heliocentric_crossing_sidereal(body, lon * _RAD, int(ayanamsa), jd_tt,
                                                     direction, max_days, byref(out))
    return out.value if rc == 0 else None


def node_crossing(body: int, jd_tt: float, direction: int = 1,
                  max_days: float = 40) -> Optional[NodeCrossing]:
    """When the body crosses the ecliptic; TT in and out."""
    jd, lon = c_double(), c_double()
    rc = _lib.eph_node_crossing(body, jd_tt, direction, max_days, byref(jd), byref(lon))
    if rc < 1:
        return None
    return NodeCrossing(jd.value, "ascending" if rc == 1 else "descending", lon.value * DEG)


def station(body: int, jd_tt: float, direction: int = 1,
            max_days: float = 800) -> Optional[Station]:
    """When the body turns retrograde or direct; TT in and out."""
    jd = c_double()
    rc = _lib.eph_station(body, jd_tt, direction, max_days, byref(jd))
    if rc < 1:
        return None
    return Station(jd.value, "retrograde" if rc == 1 else "direct")


# ---- planetary nodes and apsides --------------------------------------------------


def orbit(body: int, jd_tt: float, method: int = OrbitMethod.MEAN) -> Optional[Orbit]:
    """The orbit's elements (degrees, AU) and its nodes, perihelion and
    aphelion, heliocentric and geocentric; None if unsupported or without data."""
    out = _n.Orbit()
    if _lib.eph_orbit_tt(body, jd_tt, method, byref(out)) != 0:
        return None

    def points(arr) -> OrbitPoints:
        return OrbitPoints(*(_position(arr[i]) for i in range(4)))

    return Orbit(out.semi_major_axis, out.eccentricity, out.inclination * DEG, out.node * DEG,
                 out.perihelion * DEG, points(out.heliocentric), points(out.geocentric))


# ---- eclipses ------------------------------------------------------------------------


def solar_eclipse(jd_tt: float, direction: int = 1) -> Optional[SolarEclipse]:
    """The next (or previous, direction -1) solar eclipse anywhere; TT."""
    o = _n.SolarEclipse()
    if _lib.eph_solar_eclipse_next(jd_tt, direction, byref(o)) != 0:
        return None
    return SolarEclipse(EclipseType(o.type), o.maximum, o.gamma, o.magnitude, o.latitude * DEG,
                        o.longitude * DEG, o.saros, o.partial_begin, o.partial_end,
                        o.central_begin, o.central_end, o.center_line_begin, o.center_line_end)


def solar_eclipse_local(jd_tt: float, lon: float, lat: float, elevation: float = 0.0,
                        direction: int = 1) -> Optional[LocalSolarEclipse]:
    """The next solar eclipse seen from a place (check ``altitude``: the Sun
    may be below the horizon)."""
    o = _n.LocalSolarEclipse()
    obs = _observer(lon, lat, elevation)
    if _lib.eph_solar_eclipse_local(jd_tt, byref(obs), direction, byref(o)) != 0:
        return None
    return LocalSolarEclipse(EclipseType(o.type), o.maximum, o.first_contact, o.second_contact,
                             o.third_contact, o.fourth_contact, o.magnitude, o.obscuration,
                             o.altitude * DEG, o.saros)


def solar_eclipse_where(jd_tt: float) -> Optional[GeoPoint]:
    """Where the shadow axis meets the Earth; None if it misses."""
    lat, lon = c_double(), c_double()
    if _lib.eph_solar_eclipse_where(jd_tt, byref(lat), byref(lon)) != 0:
        return None
    return GeoPoint(lat.value * DEG, lon.value * DEG)


def solar_eclipse_how(jd_tt: float, lon: float, lat: float,
                      elevation: float = 0.0) -> Optional[LocalEclipseState]:
    """The solar eclipse seen from a place at an instant (type 0: none)."""
    o = _n.LocalEclipseState()
    obs = _observer(lon, lat, elevation)
    if _lib.eph_solar_eclipse_how(jd_tt, byref(obs), byref(o)) != 0:
        return None
    return LocalEclipseState(EclipseType(o.type), o.magnitude, o.obscuration, o.altitude * DEG)


def lunar_eclipse(jd_tt: float, direction: int = 1) -> Optional[LunarEclipse]:
    """The next (or previous) lunar eclipse, penumbral ones included; TT."""
    o = _n.LunarEclipse()
    if _lib.eph_lunar_eclipse_next(jd_tt, direction, byref(o)) != 0:
        return None
    return LunarEclipse(EclipseType(o.type), o.maximum, o.gamma, o.umbral_magnitude,
                        o.penumbral_magnitude, o.saros, o.penumbral_begin, o.penumbral_end,
                        o.partial_begin, o.partial_end, o.total_begin, o.total_end)


# ---- minor planets ------------------------------------------------------------------


def minor_body(name: Union[str, int]) -> Optional[int]:
    """The body id of a loaded minor planet: 433, "433", "Eros" or "433 Eros";
    None if none is loaded.  Works with position(), searches, phenomena(),
    orbit() and heliocentric()."""
    i = _lib.eph_minor_find(str(name).encode("utf-8"))
    return None if i < 0 else MINOR_BASE + i


def minor_count() -> int:
    return _lib.eph_minor_count()


def minor_name(body: int) -> str:
    return _str(_lib.eph_minor_name(body - MINOR_BASE))


def minor_list() -> List[MinorPlanet]:
    return [MinorPlanet(MINOR_BASE + i, _lib.eph_minor_number(i), _str(_lib.eph_minor_name(i)))
            for i in range(_lib.eph_minor_count())]


# ---- fixed stars ------------------------------------------------------------------


def star_find(name: str) -> Optional[int]:
    """Star index by name ("Regulus"), designation ("alf Leo") or "HIP 49669";
    None if unknown (or the stars are not loaded)."""
    i = _lib.eph_star_find(name.encode("utf-8"))
    return None if i < 0 else i


def star_count() -> int:
    """Stars loaded (0 without eph-stars.bin); indices 0 .. count - 1."""
    return _lib.eph_star_count()


def star_name(index: int) -> str:
    return _str(_lib.eph_star_name(index))


def star_designation(index: int) -> str:
    return _str(_lib.eph_star_designation(index))


def star_magnitude(index: int) -> float:
    return _lib.eph_star_magnitude(index)


def star_position(index: int, jd_ut: float, year: Optional[float] = None) -> Optional[Position]:
    out = _n.Position()
    rc = _lib.eph_star_position_apparent_ut(index, jd_ut, _year(jd_ut) if year is None else year,
                                            byref(out))
    return _position(out) if rc == 0 else None


# ---- ayanamsas --------------------------------------------------------------------


def ayanamsa(system: int, jd_tt: float) -> float:
    """Ayanamsa in degrees: sidereal = tropical - ayanamsa."""
    return _lib.eph_ayanamsa(system, jd_tt) * DEG


def ayanamsa_custom(jd_t0: float, value_at_t0: float, jd_tt: float) -> float:
    """A user-defined ayanamsa (degrees at jd_t0) carried by the precession."""
    return _lib.eph_ayanamsa_custom(jd_t0, value_at_t0 * _RAD, jd_tt) * DEG


def ayanamsa_lahiri(jd_tt: float) -> float:
    return _lib.eph_ayanamsa_lahiri(jd_tt) * DEG


def ayanamsa_name(system: int) -> str:
    """'Lahiri', 'Fagan/Bradley', ..."""
    return _str(_lib.eph_ayanamsa_name(system))


# ---- utilities ---------------------------------------------------------------------


def midpoint(a: float, b: float) -> float:
    """Midpoint of the shorter arc between two longitudes, degrees."""
    return _lib.eph_midpoint(a * _RAD, b * _RAD) * DEG


def split_longitude(lon: float) -> SplitLongitude:
    """Zodiac sign (0 = Aries .. 11) and the position within it."""
    sign, deg, minute, sec = c_int(), c_int(), c_int(), c_double()
    _lib.eph_split_longitude(lon * _RAD, byref(sign), byref(deg), byref(minute), byref(sec))
    return SplitLongitude(sign.value, deg.value, minute.value, sec.value)

