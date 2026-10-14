"""Smoke tests of the Python package (make python-test).  Reference values
come from the C library's tests and the npm package's."""

import datetime as dt
import math

import pytest

import ephemera as eph

JD = eph.julian_date(2026, 1, 1.5)  # 1 Jan 2026, 12:00 UT
LONDON = (-0.1278, 51.5074)


@pytest.fixture(autouse=True)
def built_in_data():
    """Each test starts with only the built-in data loaded."""
    eph.unload_data()
    assert eph.load_embedded()
    yield


def test_time():
    assert JD == 2461042.0
    assert eph.calendar_date(JD) == (2026, 1, 1.5)
    tt, ut = eph.utc_to_jd(2026, 1, 1, 12)
    assert abs((tt - ut) * 86400 - eph.delta_t(2026.0)) < 1.0
    aware = dt.datetime(2026, 1, 1, 13, tzinfo=dt.timezone(dt.timedelta(hours=1)))
    assert eph.utc_to_jd(aware) == eph.utc_to_jd(2026, 1, 1, 12)
    assert eph.utc_to_jd(2026, 2, 30) is None


def test_built_in_positions():
    moon = eph.position(eph.Body.MOON, JD)
    assert abs(moon.longitude - 74.2252) < 1e-3  # the C library, eph-1as.bin
    assert 11 < moon.longitude_speed < 16
    for body in (eph.Body.SUN, eph.Body.MERCURY, eph.Body.PLUTO, eph.Body.CHIRON,
                 eph.Body.TRUE_NODE, eph.Body.MEAN_LILITH, eph.Body.CUPIDO):
        assert eph.position(body, JD) is not None, body
    assert eph.position(eph.Body.CERES, JD) is None  # needs eph-asteroids.bin
    assert eph.heliocentric(eph.Body.EARTH, eph.ut_to_tt(JD)).distance == pytest.approx(0.983, abs=0.001)
    assert eph.body_name(eph.Body.MOON) == "Moon"


def test_houses_and_ayanamsa():
    cusps = eph.houses(JD, *LONDON, eph.HouseSystem.PLACIDUS)
    assert len(cusps) == 12
    assert cusps[0] == pytest.approx(25.536293992570464, abs=1e-9)  # the npm package
    assert cusps[9] == pytest.approx(280.137101405305, abs=1e-9)
    assert eph.asc_mc(JD, *LONDON).asc == pytest.approx(cusps[0], abs=1e-9)
    assert eph.houses_info(JD, 18.96, 69.65, eph.HouseSystem.PLACIDUS) is not None
    assert eph.house_position(cusps, cusps[0] + 1e-9) == pytest.approx(1.0, abs=1e-6)
    assert eph.ayanamsa(eph.Ayanamsa.LAHIRI, 2451545.0) == pytest.approx(23.857092325, abs=1e-6)
    assert eph.ayanamsa_name(eph.Ayanamsa.FAGAN_BRADLEY) == "Fagan/Bradley"


def test_searches_and_eclipses():
    sun = eph.position(eph.Body.SUN, JD).longitude
    ret = eph.longitude_crossing(eph.Body.SUN, sun, JD + 10)
    assert 365 < ret - JD < 366
    sidereal = eph.longitude_crossing(eph.Body.SUN, 0.0, JD, ayanamsa=eph.Ayanamsa.LAHIRI)
    assert sidereal is not None
    assert eph.station(eph.Body.MERCURY, eph.ut_to_tt(JD)).turns in ("retrograde", "direct")
    eclipse = eph.solar_eclipse(2460400.0)
    assert eclipse.saros == 139 and eclipse.type & eph.EclipseType.TOTAL
    assert eph.lunar_eclipse(2460400.0) is not None
    orbit = eph.orbit(eph.Body.MARS, eph.ut_to_tt(JD))
    assert orbit.semi_major_axis == pytest.approx(1.5237, abs=1e-3)


def test_package_data():
    assert eph.star_count() == 0 and eph.star_find("Regulus") is None
    eph.load_package_data("stars", "asteroids")
    assert eph.position(eph.Body.CERES, JD) is not None
    n = eph.star_count()
    assert n == 411
    assert all(eph.star_find(eph.star_name(i)) == i for i in range(n))
    reg = eph.star_find("Regulus")
    assert eph.star_designation(reg) == "alf Leo"
    assert abs(eph.star_magnitude(reg) - 1.4) < 0.1
    assert eph.star_position(reg, JD).longitude == pytest.approx(150.1, abs=0.2)
    with pytest.raises(ValueError):
        eph.data_path("nope")


def test_utilities():
    assert eph.midpoint(350, 10) == pytest.approx(0.0, abs=1e-9) or eph.midpoint(350, 10) == pytest.approx(360.0)
    s = eph.split_longitude(123.5)
    assert (s.sign, s.degrees, s.minutes) == (4, 3, 30)
    assert math.isfinite(eph.sidereal_time(JD))
