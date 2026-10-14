<p align="center"><img alt="Ephemera" src="https://raw.githubusercontent.com/intellecat/ephemera/main/assets/logo/lockup-horizontal-light.svg" width="420"></p>

# pyephemera

**A pocket-sized, MIT-licensed ephemeris for astrology and astronomy.**
JPL-grade planetary positions, houses, ayanamsas, eclipses, fixed stars and
1,000+ asteroids, from a small C library with no dependencies. Wheels for
Linux, macOS and Windows; Python 3.9+.

- **MIT licensed:** use it in closed-source and commercial apps.
- **Accurate:** within 0.6″ of NASA JPL's DE440 over 1550–2650 with the
  built-in data, or 0.03″ with the high-precision tables.
- **Small:** the planets, the Moon and Chiron are built in (204 KB); more
  data files come with the package.
- **Complete for astrology:** 14 house systems, 42 ayanamsas, nodes and
  Lilith, Uranian points, eclipses, returns and stations, sidereal
  returns, planetary nodes and apsides, and 411 fixed stars.
- **Independently verified** against JPL, NASA's eclipse canon and
  India's official ephemeris. See the
  [accuracy page](https://github.com/intellecat/ephemera/blob/main/docs/accuracy.md).

```sh
pip install pyephemera
```

```python
import ephemera as eph

jd = eph.julian_date(2026, 1, 1.5)                    # 1 Jan 2026, 12:00 UT
moon = eph.position(eph.Body.MOON, jd)                # degrees, AU, per day
cusps = eph.houses(jd, -0.1278, 51.5074, eph.HouseSystem.PLACIDUS)  # London
lahiri = eph.ayanamsa(eph.Ayanamsa.LAHIRI, eph.ut_to_tt(jd))        # sidereal = tropical − ayanamsa
eclipse = eph.solar_eclipse(eph.ut_to_tt(jd))         # the next solar eclipse
```

## Times

Positions take a **UT** Julian day; searches, eclipses and orbits take
**TT**. `utc_to_jd` converts civil time (with leap seconds) and accepts a
`datetime`:

```python
from datetime import datetime
from zoneinfo import ZoneInfo

jd_tt, jd_ut = eph.utc_to_jd(datetime(1990, 6, 15, 14, 30, tzinfo=ZoneInfo("Europe/London")))
```

## Data files

The 1″ series with Chiron is built in. The package also carries:

```python
eph.load_package_data("stars")       # 411 named stars, star-based ayanamsas
eph.load_package_data("asteroids")   # Ceres, Pallas, Juno, Vesta (1600–2500)
eph.load_package_data("01as")        # the 0.2″ series
eph.load_data("de440-2000.bin")      # any file from the repository, by path or bytes
```

The DE440 tables (0.002″) and the minor-planet files (1,000+ asteroids,
1900–2100) are downloads from the
[repository](https://github.com/intellecat/ephemera/blob/main/docs/data.md).
Loaded data is shared by the whole process.

## Documentation

The [guide](https://github.com/intellecat/ephemera/blob/main/docs/guide.md)
covers every feature. The Python names follow the JavaScript API in
snake_case: `houses_info`, `longitude_crossing`, `solar_eclipse_local`,
`minor_body`; calls that find nothing return `None`.

## License

MIT. Built from public data: NASA JPL, NASA's eclipse catalogues and
Hipparcos. No third-party code is included.
