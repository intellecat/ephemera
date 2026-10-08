# ephemeris

A small, dependency-free C/C++ library that computes positions of the Sun,
Moon, Pluto and the eight planets from pure analytical series.  No JPL data
file, no SPICE kernel, no 30 MB download at runtime — every coefficient is
compiled straight into the library.

Internally it evaluates:

* **VSOP87D** (Bretagnon & Francou, 1988) for Mercury…Neptune — the full
  IMCCE coefficient tables (31,577 terms) for heliocentric spherical
  coordinates in the *mean ecliptic and equinox of date*.
* **Meeus / ELP2000** truncated lunar theory (Astronomical Algorithms,
  ch. 47) for the Moon, good to a few arcseconds.
* **Meeus** analytic Pluto series (ch. 37) for Pluto.
* **IAU 1980** nutation (63 terms), the IAU 1980 mean obliquity, and the
  Espenak–Meeus ΔT polynomials for the UT↔TT conversion.

Angles are radians, distances are AU, and time is Julian Date (TT or UT).

## Accuracy

Measured against JPL DE421 over 1900–2050 (see `tests/`), maximum error in
geocentric ecliptic longitude/latitude:

| Body    | geometric | apparent | notes                         |
|---------|-----------|----------|-------------------------------|
| Mercury | 0.21″     | 0.20″    |                               |
| Venus   | 0.21″     | 0.28″    |                               |
| Mars    | 0.21″     | 0.29″    |                               |
| Jupiter | 0.57″     | 0.57″    | DE barycenter reference       |
| Saturn  | 0.44″     | 0.41″    | DE barycenter reference       |
| Uranus  | 1.16″     | 1.18″    | DE barycenter reference       |
| Neptune | 2.45″     | 2.43″    | DE barycenter reference       |
| Pluto   | 53.9″     | 53.9″    | Meeus fit to DE200, ~1 arcmin |
| Sun     | 0.21″     | 0.21″    |                               |
| Moon    | 4.2″      | 4.9″     | Meeus truncated ELP2000       |

Inner planets are one to three arcseconds as advertised; the outer planets
are better than that for 1900–2050, and Pluto stays within about one
arcminute.

## Building

With CMake:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build            # run the DE421 validation
./build/ephcliapp 2026 1 1 12     # year month day [UT fraction]
```

Or with the convenience Makefile:

```sh
make -j && make test
./build/ephcliapp 2000 1 1.5
```

Link against `libephemeris.a` and include `ephemeris/ephemeris.h`.  The
public API is plain C, so it can be used from C or C++.

## Using the API

```c
#include "ephemeris/ephemeris.h"

double jd_ut = eph_julian_date(1992, 4, 12.0);   /* 1992-04-12 00:00 UT */
double year  = 1992.25;                          /* for Delta-T         */

eph_position moon;
eph_position_apparent_ut(EPH_MOON, jd_ut, year, &moon);

printf("Moon ecliptic longitude: %.6f deg\n", moon.longitude * 180.0 / M_PI);
printf("Moon daily motion:       %.6f deg/day\n",
       moon.longitude_speed * 180.0 / M_PI);
```

`eph_position_tt` / `eph_position_apparent_tt` take Terrestrial Time
directly.  The `_ut` variants apply ΔT.  Every position carries its own
longitude/latitude/distance and their time derivatives (negative longitude
speed means retrograde).

For real equatorial coordinates:

```c
double ra, dec;
eph_ecliptic_to_equatorial(moon.longitude, moon.latitude, jd_tt, &ra, &dec);
```

Heliocentric coordinates are available for the planets and Pluto via
`eph_heliocentric_tt`.

## Conventions

* VSOP87D is referred to the **mean dynamical ecliptic and equinox of
  date**, which is exactly the tropical frame used in astrology.  The
  `_apparent_` functions additionally apply light-time (iterated), annual
  aberration from the Earth's velocity, and nutation in longitude.
* The Moon uses the mean-equinox-of-date series; nutation is applied for
  apparent positions.
* The `year` argument to the `_ut` functions is only used to evaluate the
  ΔT polynomial; pass a fractional calendar year.

## Coefficient provenance and regeneration

The VSOP87D tables are downloaded from IMCCE (Observatoire de Paris) and
compiled into `data/vsop87d_data.cpp`:

```sh
python3 tools/gen_vsop87.py        # downloads data/vsop87d/* and regenerates the tables
```

The test reference file `tests/reference_data.h` is generated once from a
JPL DE421 kernel with Skyfield; the kernel is **not** part of the library
and is never shipped.  See `tools/make_reference.py` and
`tools/gen_reference.py` if you want to regenerate it.

## Layout

```
include/ephemeris/ephemeris.h   public C API
src/                            VSOP87 evaluator, Moon, Pluto, frames, time
data/vsop87d_data.{h,cpp}       generated coefficient tables
tools/                          coefficient and reference generators
tests/                          DE421 validation + time/calendar checks
examples/positions.c            CLI that prints a chart
```

## References

* P. Bretagnon, G. Francou, *Astron. Astrophys.* **202**, 309 (1988).
* J. Meeus, *Astronomical Algorithms*, 2nd ed. (Willmann-Bell, 1998).
* F. Espenak, J. Meeus, *Five Millennium Canon of Solar Eclipses*, NASA/TP-2006-214141.
