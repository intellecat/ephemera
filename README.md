# ephemeris

A small, dependency-free C/C++ library that computes positions of the Sun,
Moon, Pluto and the eight planets from pure analytical series.  No JPL data
file, no SPICE kernel, no 30 MB download at runtime — every coefficient is
compiled straight into the library.

Internally it evaluates:

* **VSOP87D** (Bretagnon & Francou, 1988) for Mercury…Neptune — the full
  IMCCE coefficient tables (31,577 terms) for heliocentric spherical
  coordinates in the *mean ecliptic and equinox of date*.
* **ELP2000-82B** (Chapront-Touze & Chapront) full lunar theory for the
  Moon; ~0.1" near J2000, growing to a few arcseconds at the range edges
  because its constants are fitted to DE200/LE200.
* A compact Poisson series for Pluto fitted to JPL DE440s over 1850-2140,
  and small Chebyshev DE440 corrections for Uranus and Neptune.
* **IAU 1980** nutation (63 terms), the IAU 1980 mean obliquity, and the
  Espenak–Meeus ΔT polynomials for the UT↔TT conversion.

Angles are radians, distances are AU, and time is Julian Date (TT or UT).

## Accuracy

Measured against JPL DE440s over 1850–2140 (see `tests/`), maximum error in
geocentric ecliptic longitude/latitude.  The reference is reduced to the
library's frame (mean ecliptic/equinox of date, IAU 1976), so these numbers
reflect the theories, not a frame-model difference:

| Body    | geometric | apparent | notes                          |
|---------|-----------|----------|--------------------------------|
| Mercury | 0.12″     | 0.13″    |                                |
| Venus   | 0.14″     | 0.21″    |                                |
| Mars    | 0.09″     | 0.10″    |                                |
| Jupiter | 0.49″     | 0.49″    | DE barycenter reference        |
| Saturn  | 0.37″     | 0.35″    | DE barycenter reference        |
| Uranus  | 0.01″     | 0.05″    | VSOP87 + DE440 correction      |
| Neptune | 0.01″     | 0.05″    | VSOP87 + DE440 correction      |
| Pluto   | 0.04″     | 0.07″    | Poisson fit to DE440s          |
| Sun     | 0.13″     | 0.13″    |                                |
| Moon    | 0.20″     | 0.94″    | full ELP2000-82B + DE440 corr  |

Inner planets are one to three arcseconds as advertised; the outer planets
are better than that over the fitted span, Pluto is sub-arcsecond over 1850–2140, and the Moon is ~0.1" near
J2000 (up to a few arcseconds at the edges, reflecting its DE200 fit).

### Where the remaining error comes from

* **Uranus/Neptune.**  VSOP87 is fitted to DE200 and drifts from DE440 by
  ~1-2" over 1850-2140.  This is absorbed by a small Chebyshev correction
  (`tools/fit_correction.py`), so both are now ~0.01-0.05".  The
  correction is clamped to its 1850-2140 fit window; outside it freezes at
  the boundary value.
* **Moon.**  ELP2000-82B is fitted to DE200/LE200 and drifts ~1.2"/
  century from DE440; this is absorbed by a Chebyshev correction
  (`tools/fit_moon_correction.py`), leaving ~0.2" at the 1850/2140 edges
  and ~0.01" near J2000.  The residual apparent figure (~0.9") is set by
  the first-order nutation handling, not the lunar theory.
* **Apparent place.**  Uses IAU-1976 precession, IAU-1980 nutation, the
  IAU-1980 obliquity, iterated light-time and first-order annual
  aberration.  Relativistic light deflection (<0.02" except within a few
  degrees of the Sun) and higher-order aberration are omitted.  The
  frames deliberately follow VSOP87D (IAU-1976); migrating the whole
  library to IAU-2000/2006 is possible but would require transforming the
  VSOP87 series as well.  `eph_tdb_minus_tt()` is provided for callers
  that want the TDB timescale used by the dynamical theories.

## Building

With CMake:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build            # run the DE440s validation
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
JPL DE440s kernel with Skyfield (reduced with ERFA/IAU-1976); the kernel is **not** part of the library
and is never shipped.  See `tools/make_reference.py` and
`tools/gen_reference.py` if you want to regenerate it.

## Layout

```
include/ephemeris/ephemeris.h   public C API
src/                            VSOP87 evaluator, Moon, Pluto, frames, time
data/vsop87d_data.{h,cpp}       generated coefficient tables
tools/                          coefficient and reference generators
tests/                          DE440s validation + time/calendar checks
examples/positions.c            CLI that prints a chart
```

## References

* P. Bretagnon, G. Francou, *Astron. Astrophys.* **202**, 309 (1988).
* J. Meeus, *Astronomical Algorithms*, 2nd ed. (Willmann-Bell, 1998).
* F. Espenak, J. Meeus, *Five Millennium Canon of Solar Eclipses*, NASA/TP-2006-214141.
