# ephemeris

A small, dependency-free C/C++ library that computes positions of the Sun,
Moon, Pluto, the eight planets and the four largest asteroids (Ceres,
Pallas, Juno, Vesta) from pure analytical series.  No JPL data
file, no SPICE kernel, no 30 MB download at runtime — every coefficient is
compiled straight into the library.

Internally it evaluates:

* **VSOP87D** (Bretagnon & Francou, 1988) for Mercury…Neptune — the full
  IMCCE coefficient tables (31,577 terms) for heliocentric spherical
  coordinates in the *mean ecliptic and equinox of date*.
* **ELP2000-82B** (Chapront-Touze & Chapront) full lunar theory for the
  Moon; ~0.1" near J2000, growing to a few arcseconds at the range edges
  because its constants are fitted to DE200/LE200.
* Compact Chebyshev series for Ceres/Pallas/Juno/Vesta fitted to JPL
  Horizons (public data), heliocentric fit error <0.02" over 1850-2140.
  Apparent positions from UT match Horizons to ~0.04".
* A compact Poisson series for Pluto fitted to JPL DE440s over 1850-2140,
  and small Chebyshev DE440 corrections for Uranus and Neptune.
* **IAU 2000B** nutation (77 terms), a rotation from the series' native
  IAU-1976 frame to the **IAU-2006** ecliptic of date (the frame Swiss
  Ephemeris and JPL Horizons use), and a yearly ΔT table for UT↔TT
  (IERS observations; Stephenson, Morrison & Hohenkerk 2016 before 1973
  and as the forecast to 2150; Espenak–Meeus polynomials outside
  1620–2150).

Angles are radians, distances are AU, and time is Julian Date (TT or UT).

## Accuracy

Measured against JPL DE440s over 1850–2140 (see `tests/`), maximum error in
geocentric ecliptic longitude/latitude.  The reference is DE440s via
Skyfield, reduced with ERFA to the IAU-2006/2000A ecliptic of date -- the
same frame Swiss Ephemeris and Horizons report -- and includes dates
within ~1° of the Sun:

| Body    | geometric | apparent | compressed | notes                         |
|---------|-----------|----------|------------|-------------------------------|
| Mercury | 0.05″     | 0.05″    | 0.002″     | VSOP87 + DE440 correction     |
| Venus   | 0.06″     | 0.06″    | 0.002″     | VSOP87 + DE440 correction     |
| Mars    | 0.05″     | 0.05″    | 0.002″     | VSOP87 + DE440 correction     |
| Jupiter | 0.07″     | 0.08″    | 0.002″     | VSOP87 + DE440 correction     |
| Saturn  | 0.02″     | 0.03″    | 0.003″     | VSOP87 + DE440 correction     |
| Uranus  | 0.01″     | 0.02″    | 0.002″     | VSOP87 + DE440 correction     |
| Neptune | 0.01″     | 0.01″    | 0.002″     | VSOP87 + DE440 correction     |
| Pluto   | 0.04″     | 0.04″    | 0.002″     | Poisson fit to DE440s         |
| Sun     | 0.05″     | 0.05″    | 0.002″     | VSOP87 + DE440 correction     |
| Moon    | 0.20″     | 0.20″    | 0.002″     | full ELP2000-82B + DE440 corr |

Every planet, the Sun and Pluto are sub-0.1" over 1850–2140; the Moon is
~0.05" typical and ~0.2" at the range edges.  On 3000 random dates the
analytic median is 0.003–0.05"; with the optional compressed backend it is
**~0.001" (max 0.004")**, on par with Swiss Ephemeris' .se1 files
(0.0005–0.0007" median, 0.002" for the Moon).

### Where the remaining error comes from

* **Uranus/Neptune.**  VSOP87 is fitted to DE200 and drifts from DE440 by
  ~1-2" over 1850-2140.  This is absorbed by a small Chebyshev correction
  (`tools/fit_correction.py`), so both are now ~0.01-0.05".  The
  correction is clamped to its 1850-2140 fit window; outside it freezes at
  the boundary value.
* **Moon.**  ELP2000-82B is fitted to DE200/LE200 and drifts ~1.2"/
  century from DE440; this is absorbed by a Chebyshev correction
  (`tools/fit_moon_correction.py`), leaving ~0.2" at the 1850/2140 edges
  and ~0.01" near J2000.  The Moon's apparent place applies light-time
  (~0.7") and nutation but no annual aberration, matching the reference to
  ~0.001".
* **Frame.**  The series are evaluated in their native IAU-1976 frame and
  rotated to the IAU-2006 ecliptic of date (`frame_1976_to_2006`, checked
  against ERFA to <0.001 mas).  Without it longitudes drift by 0.3"/century
  and the ecliptic is tilted by 0.04".  Nutation is IAU 2000B, within
  ~1 mas of IAU 2000A near J2000 and ~3 mas at the 1850/2140 ends (IAU
  1980 was off by up to 18 mas).
* **Apparent place.**  Iterated light-time, relativistic light deflection
  by the Sun (finite-distance form: up to ~1.7" for a body behind the Sun,
  ~0 for one in front of it), first-order annual aberration and nutation.
  The Moon uses light-time without the annual-aberration term (the two
  cancel to first order for a body moving with the Earth).
  `eph_tdb_minus_tt()` is provided for callers that want the TDB timescale
  used by the dynamical theories.

## Performance

Each position evaluates the body's series three times (position, velocity
and acceleration together); light-time uses that expansion instead of
re-evaluating, speeds come from the same velocities, and the Earth is
cached per date and thread.  A full chart (14 bodies with speeds, Asc/MC,
houses) takes ~1.4 ms with the analytical series and ~0.11 ms with the
compressed table, on one core of a 2024 desktop CPU.

## Building

With CMake:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build            # run the DE440s validation
./build/ephcliapp 2026 1 1 12     # year month day [UT hours]
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

For an observer at a specific location, apply diurnal parallax:

```c
eph_observer obs = { lon_rad, lat_rad, elevation_m };  /* east-positive */
eph_position p;
eph_position_topocentric_ut(EPH_MOON, jd_ut, year, &obs, &p);
```

Validated to ~0.3" vs Skyfield/DE440s (1987-2005).  Uses apparent
sidereal time.  The Moon's ~1 degree parallax is the dominant effect.  Refraction
is not applied (it affects altitude, not the ecliptic longitude used in
charts).

## Sidereal zodiac and houses

```c
double ayan = eph_ayanamsa_lahiri(jd_tt);      /* sidereal = tropical - ayan */
double asc, mc;
eph_asc_mc(jd_ut, &obs, &asc, &mc);
double cusps[12];
eph_houses(jd_ut, &obs, EPH_HOUSES_PORPHYRY, cusps);
```

House systems: Whole Sign, Equal, Porphyry (Placidus/Koch are future work).
Asc/MC use apparent sidereal time and match Swiss Ephemeris to <1"
(median 0.07"); the Lahiri ayanamsa uses the IAU-2006 general precession,
consistent with the output frame.

## Conventions

* Outputs are referred to the **ecliptic and equinox of date** in the
  IAU-2006 frame, which is the tropical frame used in astrology.  The
  `_apparent_` functions additionally apply light-time (iterated), solar
  light deflection, annual aberration and nutation in longitude.
* The Moon uses the mean-equinox-of-date series; nutation is applied for
  apparent positions.
* The `year` argument to the `_ut` functions is only used to look up ΔT;
  pass a fractional calendar year.  ΔT is observed up to the last IERS
  data and a forecast after it (a few seconds off within a decade, tens
  of seconds by 2100); each second of ΔT error moves the Moon ~0.5" in
  UT-based positions.  Re-run `tools/gen_delta_t.py`
  occasionally to pick up new IERS data.

## Coefficient provenance and regeneration

The VSOP87D tables are downloaded from IMCCE (Observatoire de Paris) and
compiled into `data/vsop87d_data.cpp`:

```sh
python3 tools/gen_vsop87.py        # downloads data/vsop87d/* and regenerates the tables
```

The test reference file `tests/reference_data.h` is generated once from a
JPL DE440s kernel with Skyfield (reduced with ERFA/IAU-2006); the kernel is **not** part of the library
and is never shipped.  See `tools/make_reference.py` and
`tools/gen_reference.py` if you want to regenerate it.

## Compact encoding

The compiled tables use a quantized layout to keep the binary small
without measurable loss:

* **VSOP87D** — the secular (`C == 0`) terms are folded into an exact
  `double` polynomial (their amplitudes reach ~10⁴ rad/millennium and must
  not be float32); the small periodic terms are stored as
  `{float amplitude, uint32 phase, double frequency}` — 24 → 16 bytes/term.
* **ELP2000-82B** — the bounded phase (mod 2π) is stored in 30 bits, with
  the time power packed into the top 2 bits and only the rate kept as
  `double` — 48/24 → 40/16 bytes/term.
* **Pluto** — kept as `double`; its least-squares fit has large cancelling
  coefficients that float32 cannot represent.

Together these cut the tables roughly in half: a stripped, self-contained
binary is now **~670 KB** (down from ~1.15 MB) with identical accuracy.

## Optional compressed backend (~0.001 arcsec)

For near-mas accuracy the library can load a table generated from JPL
DE440s:

```c
eph_load_compressed("data/ephem.bin");
```

When loaded, positions for dates in 1850-2140 come from the table;
outside that range, or if it is not loaded, the analytical engine is used.
The table (`data/ephem.bin`, ~830 KB) stores **DE440 minus the library's
own series**: the residual is small and smooth, so it packs into
fixed-length Chebyshev segments with one absolute quantization step per
body and a fixed bit width per coefficient.  To keep it fast, the series
under the Moon and the inner planets are truncated to terms above 1e-7 rad
(~0.02"), and the outer planets and the Earth-Moon barycentre are stored
without a base series.  Each body is fitted to 0.002" (seen from the Earth
at its closest) and its speed to 0.01"/day; the Earth comes from the
barycentre and the Moon.  The format is our own, built from public-domain
JPL data, so no third-party license applies.  Regenerate with
`python3 tools/pack_ephemeris.py` (needs numpy + skyfield + pyerfa + a
DE440s kernel and a built `build/libephemeris.a`; `--target` sets the
accuracy).  A table must be regenerated whenever the series change.

## Layout

```
include/ephemeris/ephemeris.h   public C API
src/                            VSOP87, ELP2000 Moon, Pluto, asteroids,
                                DE440 corrections, compressed backend,
                                frames, time, houses
data/*_data.{h,cpp}             generated coefficient tables
data/ephem.bin                  optional DE440 residual table
tools/                          coefficient, fit and reference generators
tests/                          DE440s / Skyfield / Horizons validation,
                                topocentric, houses, compressed backend
examples/positions.c            CLI that prints a chart
```

## References

* P. Bretagnon, G. Francou, *Astron. Astrophys.* **202**, 309 (1988).
* J. Meeus, *Astronomical Algorithms*, 2nd ed. (Willmann-Bell, 1998).
* F. Espenak, J. Meeus, *Five Millennium Canon of Solar Eclipses*, NASA/TP-2006-214141.
