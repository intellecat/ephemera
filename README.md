# ephemeris

A small, dependency-free C/C++ library that computes positions of the Sun,
Moon, Pluto, the eight planets and the four largest asteroids (Ceres,
Pallas, Juno, Vesta) from analytical series.  No JPL kernel and no SPICE:
a ~44 KB core (as WebAssembly; ~20 KB brotli) loads one small data file
with the series coefficients at the accuracy you need -- from 145 KB for
~0.5" to DE440 tables at ~320 KB per century for ~0.002", over 1550–2650
(see [Data files](#data-files)).

Internally it evaluates:

* **VSOP87D** (Bretagnon & Francou, 1988) for Mercury…Neptune — the
  IMCCE coefficient tables (23,189 terms in full) for heliocentric spherical
  coordinates in the *mean ecliptic and equinox of date*.
* **ELP2000-82B** (Chapront-Touze & Chapront) full lunar theory for the
  Moon; ~0.1" near J2000, growing to a few arcseconds at the range edges
  because its constants are fitted to DE200/LE200.
* Compact Chebyshev series for Ceres/Pallas/Juno/Vesta fitted to JPL
  Horizons (public data), heliocentric fit error <0.03" over 1600-2500
  (Horizons' span for them).  Apparent positions match Horizons to ~0.02".
* Chebyshev intervals for Pluto fitted to JPL DE440 (0.03" over 1550-2650),
  and DE440 corrections for VSOP87D and ELP2000-82B over 1550-2650.
* **IAU 2000B** nutation (77 terms), a rotation from the series' native
  IAU-1976 frame to the **IAU-2006** ecliptic of date (the frame Swiss
  Ephemeris and JPL Horizons use), and a yearly ΔT table for UT↔TT
  (IERS observations; Stephenson, Morrison & Hohenkerk 2016 before 1973
  and as the forecast to 2150; Espenak–Meeus polynomials outside
  1620–2150).

Angles are radians, distances are AU, and time is Julian Date (TT or UT).

## Accuracy

Maximum error in apparent geocentric ecliptic longitude/latitude against
JPL DE440 on 3000 random dates in **1550–2650**, for each data file
(`tools/check_accuracy.py`).  The reference is DE440 via Skyfield in the
IAU-2006/2000A ecliptic of date -- the frame Swiss Ephemeris and Horizons
report:

| Body    | full series | 0.1″ file | 1″ file | DE440 tables |
|---------|-------------|-----------|---------|--------------|
| Mercury | 0.04″       | 0.08″     | 0.48″   | 0.04″        |
| Venus   | 0.05″       | 0.10″     | 0.58″   | 0.03″        |
| Mars    | 0.11″       | 0.13″     | 0.55″   | 0.03″        |
| Jupiter | 0.09″       | 0.10″     | 0.32″   | 0.03″        |
| Saturn  | 0.04″       | 0.05″     | 0.49″   | 0.03″        |
| Uranus  | 0.04″       | 0.07″     | 0.46″   | 0.03″        |
| Neptune | 0.03″       | 0.04″     | 0.47″   | 0.03″        |
| Pluto   | 0.04″       | 0.04″     | 0.04″   | 0.03″        |
| Sun     | 0.03″       | 0.04″     | 0.11″   | 0.03″        |
| Moon    | 0.12″       | 0.17″     | 0.47″   | 0.03″        |
| median  | 0.004–0.03″ | 0.005–0.04″ | 0.01–0.11″ | 0.002″   |

Within 1850–2140 the errors are smaller still: the DE440 tables reach
0.004" max / ~0.001" median there, on par with Swiss Ephemeris' .se1
files (0.0005–0.0007" median, 0.002" for the Moon); towards 1550 and 2650
their maximum grows to ~0.03" with the nutation model (below).  The 1″
file is comparable to Swiss Ephemeris' Moshier mode in size and better
for the Moon and the outer planets (Moshier, 1850–2140: Moon 3.2",
Neptune 2.1", Pluto 2.5" max).  Outside 1550–2650 positions come from the
series alone and degrade gradually (Pluto has no data there).

### Where the remaining error comes from

* **Planets and Moon.**  VSOP87D and ELP2000-82B are fitted to DE200/
  LE200.  Against DE440 they drift slowly (Uranus/Neptune by ~10" over
  1550–2650) and, because their orbital elements' rates differ slightly,
  periodically with growing amplitude (the Moon by 49", Mars by 2.6").
  `tools/fit_correction.py` fits both: a Chebyshev polynomial in time plus
  slowly varying amplitudes of each body's fundamental arguments
  (harmonics of a planet's mean longitude; 20 combinations of the Moon's
  Delaunay arguments).  The corrections leave 0.01–0.07" heliocentric and
  0.16" for the Moon, and are clamped to 1550–2650.  The Moon's apparent
  place applies light-time (~0.7") and nutation but no annual aberration,
  matching the reference to ~0.001".
* **Frame.**  The series are evaluated in their native IAU-1976 frame and
  rotated to the IAU-2006 ecliptic of date (`frame_1976_to_2006`, checked
  against ERFA to <0.001 mas).  Without it longitudes drift by 0.3"/century
  and the ecliptic is tilted by 0.04".  Nutation is IAU 2000B, within
  ~1 mas of IAU 2000A near J2000, ~3 mas at 1850/2140 and ~30 mas at
  1550/2650 (IAU 1980 was off by up to 18 mas over 1850–2140).
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
houses) takes, on one core of a 2024 desktop CPU, ~1.4 ms with the full
series, ~0.6 ms with the 0.1″ file, ~0.22 ms with the 1″ file and
~0.13 ms with the DE440 tables.

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
public API is plain C, so it can be used from C or C++.  The CLI loads
`data/eph-full.bin`, `data/eph-asteroids.bin` and the DE440 tables in
`data/de440/` from the current directory, or the files listed in
`$EPH_DATA` (`:`-separated).

The data files are generated from the coefficient tables in `data/` (see
[Coefficient provenance](#coefficient-provenance-and-regeneration)):

```sh
make tools && make data     # or cmake -DEPHEMERIS_BUILD_TOOLS=ON
```

**Built-in data.**  To get a library that works without any data file,
compile one in; it is loaded (in place, without a copy) when the program
starts:

```sh
cmake -S . -B build-1as -DEPHEMERIS_EMBED_DATA=data/eph-1as.bin
make EMBED=data/eph-1as.bin BUILD=build-1as      # Makefile equivalent
```

Other files can still be loaded on top (e.g. the asteroids or the DE440
table), and `eph_load_embedded()` reloads the built-in data after
`eph_unload_data()`.

## WebAssembly

```sh
make wasm          # needs zig, or WASM_CXX="clang++ --target=wasm32-wasi --sysroot=..."
make wasm-test     # compares the modules with the native library (needs node)
```

builds two modules that need **no imports** (no WASI, no Emscripten
runtime: `WebAssembly.instantiate(bytes, {})` is enough):

| Module                          | Size   | brotli | Data                         |
|---------------------------------|--------|--------|------------------------------|
| `build/wasm/ephemeris.wasm`     | 44 KB  | ~20 KB | none: load a file            |
| `build/wasm/ephemeris-1as.wasm` | 159 KB | ~92 KB | `eph-1as.bin` built in (~1″) |

`examples/wasm/ephemeris.mjs` wraps them for browsers and Node (angles in
degrees), and `examples/wasm/index.html` is a small demo page:

```js
import { createEphemeris, Body } from './ephemeris.mjs';
const eph = await createEphemeris(await (await fetch('ephemeris-1as.wasm')).arrayBuffer());
const jd = eph.julianDate(2026, 1, 1.5);            // UT
const moon = eph.position(Body.Moon, jd);           // { longitude, latitude, ..., longitudeSpeed }
const { asc, mc } = eph.ascMc(jd, -0.1278, 51.5074); // observer lon, lat (degrees)
```

With the core module, fetch a data file and pass it to
`eph.loadData(bytes)` first; the asteroid file and the DE440 table load
the same way into either module.

## Using the API

```c
#include "ephemeris/ephemeris.h"

eph_load_file("data/eph-1as.bin");               /* once, at start-up   */

double jd_ut = eph_julian_date(1992, 4, 12.0);   /* 1992-04-12 00:00 UT */
double year  = 1992.25;                          /* for Delta-T         */

eph_position moon;
eph_position_apparent_ut(EPH_MOON, jd_ut, year, &moon);

printf("Moon ecliptic longitude: %.6f deg\n", moon.longitude * 180.0 / M_PI);
printf("Moon daily motion:       %.6f deg/day\n",
       moon.longitude_speed * 180.0 / M_PI);
```

In a browser, fetch the file and pass the bytes to
`eph_load_data(ptr, size)` (the buffer is copied).  Position functions
return -1 for a body whose data is not loaded; time, nutation, sidereal
time, Asc/MC and houses need no data.

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
written to `data/vsop87d_data.cpp`; ELP2000-82B, Pluto, the DE440
corrections and the asteroid fits have their own generators in `tools/`
(the DE440 fits use the full `de440.bsp` kernel, 1550–2650, which is a
development download and never shipped).
These generated sources are not compiled into the library:
`tools/pack_series.cpp` packs them into the data files.

```sh
python3 tools/gen_vsop87.py        # downloads data/vsop87d/* and regenerates the tables
make tools && make data            # rebuilds data/eph-*.bin
```

The test reference file `tests/reference_data.h` is generated once from a
JPL DE440s kernel with Skyfield (reduced with ERFA/IAU-2006); the kernel is **not** part of the library
and is never shipped.  See `tools/make_reference.py` and
`tools/gen_reference.py` if you want to regenerate it.

## Compact encoding

The coefficient tables use a quantized layout to keep the files small
without measurable loss:

* **VSOP87D** — the secular (`C == 0`) terms are folded into an exact
  `double` polynomial (their amplitudes reach ~10⁴ rad/millennium and must
  not be float32); the small periodic terms are stored as
  `{float amplitude, uint32 phase, double frequency}` — 24 → 16 bytes/term.
* **ELP2000-82B** — the bounded phase (mod 2π) is stored in 30 bits, with
  the time power packed into the top 2 bits and only the rate kept as
  `double` — 48/24 → 40/16 bytes/term.
* **Pluto, asteroids, corrections** — float32 Chebyshev coefficients.

## Data files

| File                       | Size   | brotli | Contents                                      |
|----------------------------|--------|--------|-----------------------------------------------|
| `data/eph-1as.bin`         | 145 KB | 112 KB | series truncated to ~0.5″ (4,719 terms)        |
| `data/eph-01as.bin`        | 303 KB | 213 KB | series truncated to ~0.15″ (14,176 terms)      |
| `data/eph-full.bin`        | 649 KB | 420 KB | the full series (34,530 terms)                 |
| `data/eph-asteroids.bin`   | 365 KB | 341 KB | Ceres, Pallas, Juno, Vesta, 1600–2500          |
| `data/de440/de440-YEAR.bin`| ~325 KB each | ~305 KB | DE440 table for YEAR..YEAR+100, ~0.002″ |

The series files cover 1550–2650 (Pluto and the corrections) and work,
less accurately, outside it.  Load one series file, plus the asteroids if
you need them.  The series files include Pluto and the DE440 corrections;
truncation drops the VSOP87/ELP terms whose amplitude (times T^α at
1550/2650) is below a per-planet threshold chosen so the level holds over
1550–2650 (see `LEVEL_*` in the `Makefile`).

Several files can be loaded at once; a later file overrides what an
earlier one provides.  Each file is a small container (`src/data.h`):
typed sections used in place after loading, with a CRC-32, so a damaged
or mismatched file is rejected rather than giving wrong positions.

**DE440 tables.**  `data/de440/` holds one file per century (1550–1600,
1600–1700, …, 2600–2650; 3.7 MB in all): load the ones you need, side by
side.  Each stores **DE440 minus the library's own series**: the residual is small and smooth, so it packs into fixed-length
Chebyshev segments with one absolute quantization step per body and a
fixed bit width per coefficient.  The series under the Moon and the inner
planets are truncated to terms above 1e-7 rad (~0.02") and stored in the
same file (41 KB), so a table needs no other file, and the table records
a checksum of that base, so tables built against different series are
rejected; the outer planets and
the Earth-Moon barycentre are stored without a base series.  Each body is
fitted to 0.002" (seen from the Earth at its closest) and its speed to
0.01"/day; the Earth comes from the barycentre and the Moon.  For dates
no loaded table covers, a series file is used if one is loaded.
Regenerate with `python3 tools/pack_ephemeris.py` (needs numpy + skyfield
+ pyerfa, the DE440 kernel, and `make && make tools`; `--target` sets the
accuracy, `--chunk` the years per file).
The format is our own, built from public-domain JPL data, so no
third-party license applies.

## Layout

```
include/ephemeris/ephemeris.h   public C API
src/                            VSOP87, ELP2000 Moon, Pluto, asteroids,
                                DE440 corrections, data loader, DE440
                                table, frames, time, houses
data/*_data.{h,cpp}             generated coefficient tables (sources of
                                the data files; only delta_t is compiled in)
data/*.bin                      data files
tools/                          coefficient, fit, data-file and reference
                                generators
tests/                          DE440s / Skyfield / Horizons validation,
                                topocentric, houses, data files, DE440 table
examples/positions.c            CLI that prints a chart
```

## References

* P. Bretagnon, G. Francou, *Astron. Astrophys.* **202**, 309 (1988).
* J. Meeus, *Astronomical Algorithms*, 2nd ed. (Willmann-Bell, 1998).
* F. Espenak, J. Meeus, *Five Millennium Canon of Solar Eclipses*, NASA/TP-2006-214141.
