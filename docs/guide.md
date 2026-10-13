# Guide

[README](../README.md) · [Guide](guide.md) · [Accuracy](accuracy.md) · [Data files](data.md) · [Methods](methods.md)

How to build Ephemera and use each feature. Angles are radians, distances are AU, and time is Julian Date (TT or UT) in the C API; the JS wrapper uses degrees.

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
[Coefficient provenance](data.md#coefficient-provenance-and-regeneration)):

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

builds modules that need **no imports** (no WASI, no Emscripten
runtime: `WebAssembly.instantiate(bytes, {})` is enough):

| Module                                 | Size   | brotli | Data built in               |
|----------------------------------------|--------|--------|-----------------------------|
| `build/wasm/ephemeris.wasm`            | 118 KB | 49 KB  | none: load a file           |
| `build/wasm/ephemeris-lite.wasm`       | 100 KB | 42 KB  | none; no eclipses or stars  |
| `build/wasm/ephemeris-1as.wasm`        | 263 KB | 162 KB | `eph-1as.bin`               |
| `build/wasm/ephemeris-1as-chiron.wasm` | 322 KB | 217 KB | `eph-1as-chiron.bin`        |

They ship as the npm package **`@astrologer/ephemera`** (`bindings/js/`:
the wrapper `ephemeris.mjs`, a `load()` helper, TypeScript types and a
demo page; angles in degrees):

```sh
make npm           # stages the package in build/npm (wasm modules, data files)
make npm-test      # smoke-tests it with node and shows the tarball size (1.8 MB)
```

```js
import { load, Body } from '@astrologer/ephemera';
const eph = await load();                            // 1" data built in
const jd = eph.julianDate(2026, 1, 1.5);            // UT
const moon = eph.position(Body.Moon, jd);           // { longitude, latitude, ..., longitudeSpeed }
const { asc, mc } = eph.ascMc(jd, -0.1278, 51.5074); // observer lon, lat (degrees)
```

The lite core leaves out eclipses and fixed stars: those functions stay
but return -1 (null in the wrapper).  Native builds do the same with
`make NO_ECLIPSES=1 NO_STARS=1` or CMake's `-DEPHEMERIS_ECLIPSES=OFF
-DEPHEMERIS_STARS=OFF`.

With the core module, pass data files: `load({ module: 'core', data:
['01as', 'stars'] })`, or raw bytes such as a DE440 table; any module
can take more data the same way (or with `eph.loadData(bytes)`).

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

## Sidereal zodiac and houses

```c
double ayan = eph_ayanamsa(EPH_AYANAMSA_LAHIRI, jd_tt); /* sidereal = tropical - ayan */
double asc, mc;
eph_asc_mc(jd_ut, &obs, &asc, &mc);
double cusps[12];
eph_houses(jd_ut, &obs, EPH_HOUSES_PORPHYRY, cusps);
```

House systems: Whole Sign, Equal (from the Ascendant, the MC or 0°
Aries), Vehlow, Porphyry, Placidus, Koch, Regiomontanus, Campanus,
Alcabitius, Topocentric (Polich–Page), Meridian (axial rotation) and
Morinus, from their standard definitions; no data file is needed.
`eph_house_angles` adds the Vertex, the equatorial Ascendant (East
point), Koch's and Munkasey's co-Ascendants and the polar Ascendant.  They
match an independent geometric implementation (ERFA sidereal time and
IAU 2000A nutation; `tools/make_houses_reference.py`) to 0.011" for
1900–2030 at latitudes up to ±65° and in the polar cases, and Asc/MC use
the IAU 2006 sidereal time.  Placidus and Koch are undefined inside the
polar circles (|latitude| > 90° − obliquity): `eph_houses` then returns
Porphyry cusps and 1.  There the great-circle systems (Regiomontanus,
Campanus, Topocentric) may run backwards around the ecliptic, the MC
taking the opposite meridian point so the houses stay in order (as Swiss
Ephemeris does; for Topocentric, which has no accepted polar form, other
software may differ).  `eph_house_position(cusps, longitude)` gives the
house of a longitude as 1.0–12.999 (by longitude only).

Ayanamsas: 35 systems (`eph_ayanamsa`, names from `eph_ayanamsa_name`) --
Lahiri (and its 1940, VP285 and ICRC variants), Fagan/Bradley, Raman,
Krishnamurti (and Senthilathiban), Sri Yukteshwar, De Luce, Usha/Shashi,
Djwhal Khul, J.N. Bhasin, the Babylonian ones (Kugler 1-3, Huber, Eta
Piscium, Aldebaran, Britton), Hipparchos, Sassanian, Suryasiddhanta and
Aryabhata (true and mean Sun, SS Revati/Citra, Aryabhata 522), Vettius
Valens, Skydram, Fiorenza's galactic equator and J2000/J1900/B1950, and
seven star-based ones (see Fixed stars) -- plus user-defined ones
(`eph_ayanamsa_custom`).  Each is its defining value
carried by the IAU-2006 general precession (exactly ERFA's `p06e`).  The
J2000 values are the common realisations used in astrological software:
J2000, J1900 and B1950 match their definitions within 0.001".  Lahiri
matches the official tables of the *Indian Astronomical Ephemeris*
(Positional Astronomy Centre, IMD) within their 0.1" rounding: the 2026
edition to 0.05" (the 2019 edition, which still used the IAU 1976
precession rate, to 0.12").  Fagan/Bradley is 0.41" from its literal 1950
definition carried by IAU 2006 precession.  (Earlier versions of this
library used 23°51'11.6" for Lahiri at J2000; the official value is
23°51'25.53".)

## Lunar nodes and Lilith

`EPH_MEAN_NODE`, `EPH_TRUE_NODE`, `EPH_MEAN_LILITH` and `EPH_OSCU_LILITH`
work with the position functions like the bodies (geocentric ecliptic of
date, distance 0, no parallax; the apparent functions add nutation):

* **Mean node, mean Lilith** (mean lunar apogee): the ELP2000-82B
  polynomials as given by Meeus (*Astronomical Algorithms*, ch. 47); the
  mean apogee lies on the inclined mean orbit, so it has a latitude.  No
  data file is needed.  They agree with the full ELP2000-82B mean
  elements (IAU 2006 precession) to 0.12" / 0.54".
* **True node, osculating Lilith**: from the osculating orbit of the
  Moon's position and velocity (inertial, G(M_Earth + M_Moon) of DE440),
  with speeds from its acceleration.  Against the same elements of the
  DE440 Moon (1850–2140, max):

  | Data file    | true node | osculating Lilith |
  |--------------|-----------|-------------------|
  | DE440 tables | 0.15″     | 1.9″              |
  | full series  | 0.53″     | 2.0″              |
  | 0.1″ file    | 1.6″      | 4.8″              |
  | 1″ file      | 3.2″      | 17″               |

  The node magnifies the Moon's latitude error ~11x and the apogee its
  velocity error ~20x.  (Swiss Ephemeris' Moshier mode is 37" / 233" off.)

## Event searches and time

* `eph_longitude_crossing(_ut)`: when a body reaches a longitude --
  solar and lunar returns, ingresses, transits to natal points (retrograde
  loops are followed).  `eph_heliocentric_crossing` does the same for the
  heliocentric longitude.
* `eph_longitude_crossing_sidereal` (and `_custom`, and the heliocentric
  forms) give sidereal returns and ingresses: the longitude minus the
  ayanamsa *at the crossing* reaches the target, in one search (the target
  moves ~50"/year with the precession).
* `eph_node_crossing`: when a body crosses the ecliptic (the Moon's node
  passages); `eph_station`: when a planet turns retrograde or direct.
* Results are refined to ~1e-8 day; with the full series, solar and lunar
  crossings agree with JPL DE440 (via Skyfield) within 0.13 s / 0.02 s.
* `eph_phenomena_tt/_ut`: phase angle, illuminated fraction, elongation,
  apparent diameter and visual magnitude (Mallama & Hilton 2018 for the
  planets, Allen for the Moon to phase angle ~150°, the H-G system for the
  minor bodies); against DE440 with Skyfield's Mallama & Hilton
  magnitudes: within 0.05" in angle and 0.002 mag for the planets.
* `eph_utc_to_jd` (leap seconds; exact TT from 1972), `eph_sidereal_time`
  (IAU 2006), `eph_equation_of_time`, `eph_midpoint` and
  `eph_split_longitude`.

## Eclipses

Solar and lunar eclipses from the Sun and Moon positions, so any data file
works (no extra data):

```c
eph_solar_eclipse e;
eph_solar_eclipse_next(jd_tt, +1, &e);   /* -1: previous */
/* e.type: EPH_ECLIPSE_TOTAL / ANNULAR / HYBRID / PARTIAL (+ NONCENTRAL),
 * e.maximum, e.gamma, e.magnitude, e.saros, where greatest eclipse is
 * seen, and when the penumbra, umbra and axis touch the Earth. */
```

* `eph_solar_eclipse_next`: global circumstances, from Besselian elements
  computed from the apparent positions (Explanatory Supplement, 1992),
  with the Earth's flattening and NASA's lunar radii (k = 0.272488
  penumbra, 0.272281 umbra).
* `eph_solar_eclipse_local`: the next eclipse seen from a place: contacts,
  magnitude, obscuration and the Sun's altitude; `eph_solar_eclipse_how`
  gives the same at any instant, and `eph_solar_eclipse_where` the point
  under the shadow axis.
* `eph_lunar_eclipse_next`: total, partial and penumbral eclipses, with
  the umbral and penumbral magnitudes and the contacts; the shadow is
  enlarged for the atmosphere by Danjon's rule.

Against NASA's *Five Millennium Canon* for all 909 eclipses of 1901–2100,
types and saros series agree everywhere. Gamma and magnitudes agree within
0.0002, the greatest eclipse within 0.6 s (lunar 3 s), and lunar phase
durations within 0.5 min.  Against a brute-force reference from DE440
(Skyfield; `tools/make_eclipse_reference.py`), contacts agree within 0.3 s
(local, and global for non-grazing eclipses; up to 2 s where the shadow
only grazes the Earth), local magnitudes within 0.0003 and obscuration
within 0.002. With the 1″ file the only change is 1986 Oct 3, a hybrid of
magnitude 1.0000 that becomes annular.  A search takes ~8 ms with the full
series and ~1 ms with the 1″ file.

## Fixed stars

`data/eph-stars.bin` (26 KB) holds the 410 stars with IAU-approved names
that have a Hipparcos number, plus the Galactic Center (Sgr A*):
positions, proper motions and parallaxes from the Hipparcos new reduction
(van Leeuwen 2007), radial velocities from SIMBAD, names from the IAU
Working Group on Star Names; built by `tools/make_stars.py`.

```c
int i = eph_star_find("Regulus");        /* or "alf Leo", "HIP 49669" */
eph_position p;
eph_star_position_apparent_ut(i, jd_ut, year, &p);
```

Positions follow the stars' space motion (with parallax and radial
velocity, as ERFA's `pmsafe`) and then the same reduction as the planets
(aberration with the Earth's barycentric velocity); they agree with ERFA +
Skyfield/DE440 within 0.015" over 1600–2500.  The star-based ayanamsas
(True Citra, True Revati, True Pushya, True Mula, Galactic Center at 0°
Sagittarius, Gil Brand's and Cochrane's) use this file and agree with
the same definitions computed with ERFA + DE440 within 0.03".

## Minor planets

Two optional data files add numbered asteroids, centaurs and
trans-Neptunian objects over 1900–2100. Load them like any other data
file; they are not in the npm package. Download them from the
repository (`data/`):

| File | Bodies | Size |
|---|---|---|
| `eph-minor.bin` | 180: the classic numbered asteroids 5–100 (Astraea, Hebe, Hygiea, Psyche, Urania, Pandora, Sappho, …), near-Earth objects (Eros, Lilith 1181, Amor, Icarus, Apollo, Toro, Phaethon, Apophis, Bennu, Ryugu, …), centaurs (Chiron, Pholus, Nessus, Chariklo, Asbolus, …), TNOs (Eris, Sedna, Makemake, Haumea, Quaoar, Orcus, Varuna, Ixion, Gonggong, …) and other named favourites | 5.7 MB |
| `eph-minor-1000.bin` | all numbered asteroids 1–1000 | 17.2 MB |

```c
eph_load_file("data/eph-minor.bin");
int i = eph_minor_find("Eros");            /* or "433", "433 Eros" */
eph_position p;
eph_position_apparent_ut(EPH_MINOR_BASE + i, jd_ut, year, &p);
```

Body `EPH_MINOR_BASE + index` works with the position, heliocentric,
search, phenomena (H, G and diameter from the JPL Small-Body Database)
and orbit functions. In JS, `eph.minorBody("Eros")` returns that id;
`minorList()` lists what is loaded.

How the files are made (`tools/make_minor_planets.py`):

1. Each body starts from its osculating elements in the JPL Small-Body
   Database.
2. It is integrated over the span with the Sun, planets, Moon and Pluto
   from DE440, the 16 massive asteroids of JPL's own ephemerides, and the
   Sun's relativistic term.
3. Every body is checked against JPL Horizons. Bodies that miss, which
   are near-Earth objects with close approaches, are integrated again
   from Horizons states every 10, 2 or ½ years.
4. The positions are stored as Chebyshev pieces, accurate to 0.5″ as seen
   from the closest the Earth can get.

Against Horizons the apparent positions agree within 0.25″ (55 test
positions). The exception is Hermes (69230): its Horizons ephemeris is
consistent with itself only to ~2×10⁻⁵ AU around its long-lost years.
Apophis is chaotic after its 2029 Earth flyby, so Horizons must be
queried with dates on one side of its orbit epoch at a time; the tools
do this.

## Uranian points

`EPH_CUPIDO`, `EPH_HADES`, `EPH_ZEUS`, `EPH_KRONOS`, `EPH_APOLLON`,
`EPH_ADMETOS`, `EPH_VULKANUS` and `EPH_POSEIDON` are the Hamburg School's
hypothetical points.  Cupido to Kronos were defined by Alfred Witte
(1923–24) and the other four by Friedrich Sieggrün (1927, 1934).  They
work with all the position, heliocentric and search functions.

Their orbits are the elements Paul Schlyter fitted to the published
Witte/Sieggrün ephemeris ("Hypothetical astrological planets",
alt.astrology, 1994):

* a mean longitude that moves uniformly in the tropical zodiac;
* for the four inner points, the first term of the equation of centre;
* for Cupido and Hades, their small inclinations.

Schlyter found that these fit the ephemeris to a few arcminutes, so that
is the accuracy to expect against the printed *Transneptun-Ephemeride*
(Witte-Verlag).  Our copy of the table was checked for its mean
longitudes and rates (1900.0 to 1950.0) and its distances (Kepler's third
law from the periods).  The small eccentricity, perihelion and node terms
of Cupido to Kronos could not be cross-checked against a second copy; they
move the positions by at most ~0.5°.  Other software may use revised element sets and can
differ by a similar amount.  Geocentric positions add the Earth's position
and the same aberration and nutation as the stars.

## Planetary nodes and apsides

`eph_orbit_tt` gives a body's heliocentric orbit in the mean ecliptic and
equinox of date: semi-major axis, eccentricity, inclination, the
longitudes of the ascending node and of perihelion, and four points on
the orbit (ascending and descending node, perihelion, aphelion), each
heliocentric and geocentric, with speeds:

* `EPH_ORBIT_MEAN`: the secular mean elements of Simon et al. (1994) as
  tabulated by Meeus (*Astronomical Algorithms*, Table 31.A), for
  Mercury to Neptune.  They reproduce Meeus' worked example 31.a exactly.
* `EPH_ORBIT_OSCULATING`: from the position and velocity, for the planets
  except the Earth, plus Pluto and the minor bodies.  With the DE440 tables
  loaded, the elements and points agree with DE440 (Skyfield's osculating
  elements) within 0.8".  With the series alone they agree within 0.5" for
  Mercury to Saturn, but only 1" for Uranus, 3" for Neptune and 11" for
  Pluto, because osculating elements magnify the series' small errors.

The geocentric points are geometric, with nutation added so they compare
with the apparent positions.
