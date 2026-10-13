# Accuracy

[README](../README.md) · [Guide](guide.md) · [Accuracy](accuracy.md) · [Data files](data.md) · [Methods](methods.md)

Every figure here comes from a test or tool in this repository; see also the per-feature figures in the [guide](guide.md).

## At a glance

| What | Checked against | Agreement |
|---|---|---|
| Planets, Sun, Moon (1550–2650) | JPL DE440 | 0.03″ (DE440 tables), 0.15″ (full series), 0.6″ (145 KB file) |
| Houses (14 systems) and angles | independent geometric implementation with ERFA | 0.011″ |
| Ayanamsas | IAU 2006 precession (ERFA); *Indian Astronomical Ephemeris* (Lahiri) | exact; 0.05″ |
| Lunar nodes, Lilith | ELP2000-82B mean elements; DE440 osculating orbit | 0.12″ / 0.54″ (mean); 0.15″ / 1.9″ (true, DE440 tables) |
| Fixed stars (411) | ERFA + DE440 | 0.015″ |
| Eclipses (909, 1901–2100) | NASA *Five Millennium Canon* | types and saros exact; gamma, magnitude 0.0002; maximum 0.6 s |
| Eclipse contacts | brute-force DE440 reference | 0.3 s (2 s for grazing eclipses) |
| Returns, ingresses | JPL DE440 | Sun 0.13 s, Moon 0.02 s |
| Phenomena (phase, magnitude) | DE440 + Mallama & Hilton 2018 | 0.05″, 0.002 mag |
| Planetary nodes and apsides | VSOP87 secular terms; DE440 osculating elements | exact; 0.8″ |
| Minor planets (1,083) | JPL Horizons | 0.25″ |
| Hand-entered constants | IERS, JPL SBDB, DE440 constants, VSOP87 (`tools/audit_constants.py`) | all match |


## Planets, Sun and Moon

Maximum error in apparent geocentric ecliptic longitude/latitude against
JPL DE440 on 3000 random dates in **1550–2650**, for each data file
(`tools/check_accuracy.py`).  The reference is DE440 via Skyfield in the
IAU-2006/2000A ecliptic of date -- the frame Swiss Ephemeris and Horizons
report:

| Body    | full series | 0.1″ file | 1″ file | DE440 tables |
|---------|-------------|-----------|---------|--------------|
| Mercury | 0.04″       | 0.09″     | 0.50″   | 0.03″        |
| Venus   | 0.05″       | 0.11″     | 0.58″   | 0.03″        |
| Mars    | 0.14″       | 0.14″     | 0.55″   | 0.03″        |
| Jupiter | 0.09″       | 0.10″     | 0.32″   | 0.03″        |
| Saturn  | 0.04″       | 0.05″     | 0.49″   | 0.03″        |
| Uranus  | 0.04″       | 0.07″     | 0.46″   | 0.03″        |
| Neptune | 0.03″       | 0.04″     | 0.47″   | 0.03″        |
| Pluto   | 0.04″       | 0.04″     | 0.04″   | 0.03″        |
| Sun     | 0.04″       | 0.04″     | 0.11″   | 0.03″        |
| Moon    | 0.15″       | 0.22″     | 0.47″   | 0.03″        |
| median  | 0.004–0.03″ | 0.005–0.04″ | 0.01–0.12″ | 0.002″   |

(Maximum of two runs with different dates.)  Within 1850–2140 the errors are smaller still: the DE440 tables reach
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
