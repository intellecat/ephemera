# Methods

[README](../README.md) · [Guide](guide.md) · [Accuracy](accuracy.md) · [Data files](data.md) · [Methods](methods.md)

## Theories and models

Internally it evaluates:

* **VSOP87D** (Bretagnon & Francou, 1988) for Mercury…Neptune — the
  IMCCE coefficient tables (23,189 terms in full) for heliocentric spherical
  coordinates in the *mean ecliptic and equinox of date*.
* **ELP2000-82B** (Chapront-Touze & Chapront) full lunar theory for the
  Moon; ~0.1" near J2000, growing to a few arcseconds at the range edges
  because its constants are fitted to DE200/LE200.
* Compact Chebyshev series for Ceres/Pallas/Juno/Vesta/Chiron fitted to JPL
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

## References

* P. Bretagnon, G. Francou, *Astron. Astrophys.* **202**, 309 (1988).
* J. Meeus, *Astronomical Algorithms*, 2nd ed. (Willmann-Bell, 1998).
* F. Espenak, J. Meeus, *Five Millennium Canon of Solar Eclipses*, NASA/TP-2006-214141;
  *Five Millennium Canon of Lunar Eclipses*, NASA/TP-2009-214172 (the
  catalogues are the eclipse test reference, `tools/make_nasa_eclipse_reference.py`).
* P. K. Seidelmann (ed.), *Explanatory Supplement to the Astronomical
  Almanac* (University Science Books, 1992), ch. 8.

## Data sources

The data the library and its data files are built from:

* JPL DE440 / DE440s and Horizons, NASA's eclipse catalogues: US
  government works, public domain.
* VSOP87 and ELP2000-82B (IMCCE, Bureau des longitudes): published
  theories, free to use with citation.
* Hipparcos new reduction (ESA; van Leeuwen 2007), SIMBAD radial
  velocities (CDS, Strasbourg), IAU Working Group on Star Names: free to
  use with acknowledgement.
* The *Indian Astronomical Ephemeris* (Positional Astronomy Centre, IMD):
  used only as test reference for the Lahiri ayanamsa.

The tests compare against JPL DE440 (via Skyfield and ERFA), NASA's
eclipse catalogues, the Indian Astronomical Ephemeris and independent
implementations in `tools/`; no third-party code is included.
