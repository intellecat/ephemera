/*
 * ephemeris.h -- lightweight analytical ephemeris engine.
 *
 * Positions of the Sun, Moon, Pluto, the eight planets and
 * Ceres/Pallas/Juno/Vesta using the VSOP87D planetary theory (Bretagnon &
 * Francou 1988), the ELP2000-82B lunar theory and compact series fitted to
 * JPL data, with small DE440 corrections.  The series coefficients are
 * not compiled in: load a data file first (see "data files" below).
 * Time, nutation, frames, sidereal time and houses need no data.
 *
 * Conventions:
 *   - Angles are radians on output.
 *   - Ecliptic coordinates are heliocentric/geocentric, referred to the
 *     mean ecliptic and equinox of date (i.e. tropical), IAU 2006 frame.
 *   - Distances are astronomical units (the Moon too, for consistency;
 *     1 AU = 149597870.7 km).
 *   - Time arguments are Julian dates.  Functions suffixed _tt take
 *     Terrestrial Time; functions suffixed _ut take Universal Time and
 *     apply the Delta-T correction internally.
 *   - All time in this library is astronomical (no calendars), so the
 *     Julian date you pass for a calendar instant should be the JD of the
 *     UT or TT instant you mean.
 */
#ifndef EPHEMERIS_H
#define EPHEMERIS_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
    EPH_MERCURY = 0,
    EPH_VENUS,
    EPH_EARTH,
    EPH_MARS,
    EPH_JUPITER,
    EPH_SATURN,
    EPH_URANUS,
    EPH_NEPTUNE,
    EPH_PLUTO,
    EPH_SUN,
    EPH_MOON,
    EPH_CERES,
    EPH_PALLAS,
    EPH_JUNO,
    EPH_VESTA,
    EPH_BODY_COUNT
} eph_body;

typedef struct {
    double longitude;      /* geocentric ecliptic longitude, radians [0,2pi) */
    double latitude;       /* geocentric ecliptic latitude, radians          */
    double distance;       /* geocentric distance, AU                         */
    double longitude_speed;/* d(longitude)/dt, radians/day                    */
    double latitude_speed; /* d(latitude)/dt, radians/day                     */
    double distance_speed; /* d(distance)/dt, AU/day                          */
} eph_position;

/* ---- time ---------------------------------------------------------- */

/* Julian date from a Gregorian (or Julian before 1582-10-15) calendar date.
 * `day` carries the fractional part of the day (e.g. 1.5 = noon). */
double eph_julian_date(int year, int month, double day);

/* Inverse of eph_julian_date.  Gregorian calendar for jd >= 2299160.5. */
void eph_calendar_date(double jd, int* year, int* month, double* day);

/* Approximate Delta-T = TT - UT in seconds (Espenak & Meeus polynomials). */
double eph_delta_t(double year);

/* TDB - TT in seconds (Fairhead & Bretagnon periodic series, ~1 us).
 * Dynamical theories such as VSOP87/ELP2000 are functions of TDB. */
double eph_tdb_minus_tt(double jd_tt);

/* Convert between Universal Time and Terrestrial Time.  The argument is a
 * calendar year with fraction (used only to evaluate Delta-T). */
double eph_ut_to_tt(double jd_ut, double year);
double eph_tt_to_ut(double jd_tt, double year);

/* ---- reference frames ---------------------------------------------- */

/* Mean obliquity of the ecliptic (radians), IAU 2006. */
double eph_mean_obliquity(double jd_tt);

/* Nutation in longitude and obliquity (radians), IAU 2000B, 77 terms. */
void eph_nutation(double jd_tt, double* dpsi, double* deps);

/* Convert apparent (true-equinox-of-date) ecliptic coordinates, such as
 * those returned by eph_position_apparent_tt, to equatorial RA/Dec using
 * the true obliquity. */
void eph_ecliptic_to_equatorial(double lon, double lat, double jd_tt,
                                double* ra, double* dec);

/* ---- positions ----------------------------------------------------- */

/* Geometric geocentric position (no light-time, aberration or nutation),
 * mean ecliptic and equinox of date. */
int eph_position_tt(eph_body body, double jd_tt, eph_position* out);

/* Apparent geocentric position: light-time + annual aberration + nutation. */
int eph_position_apparent_tt(eph_body body, double jd_tt, eph_position* out);

/* UT wrappers.  `year` is the (fractional) calendar year for Delta-T. */
int eph_position_ut(eph_body body, double jd_ut, double year, eph_position* out);
int eph_position_apparent_ut(eph_body body, double jd_ut, double year,
                             eph_position* out);

/* Heliocentric ecliptic-of-date position.  Fails for Sun and Moon. */
int eph_heliocentric_tt(eph_body body, double jd_tt,
                        double* lon, double* lat, double* distance);

/* ---- topocentric (observer) positions ------------------------------- */

typedef struct {
    double longitude; /* geographic, radians, east positive */
    double latitude;  /* geodetic, radians */
    double elevation; /* metres above the reference ellipsoid */
} eph_observer;

/* Apparent topocentric position (geocentric apparent place with diurnal
 * parallax applied).  The UT is needed for the sidereal time; `year` is
 * the fractional calendar year for Delta-T.  Atmospheric refraction is not
 * applied (it affects altitude, not the ecliptic longitude used in charts). */
int eph_position_topocentric_ut(eph_body body, double jd_ut, double year,
                                const eph_observer* observer,
                                eph_position* out);

/* As above with an explicit TT for the ephemeris and UT for sidereal time. */
int eph_position_topocentric_tt(eph_body body, double jd_tt, double jd_ut,
                                const eph_observer* observer,
                                eph_position* out);

/* ---- sidereal zodiac and houses ------------------------------------ */

/* Lahiri ayanamsa (radians): sidereal longitude = tropical - ayanamsa. */
double eph_ayanamsa_lahiri(double jd_tt);

/* Ascendant and Midheaven ecliptic longitudes (radians, tropical). */
int eph_asc_mc(double jd_ut, const eph_observer* observer, double* asc,
               double* mc);

typedef enum {
    EPH_HOUSES_WHOLE_SIGN = 0,
    EPH_HOUSES_EQUAL,
    EPH_HOUSES_PORPHYRY
} eph_house_system;

/* House cusps 1..12 as ecliptic longitudes (radians, tropical). */
int eph_houses(double jd_ut, const eph_observer* observer,
               eph_house_system system, double cusps[12]);

const char* eph_body_name(eph_body body);

/* ---- data files ------------------------------------------------------ *
 * Positions need series data from a file made by tools/pack_series or
 * tools/pack_ephemeris.py.  Pick one series level:
 *   data/eph-1as.bin    ~1"   (truncated series)
 *   data/eph-01as.bin   ~0.1" (truncated series)
 *   data/eph-full.bin   the full series (~0.05")
 * and optionally add data/eph-asteroids.bin (Ceres, Pallas, Juno, Vesta)
 * and/or data/ephem.bin (DE440 residual table, ~0.001" over 1850-2140,
 * self-contained: it works without a series file inside its range).
 * Several files can be loaded; a later file overrides what an earlier one
 * provides.  Position functions return -1 for bodies with no data loaded.
 * Loading is not thread-safe: load before computing from other threads. */

/* Load a data file from memory (the buffer is copied) or from a path.
 * Returns 0 on success, -1 if the data is invalid. */
int eph_load_data(const void* data, size_t size);
int eph_load_file(const char* path);

/* Load a data file used in place, without a copy: the buffer must be
 * 8-byte aligned and stay valid and unchanged until it is unloaded (e.g. a
 * static array or a memory-mapped file). */
int eph_load_static(const void* data, size_t size);

/* Builds with EPHEMERIS_EMBED_DATA carry a data file (e.g. eph-1as.bin)
 * that is loaded when the program starts; this loads it again after
 * eph_unload_data().  Returns -1 if the library has none. */
int eph_load_embedded(void);

/* Unload all data files. */
void eph_unload_data(void);

/* Earlier names: eph_load_compressed loads any data file,
 * eph_unload_compressed unloads the files holding a residual table, and
 * eph_compressed_available tells whether a residual table is loaded. */
int eph_load_compressed(const char* path);
void eph_unload_compressed(void);
int eph_compressed_available(void);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* EPHEMERIS_H */
