/*
 * ephemeris.h -- lightweight analytical ephemeris engine.
 *
 * Pure-code positions of the Sun, Moon, Pluto and the eight planets using
 * the VSOP87D planetary theory (Bretagnon & Francou 1988) and the
 * Meeus/ELP2000 truncated lunar theory.  No external data file is needed:
 * the series coefficients are compiled into the library.
 *
 * Conventions:
 *   - Angles are radians on output.
 *   - Ecliptic coordinates are heliocentric/geocentric, referred to the
 *     mean dynamical ecliptic and equinox of date (i.e. tropical).
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

/* Mean obliquity of the ecliptic (radians), IAU 1980. */
double eph_mean_obliquity(double jd_tt);

/* Nutation in longitude and obliquity (radians), IAU 1980, 63 terms. */
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

const char* eph_body_name(eph_body body);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* EPHEMERIS_H */
