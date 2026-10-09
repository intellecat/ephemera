/*
 * ephemeris.h -- lightweight analytical ephemeris engine.
 *
 * Positions of the Sun, Moon, Pluto, the eight planets and
 * Ceres/Pallas/Juno/Vesta/Chiron using the VSOP87D planetary theory (Bretagnon &
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
    EPH_CHIRON,
    /* Lunar points (directions, distance 0; apparent adds nutation): */
    EPH_MEAN_NODE,   /* mean ascending node (no data file needed)          */
    EPH_TRUE_NODE,   /* osculating ascending node (needs the Moon's data)  */
    EPH_MEAN_LILITH, /* mean lunar apogee, Black Moon (no data needed)     */
    EPH_OSCU_LILITH, /* osculating lunar apogee (needs the Moon's data)    */
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

/* Delta-T = TT - UT in seconds: a yearly table for 1620-2150 (IERS
 * observations, Stephenson, Morrison & Hohenkerk 2016 and their forecast),
 * Espenak & Meeus polynomials outside it. */
double eph_delta_t(double year);

/* TDB - TT in seconds (Fairhead & Bretagnon periodic series, ~1 us).
 * Dynamical theories such as VSOP87/ELP2000 are functions of TDB. */
double eph_tdb_minus_tt(double jd_tt);

/* Convert between Universal Time and Terrestrial Time.  The argument is a
 * calendar year with fraction (used only to evaluate Delta-T). */
double eph_ut_to_tt(double jd_ut, double year);
double eph_tt_to_ut(double jd_tt, double year);

/* Civil UTC to TT and UT1 Julian dates.  From 1972 TT is exact (leap
 * second table; TT = TAI + 32.184 s) and UT1 comes from Delta-T (within
 * ~0.1 s); before 1972 UTC is taken as UT.  Returns -1 for an invalid date.
 * Either output pointer may be NULL. */
int eph_utc_to_jd(int year, int month, int day, int hour, int minute,
                  double second, double* jd_tt, double* jd_ut);

/* Greenwich apparent sidereal time (radians), IAU 2006.  Local sidereal
 * time = this + east longitude. */
double eph_sidereal_time(double jd_ut);

/* Equation of time (days): apparent minus mean solar time.  Local apparent
 * time = local mean time + e.  Needs the Sun's data; returns -1 without. */
int eph_equation_of_time(double jd_ut, double* e_days);

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

typedef enum {
    EPH_AYANAMSA_LAHIRI = 0,
    EPH_AYANAMSA_FAGAN_BRADLEY,
    EPH_AYANAMSA_RAMAN,
    EPH_AYANAMSA_KRISHNAMURTI,
    EPH_AYANAMSA_YUKTESHWAR,
    EPH_AYANAMSA_DE_LUCE,
    EPH_AYANAMSA_USHA_SHASHI,
    EPH_AYANAMSA_DJWHAL_KHUL,
    EPH_AYANAMSA_JN_BHASIN,
    EPH_AYANAMSA_BABYLONIAN_KUGLER1,
    EPH_AYANAMSA_BABYLONIAN_KUGLER2,
    EPH_AYANAMSA_BABYLONIAN_KUGLER3,
    EPH_AYANAMSA_BABYLONIAN_HUBER,
    EPH_AYANAMSA_BABYLONIAN_ETA_PISCIUM,
    EPH_AYANAMSA_BABYLONIAN_ALDEBARAN,
    EPH_AYANAMSA_HIPPARCHOS,
    EPH_AYANAMSA_SASSANIAN,
    EPH_AYANAMSA_J2000,
    EPH_AYANAMSA_J1900,
    EPH_AYANAMSA_B1950,
    EPH_AYANAMSA_SURYASIDDHANTA,
    EPH_AYANAMSA_SURYASIDDHANTA_MEAN_SUN,
    EPH_AYANAMSA_ARYABHATA,
    EPH_AYANAMSA_ARYABHATA_MEAN_SUN,
    EPH_AYANAMSA_SS_REVATI,
    EPH_AYANAMSA_SS_CITRA,
    EPH_AYANAMSA_SKYDRAM,
    EPH_AYANAMSA_ARYABHATA_522,
    EPH_AYANAMSA_BABYLONIAN_BRITTON,
    EPH_AYANAMSA_GALACTIC_EQUATOR_FIORENZA,
    EPH_AYANAMSA_VETTIUS_VALENS,
    EPH_AYANAMSA_LAHIRI_1940,
    EPH_AYANAMSA_LAHIRI_VP285,
    EPH_AYANAMSA_KRISHNAMURTI_SENTHILATHIBAN,
    EPH_AYANAMSA_LAHIRI_ICRC,
    EPH_AYANAMSA_COUNT
} eph_ayanamsa_system;

/* Mean ayanamsa (radians): sidereal longitude = tropical - ayanamsa.
 * NaN for an unknown system.  Needs no data file.  Each system is its
 * defining value carried by the IAU-2006 general precession; this matches
 * Swiss Ephemeris within 0.01" for the modern definitions and 0.3" (over
 * 1600-2500) for those fixed at an ancient epoch.  Star-based systems
 * (True Citra, galactic centre, ...) are not included. */
double eph_ayanamsa(eph_ayanamsa_system system, double jd_tt);

/* Display name of a system ("Lahiri", ...), or NULL. */
const char* eph_ayanamsa_name(eph_ayanamsa_system system);

/* A user-defined ayanamsa: value_t0 (radians) at jd_t0, carried by the
 * general precession to jd_tt. */
double eph_ayanamsa_custom(double jd_t0, double value_t0, double jd_tt);

/* eph_ayanamsa(EPH_AYANAMSA_LAHIRI, jd_tt). */
double eph_ayanamsa_lahiri(double jd_tt);

/* Ascendant and Midheaven ecliptic longitudes (radians, tropical). */
int eph_asc_mc(double jd_ut, const eph_observer* observer, double* asc,
               double* mc);

typedef enum {
    EPH_HOUSES_WHOLE_SIGN = 0,
    EPH_HOUSES_EQUAL,
    EPH_HOUSES_PORPHYRY,
    EPH_HOUSES_PLACIDUS,
    EPH_HOUSES_KOCH,
    EPH_HOUSES_REGIOMONTANUS,
    EPH_HOUSES_CAMPANUS,
    EPH_HOUSES_ALCABITIUS,
    EPH_HOUSES_TOPOCENTRIC, /* Polich-Page */
    EPH_HOUSES_VEHLOW,      /* equal, Ascendant in the middle of house 1   */
    EPH_HOUSES_EQUAL_MC,    /* equal from the MC (house 10)                */
    EPH_HOUSES_EQUAL_ARIES, /* equal from 0 Aries                          */
    EPH_HOUSES_MERIDIAN,    /* axial rotation: equal steps of RA from RAMC */
    EPH_HOUSES_MORINUS      /* equator steps, projected through the
                               ecliptic poles                              */
} eph_house_system;

/* House cusps 1..12 as ecliptic longitudes (radians, tropical).  Returns
 * 0, -1 on bad arguments, or 1 where the system is undefined (Placidus,
 * Koch and Alcabitius near the poles, when a needed point never rises or
 * sets): the cusps are then Porphyry's.  Needs no data file. */
int eph_houses(double jd_ut, const eph_observer* observer,
               eph_house_system system, double cusps[12]);

/* The chart's angles and sensitive points (ecliptic longitudes, radians,
 * except armc).  Needs no data file. */
typedef struct {
    double ascendant;
    double mc;
    double armc;                  /* right ascension of the MC (sidereal time) */
    double vertex;                /* ecliptic x prime vertical, in the west    */
    double equatorial_ascendant;  /* "East point": ecliptic point at RA RAMC+90 */
    double co_ascendant_koch;     /* W. Koch                                   */
    double co_ascendant_munkasey; /* M. Munkasey                               */
    double polar_ascendant;       /* M. Munkasey                               */
} eph_angles;

int eph_house_angles(double jd_ut, const eph_observer* observer,
                     eph_angles* out);

/* House position of an ecliptic longitude (radians) in the given cusps:
 * 1.0 at the start of house 1 up to 12.999..., linear in longitude within
 * each house.  (Swiss Ephemeris' swe_house_pos also uses the latitude for
 * some systems; this uses the longitude only.) */
double eph_house_position(const double cusps[12], double longitude);

const char* eph_body_name(eph_body body);

/* ---- phenomena ------------------------------------------------------- */

typedef struct {
    double phase_angle; /* Sun-body-Earth angle, radians                     */
    double phase;       /* illuminated fraction of the disc, 0..1            */
    double elongation;  /* angle from the Sun as seen from the Earth, radians */
    double diameter;    /* apparent diameter, radians (mean radius)           */
    double magnitude;   /* visual magnitude (Mallama & Hilton 2018 for the
                           planets, H-G system for the minor bodies)         */
} eph_phenomena;

/* Apparent phenomena of the Sun, Moon, planets, Pluto and minor bodies
 * (not the lunar points).  Returns -1 without the needed data. */
int eph_phenomena_tt(eph_body body, double jd_tt, eph_phenomena* out);
int eph_phenomena_ut(eph_body body, double jd_ut, double year, eph_phenomena* out);

/* ---- event searches -------------------------------------------------- *
 * Search forward (direction +1) or backward (-1) from a date for at most
 * max_days.  Return -1 if there is no event in that window or a position
 * cannot be computed.  Times are TT unless the name ends in _ut; results
 * are good to ~1e-8 day of the positions. */

/* When the body's apparent geocentric longitude reaches `longitude`
 * (radians): returns (Sun: solar return), ingresses, transits to a point.
 * Returns 0.  Retrograde loops are followed (the first crossing is found);
 * a target the body only touches at a station is not a crossing. */
int eph_longitude_crossing(eph_body body, double longitude, double jd_tt,
                           int direction, double max_days, double* jd_out);
int eph_longitude_crossing_ut(eph_body body, double longitude, double jd_ut,
                              int direction, double max_days, double* jd_out);

/* The same for the heliocentric longitude (planets, Pluto, asteroids). */
int eph_heliocentric_crossing(eph_body body, double longitude, double jd_tt,
                              int direction, double max_days, double* jd_out);

/* When the body crosses the ecliptic (latitude 0), e.g. the Moon at its
 * node.  Returns 1 at the ascending node, 2 at the descending node; the
 * longitude there is stored if `longitude` is not NULL. */
int eph_node_crossing(eph_body body, double jd_tt, int direction,
                      double max_days, double* jd_out, double* longitude);

/* When the body's longitude speed changes sign.  Returns 1 if it turns
 * retrograde, 2 if it turns direct. */
int eph_station(eph_body body, double jd_tt, int direction, double max_days,
                double* jd_out);

/* ---- utilities -------------------------------------------------------- */

/* Midpoint (radians) of the shorter arc between two longitudes. */
double eph_midpoint(double a, double b);

/* Split a longitude (radians) into zodiac sign (0 = Aries .. 11), degrees,
 * minutes and seconds within the sign.  Any output pointer may be NULL. */
void eph_split_longitude(double longitude, int* sign, int* degrees,
                         int* minutes, double* seconds);

/* ---- data files ------------------------------------------------------ *
 * Positions need series data from a file made by tools/pack_series or
 * tools/pack_ephemeris.py.  Pick one series level (fitted over 1550-2650):
 *   data/eph-1as.bin    ~0.5" (truncated series)
 *   data/eph-01as.bin   ~0.2" (truncated series)
 *   data/eph-full.bin   the full series (~0.1")
 * (each also as eph-*-chiron.bin, with Chiron), and optionally add
 * data/eph-chiron.bin (Chiron) or data/eph-asteroids.bin (Ceres, Pallas,
 * Juno, Vesta, Chiron; 1600-2500) and/or DE440 tables data/de440/de440-YEAR.bin (one century
 * each, ~0.002"; self-contained: they work without a series file).
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
