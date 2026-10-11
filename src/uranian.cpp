#include "internal.h"

#include "../include/ephemeris/ephemeris.h"

/* The Hamburg school's hypothetical "Uranian" points: Cupido, Hades, Zeus
 * and Kronos (A. Witte, 1923-24), Apollon, Admetos, Vulkanus and Poseidon
 * (F. Sieggrun, 1927 and 1934).  Their orbits are the elements P. Schlyter
 * derived from the published Witte/Sieggrun ephemeris ("Hypothetical
 * astrological planets", alt.astrology, 1994): a mean longitude moving
 * uniformly in the tropical zodiac (the rate includes precession), the
 * first term of the equation of centre for the four inner points, and the
 * inclination of Cupido and Hades; Schlyter found the ephemeris fitted to a
 * few arcminutes.  Checked here: the mean longitudes and rates (1900.0 to
 * 1950.0) and the distances (Kepler's third law); the e, perihelion and
 * node terms (< ~0.5 deg in position) are as transcribed, unverified.  The points are then seen from the Earth like a body. */

namespace eph {

namespace {

struct Uranian {
    double lon1900;  /* mean longitude at 1900.0 (JD 2415020.0), degrees */
    double rate;     /* degrees per Julian century, tropical             */
    double a;        /* AU                                               */
    double e;
    double peri;     /* longitude of perihelion, degrees                 */
    double inc;      /* degrees                                          */
    double node;     /* longitude of the ascending node, degrees         */
};

const Uranian kUranian[8] = {
    {105.026, 138.535, 40.99, 0.0045, 300.0, 1.08, 130.0}, /* Cupido   */
    {337.180, 101.219, 50.67, 0.0024, 305.0, 1.05, 342.0}, /* Hades    */
    {104.172, 80.406, 59.21, 0.0011, 295.0, 0.0, 0.0},     /* Zeus     */
    {17.298, 70.394, 64.81, 0.0028, 205.0, 0.0, 0.0},      /* Kronos   */
    {138.060, 62.474, 70.3, 0.0, 0.0, 0.0, 0.0},           /* Apollon  */
    {351.334, 58.374, 73.6, 0.0, 0.0, 0.0, 0.0},           /* Admetos  */
    {55.900, 54.407, 77.3, 0.0, 0.0, 0.0, 0.0},            /* Vulkanus */
    {165.521, 48.435, 83.7, 0.0, 0.0, 0.0, 0.0},           /* Poseidon */
};

}  // namespace

bool is_uranian(int body) { return body >= EPH_CUPIDO && body <= EPH_POSEIDON; }

/* Heliocentric position (AU), ecliptic and equinox of date. */
bool uranian_heliocentric(int body, double jd_tt, Vec3& v) {
    if (!is_uranian(body)) return false;
    const Uranian& u = kUranian[body - EPH_CUPIDO];
    const double T = (jd_tt - 2415020.0) / DAYS_PER_CENTURY;
    const double L = (u.lon1900 + u.rate * T) * DEG;
    const double M = L - u.peri * DEG;
    const double lon = L + 2.0 * u.e * std::sin(M);
    const double lat = u.inc * DEG * std::sin(lon - u.node * DEG);
    v = spherical_to_cartesian(lon, lat, u.a * (1.0 - u.e * std::cos(M)));
    return true;
}

/* Geocentric position: geometric, or apparent (aberration with the Earth's
 * velocity and nutation, as for the stars; the light-time, ~6 h, moves
 * them by < 0.1"). */
bool uranian_position(int body, double jd_tt, bool apparent, double& lon, double& lat,
                      double& dist) {
    Vec3 h, e, ev;
    if (!uranian_heliocentric(body, jd_tt, h) || !earth_heliocentric_state(jd_tt, e, ev))
        return false;
    const Vec3 g = {h.x - e.x, h.y - e.y, h.z - e.z};
    dist = length(g);
    if (apparent) {
        reduce_direction(jd_tt, true, g, e, ev, lon, lat);
    } else {
        cartesian_to_spherical(frame_1976_to_2006(jd_tt, g), lon, lat, dist);
    }
    return true;
}

}  // namespace eph
