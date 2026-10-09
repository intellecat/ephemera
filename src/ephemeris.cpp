#include "internal.h"

#include "../include/ephemeris/ephemeris.h"

namespace eph {
namespace {

Vec3 planet_heliocentric(int planet, double jd) {
    double l, b, r;
    vsop87_heliocentric(planet, jd, l, b, r);
    return spherical_to_cartesian(l, b, r);
}

Vec3 pluto_heliocentric(double jd) {
    double l, b, r;
    pluto_heliocentric_j2000(jd, l, b, r);
    precess_ecliptic_j2000_to_date(jd, l, b);
    return spherical_to_cartesian(l, b, r);
}

Vec3 j2000_to_date(double jd, const Vec3& v) {
    double l, b, r;
    cartesian_to_spherical(v, l, b, r);
    precess_ecliptic_j2000_to_date(jd, l, b);
    return spherical_to_cartesian(l, b, r);
}

bool body_heliocentric(eph_body body, double jd, Vec3& out) {
    if (body == EPH_PLUTO) {
        out = pluto_heliocentric(jd);
        return true;
    }
    if (body >= EPH_CERES && body <= EPH_VESTA) {
        Vec3 v;
        if (!asteroid_heliocentric((int)body - EPH_CERES, jd, v)) return false;
        out = j2000_to_date(jd, v);
        return true;
    }
    out = planet_heliocentric(static_cast<int>(body), jd);
    return true;
}

Vec3 normalize_vec(const Vec3& v) {
    const double n = length(v);
    return (n > 0.0) ? Vec3{v.x / n, v.y / n, v.z / n} : v;
}

Vec3 add(const Vec3& a, const Vec3& b) { return {a.x + b.x, a.y + b.y, a.z + b.z}; }

Vec3 scale(const Vec3& v, double k) { return {v.x * k, v.y * k, v.z * k}; }

double dot(const Vec3& a, const Vec3& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

double angle_diff(double a, double b) { return normalize(a - b + PI) - PI; }

/* Position of a body in the IAU-1976 ecliptic of date: heliocentric for the
 * planets and the Earth, geocentric for the Moon, zero for the Sun.  From
 * the compressed table when `table` is set, else the analytical series. */
bool source(eph_body body, double jd, bool table, Vec3& v) {
    if (body == EPH_SUN) {
        v = {0.0, 0.0, 0.0};
        return true;
    }
    if (table) {
        Vec3 j;
        if (!compressed_raw((int)body, jd, j)) return false;
        v = j2000_to_date(jd, j);
        return true;
    }
    if (body == EPH_MOON) {
        double l, b, r;
        moon_position(jd, l, b, r);
        v = spherical_to_cartesian(l, b, r);
        return true;
    }
    return body_heliocentric(body, jd, v);
}

/* Position, velocity and acceleration from three evaluations. */
struct State {
    Vec3 p, v, a;
};

bool state(eph_body body, double jd, bool table, State& s) {
    /* Small enough for the fast inner bodies, large enough not to amplify
     * rounding in the slow outer ones (Pluto's series cancels heavily). */
    const double h = (body == EPH_MOON)                         ? 0.001
                     : (body >= EPH_JUPITER && body <= EPH_PLUTO) ||
                               body >= EPH_CERES
                         ? 0.1
                         : 0.01;
    Vec3 m, c, pl;
    if (table && body != EPH_SUN) {
        /* Table segments join with small steps, so take the derivatives
         * from the segment containing jd and expand it to jd +- h. */
        Vec3 p, v, a;
        if (!compressed_raw((int)body, jd, p, &v, &a)) return false;
        m = j2000_to_date(jd - h, add(p - scale(v, h), scale(a, 0.5 * h * h)));
        c = j2000_to_date(jd, p);
        pl = j2000_to_date(jd + h, add(add(p, scale(v, h)), scale(a, 0.5 * h * h)));
    } else if (!source(body, jd - h, table, m) || !source(body, jd, table, c) ||
               !source(body, jd + h, table, pl)) {
        return false;
    }
    s.p = c;
    s.v = scale(pl - m, 1.0 / (2.0 * h));
    s.a = scale(add(pl, m) - scale(c, 2.0), 1.0 / (h * h));
    return true;
}

/* The Earth is needed for every body at the same instant; keep the last one
 * (per thread, keyed on the table generation). */
struct EarthCache {
    double jd = -1.0;
    int generation = -1;
    State s;
};
thread_local EarthCache g_earth;

bool earth_state(double jd, bool table, State& s) {
    const int gen = table ? compressed_generation() : 0;
    if (g_earth.jd == jd && g_earth.generation == gen) {
        s = g_earth.s;
        return true;
    }
    if (!state(EPH_EARTH, jd, table, s)) return false;
    g_earth.jd = jd;
    g_earth.generation = gen;
    g_earth.s = s;
    return true;
}

/* Position and velocity at jd - tau from the expansion about jd; the
 * neglected third-order term is < 1e-10 AU for the light-times involved. */
void retard(const State& s, double tau, Vec3& p, Vec3& v) {
    p = add(s.p - scale(s.v, tau), scale(s.a, 0.5 * tau * tau));
    v = s.v - scale(s.a, tau);
}

/* Rotate about the ecliptic pole (longitude += a). */
Vec3 rot_z(const Vec3& v, double a) {
    const double c = std::cos(a), s = std::sin(a);
    return {c * v.x - s * v.y, s * v.x + c * v.y, v.z};
}

/* General precession in longitude, IAU 1976, radians per day. */
constexpr double kPrecessionRate = 5029.0966 * ARCSEC / DAYS_PER_CENTURY;

/* Geometric geocentric position and velocity, optionally retarded for
 * light-time.  The Moon's source is already geocentric. */
void geocentric(eph_body body, const State& b, const State& e,
                bool light_time, Vec3& g, Vec3& gv) {
    const bool geo = (body == EPH_MOON);
    const Vec3 e0 = geo ? Vec3{0.0, 0.0, 0.0} : e.p;
    const Vec3 e1 = geo ? Vec3{0.0, 0.0, 0.0} : e.v;
    Vec3 p = b.p, v = b.v;
    g = p - e0;
    if (light_time) {
        for (int i = 0; i < 3; ++i) {
            const double tau = length(g) / AU_PER_LIGHT_DAY;
            retard(b, tau, p, v);
            /* The retarded position is in the ecliptic of date t - tau;
             * precess it to the observer's date t (0.03" for Pluto). */
            p = rot_z(p, kPrecessionRate * tau);
            v = rot_z(v, kPrecessionRate * tau);
            g = p - e0;
        }
        /* d(tau)/dt = u.(dg/dt)/c with dg/dt = v (1 - dtau/dt) - e1. */
        const Vec3 u = normalize_vec(g);
        const double tau_dot =
            dot(u, v - e1) / (AU_PER_LIGHT_DAY + dot(u, v));
        gv = scale(v, 1.0 - tau_dot) - e1;
    } else {
        gv = v - e1;
    }
}

enum class Kind { Geometric, Apparent };

/* Relativistic light deflection by the Sun for a body at finite distance
 * (Explanatory Supplement):
 *   p' = p + 2GM/(c^2 E) [(p.q) e - (e.p) q] / (1 + q.e)
 * p: geocentric, q: heliocentric body, e: heliocentric Earth (unit vectors).
 * Far behind the Sun this is 2GM/(c^2 E) cot(psi/2); it vanishes for a body
 * in front of the Sun (the Moon, inferior conjunctions). */
Vec3 deflect(const Vec3& p, const Vec3& q, const Vec3& e, double earth_sun_au) {
    const double denom = 1.0 + dot(q, e);
    if (denom < 1e-9) return p;
    const double gm = 1.32712440018e20; /* m^3/s^2 */
    const double c = 2.99792458e8;      /* m/s     */
    const double f = 2.0 * gm / (c * c * earth_sun_au * 1.495978707e11);
    const double pq = dot(p, q), ep = dot(e, p);
    return normalize_vec(Vec3{p.x + f * (pq * e.x - ep * q.x) / denom,
                              p.y + f * (pq * e.y - ep * q.y) / denom,
                              p.z + f * (pq * e.z - ep * q.z) / denom});
}

/* Geocentric vector -> ecliptic longitude/latitude of date (IAU-2006).
 * Apparent adds deflection, annual aberration (Earth velocity ev) and
 * nutation; the Moon gets no annual aberration (it moves with the Earth). */
void reduce(eph_body body, double jd, bool apparent, const Vec3& g,
            const Vec3& e, const Vec3& ev, double& lon, double& lat) {
    Vec3 u = normalize_vec(g);
    if (apparent) {
        if (body != EPH_SUN) {
            const double em = length(e);
            u = deflect(u, normalize_vec(add(e, g)), scale(e, 1.0 / em), em);
        }
        if (body != EPH_MOON)
            u = normalize_vec(add(u, scale(ev, 1.0 / AU_PER_LIGHT_DAY)));
    }
    double r;
    cartesian_to_spherical(frame_1976_to_2006(jd, u), lon, lat, r);
    if (apparent) {
        double dpsi, deps;
        nutation(jd, dpsi, deps);
        lon = normalize(lon + dpsi);
    }
}

bool compute(eph_body body, double jd_tt, Kind kind, eph_position* out) {
    if (body < 0 || body >= EPH_BODY_COUNT || body == EPH_EARTH) return false;
    const bool apparent = (kind == Kind::Apparent);
    /* Light-time applies to everything except the Sun (at the origin). */
    const bool light_time = apparent && body != EPH_SUN;

    /* The compressed table when it covers this date, else the analytical
     * series; bodies the table lacks (asteroids) fall back individually. */
    const bool table = eph_compressed_available() && compressed_covers(jd_tt);
    State e, b;
    if (!earth_state(jd_tt, table, e) && !(table && earth_state(jd_tt, false, e)))
        return false;
    if (!state(body, jd_tt, table, b) && !(table && state(body, jd_tt, false, b)))
        return false;

    Vec3 g, gv;
    geocentric(body, b, e, light_time, g, gv);
    reduce(body, jd_tt, apparent, g, e.p, e.v, out->longitude, out->latitude);
    out->distance = length(g);

    /* Speeds: central difference of the reduction along the velocities
     * (exact to first order, so the step only has to be small). */
    const double d = 1e-3;
    double l1, b1, l2, b2;
    reduce(body, jd_tt - d, apparent, g - scale(gv, d), e.p - scale(e.v, d),
           e.v - scale(e.a, d), l1, b1);
    reduce(body, jd_tt + d, apparent, add(g, scale(gv, d)),
           add(e.p, scale(e.v, d)), add(e.v, scale(e.a, d)), l2, b2);
    out->longitude_speed = angle_diff(l2, l1) / (2 * d);
    out->latitude_speed = (b2 - b1) / (2 * d);
    out->distance_speed = dot(normalize_vec(g), gv);
    return true;
}

/* Diurnal parallax: geocentric apparent RA/Dec -> topocentric. */
void topocentric_correct(double ra, double dec, double dist_au, double jd_ut,
                         const eph_observer* obs, double* ra2, double* dec2,
                         double* dist2) {
    const double phi = obs->latitude;
    const double u = std::atan(0.99664719 * std::tan(phi));
    const double rs =
        0.99664719 * std::sin(u) + obs->elevation / 6378140.0 * std::sin(phi);
    const double rc =
        std::cos(u) + obs->elevation / 6378140.0 * std::cos(phi);
    const double sinpi = 6378.14 / (dist_au * AU_KM);
    const double theta = gast_degrees(jd_ut) * DEG + obs->longitude;
    const double H = theta - ra;
    const double da = std::atan2(-rc * sinpi * std::sin(H),
                                 std::cos(dec) - rc * sinpi * std::cos(H));
    *ra2 = ra + da;
    *dec2 = std::atan2((std::sin(dec) - rs * sinpi) * std::cos(da),
                       std::cos(dec) - rc * sinpi * std::cos(H));

    /* Topocentric distance: body minus the observer's geocentric vector
     * (Earth equatorial radius 6378.14 km, in AU). */
    const double re = 6378.14 / AU_KM;
    const double xo = rc * re * std::cos(theta), yo = rc * re * std::sin(theta);
    const double zo = rs * re;
    const double xb = dist_au * std::cos(dec) * std::cos(ra);
    const double yb = dist_au * std::cos(dec) * std::sin(ra);
    const double zb = dist_au * std::sin(dec);
    const double dx = xb - xo, dy = yb - yo, dz = zb - zo;
    *dist2 = std::sqrt(dx * dx + dy * dy + dz * dz);
}

bool topocentric_eval(eph_body body, double jd_tt, double jd_ut,
                      const eph_observer* obs, double* lon, double* lat,
                      double* dist) {
    eph_position g;
    if (eph_position_apparent_tt(body, jd_tt, &g) != 0) return false;
    double ra, dec;
    ecliptic_to_equatorial(g.longitude, g.latitude, jd_tt, &ra, &dec);
    double ra2, dec2;
    topocentric_correct(ra, dec, g.distance, jd_ut, obs, &ra2, &dec2, dist);
    equatorial_to_ecliptic(ra2, dec2, jd_tt, lon, lat);
    return true;
}

}  // namespace

extern "C" {

int eph_position_tt(eph_body body, double jd_tt, eph_position* out) {
    if (!out) return -1;
    return compute(body, jd_tt, Kind::Geometric, out) ? 0 : -1;
}

int eph_position_apparent_tt(eph_body body, double jd_tt, eph_position* out) {
    if (!out) return -1;
    return compute(body, jd_tt, Kind::Apparent, out) ? 0 : -1;
}

int eph_position_ut(eph_body body, double jd_ut, double year,
                    eph_position* out) {
    return eph_position_tt(body, eph_ut_to_tt(jd_ut, year), out);
}

int eph_position_apparent_ut(eph_body body, double jd_ut, double year,
                             eph_position* out) {
    return eph_position_apparent_tt(body, eph_ut_to_tt(jd_ut, year), out);
}

int eph_heliocentric_tt(eph_body body, double jd_tt, double* lon, double* lat,
                        double* distance) {
    const bool planet = (body >= EPH_MERCURY && body <= EPH_PLUTO);
    const bool asteroid = (body >= EPH_CERES && body <= EPH_VESTA);
    if (!planet && !asteroid) return -1;
    Vec3 h;
    if (!body_heliocentric(body, jd_tt, h)) return -1;
    const Vec3 v = frame_1976_to_2006(jd_tt, h);
    double r;
    double l, b;
    cartesian_to_spherical(v, l, b, r);
    if (lon) *lon = l;
    if (lat) *lat = b;
    if (distance) *distance = r;
    return 0;
}

int eph_position_topocentric_tt(eph_body body, double jd_tt, double jd_ut,
                                const eph_observer* observer,
                                eph_position* out) {
    if (!observer || !out) return -1;
    double lon, lat, dist;
    if (!topocentric_eval(body, jd_tt, jd_ut, observer, &lon, &lat, &dist))
        return -1;
    out->longitude = lon;
    out->latitude = lat;
    out->distance = dist;

    const double dt = (body == EPH_MOON) ? 0.001 : 0.1;
    double l1, b1, d1, l2, b2, d2;
    if (topocentric_eval(body, jd_tt - dt, jd_ut - dt, observer, &l1, &b1, &d1) &&
        topocentric_eval(body, jd_tt + dt, jd_ut + dt, observer, &l2, &b2, &d2)) {
        out->longitude_speed = angle_diff(l2, l1) / (2 * dt);
        out->latitude_speed = (b2 - b1) / (2 * dt);
        out->distance_speed = (d2 - d1) / (2 * dt);
    } else {
        out->longitude_speed = out->latitude_speed = out->distance_speed = 0.0;
    }
    return 0;
}

int eph_position_topocentric_ut(eph_body body, double jd_ut, double year,
                                const eph_observer* observer,
                                eph_position* out) {
    return eph_position_topocentric_tt(body, eph_ut_to_tt(jd_ut, year), jd_ut,
                                       observer, out);
}

}  // extern "C"

}  // namespace eph
