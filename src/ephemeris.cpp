#include "internal.h"

#include "../include/ephemeris/ephemeris.h"
#include "data.h"

namespace eph {
namespace {

bool planet_heliocentric(int planet, double jd, Vec3& out) {
    double l, b, r;
    if (!vsop87_heliocentric(planet, jd, l, b, r)) return false;
    out = spherical_to_cartesian(l, b, r);
    return true;
}

bool pluto_heliocentric(double jd, double anchor, Vec3& out) {
    double l, b, r;
    if (!pluto_heliocentric_j2000(jd, anchor, l, b, r)) return false;
    precess_ecliptic_j2000_to_date(jd, l, b);
    out = spherical_to_cartesian(l, b, r);
    return true;
}

Vec3 j2000_to_date(double jd, const Vec3& v) {
    double l, b, r;
    cartesian_to_spherical(v, l, b, r);
    precess_ecliptic_j2000_to_date(jd, l, b);
    return spherical_to_cartesian(l, b, r);
}

/* anchor: the date whose Chebyshev interval is used for Pluto and the
 * asteroids (see chebyshev_body). */
bool body_heliocentric(eph_body body, double jd, Vec3& out, double anchor) {
    if (body == EPH_PLUTO) return pluto_heliocentric(jd, anchor, out);
    if (body >= EPH_CERES && body <= EPH_CHIRON) {
        Vec3 v;
        if (!asteroid_heliocentric((int)body - EPH_CERES, jd, anchor, v))
            return false;
        out = j2000_to_date(jd, v);
        return true;
    }
    if (is_minor(body)) {
        Vec3 v;
        if (!minor_heliocentric(body, jd, v)) return false;
        out = j2000_to_date(jd, v);
        return true;
    }
    return planet_heliocentric(static_cast<int>(body), jd, out);
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

/* Position of a body in the IAU-1976 ecliptic of date from the analytical
 * series: heliocentric for the planets and the Earth, geocentric for the
 * Moon, zero for the Sun. */
bool source(eph_body body, double jd, double anchor, Vec3& v) {
    if (body == EPH_SUN) {
        v = {0.0, 0.0, 0.0};
        return true;
    }
    if (body == EPH_MOON) {
        double l, b, r;
        if (!moon_position(jd, l, b, r)) return false;
        v = spherical_to_cartesian(l, b, r);
        return true;
    }
    return body_heliocentric(body, jd, v, anchor);
}

/* The truncated series under the residual table, same frame. */
bool base(eph_body body, double jd, Vec3& v) {
    if (body == EPH_MOON) {
        double x, y, z;
        if (!elp_base_j2000(jd, x, y, z)) return false;
        v = j2000_to_date(jd, {x / AU_KM, y / AU_KM, z / AU_KM});
        return true;
    }
    double l, b, r;
    if (!vsop87_base((int)body, jd, l, b, r)) return false;
    v = spherical_to_cartesian(l, b, r);
    return true;
}

/* Position, velocity and acceleration. */
struct State {
    Vec3 p, v, a;
};

/* Finite differences of a function over jd - h, jd, jd + h. */
void differentiate(const Vec3& m, const Vec3& c, const Vec3& pl, double h,
                   State& s) {
    s.p = c;
    s.v = scale(pl - m, 1.0 / (2.0 * h));
    s.a = scale(add(pl, m) - scale(c, 2.0), 1.0 / (h * h));
}

/* Moon/Earth mass ratio + 1 (DE440). */
constexpr double kEarthMoonRatio1 = 82.3005682214972154;

bool state(eph_body body, double jd, bool table, State& s);

/* The Earth and the Moon are needed for every body at the same instant;
 * keep the last ones (per thread, keyed on the loaded data). */
struct Cache {
    double jd = -1.0;
    int generation = -1;
    bool table = false;
    State s;
};
thread_local Cache g_earth, g_moon;

bool cached_state(eph_body body, double jd, bool table, State& s) {
    Cache& c = (body == EPH_EARTH) ? g_earth : g_moon;
    const int gen = data_generation();
    if (c.jd == jd && c.generation == gen && c.table == table) {
        s = c.s;
        return true;
    }
    if (!state(body, jd, table, s)) return false;
    c.jd = jd;
    c.generation = gen;
    c.table = table;
    c.s = s;
    return true;
}

bool state(eph_body body, double jd, bool table, State& s) {
    /* Small enough for the fast inner bodies, large enough not to amplify
     * rounding in the slow outer ones (Pluto's series cancels heavily). */
    const double h = (body == EPH_MOON)                         ? 0.001
                     : (body >= EPH_JUPITER && body <= EPH_PLUTO) ||
                               body >= EPH_CERES
                         ? 0.1
                         : 0.01;
    if (table && body != EPH_SUN) {
        /* Residual derivatives come from its own Chebyshev series. */
        bool has_base;
        Vec3 p, v, a;
        if (!compressed_residual((int)body, jd, p, v, a, has_base))
            return false;
        s = {p, v, a};
        if (has_base) {
            Vec3 m, c, pl;
            if (!base(body, jd - h, m) || !base(body, jd, c) ||
                !base(body, jd + h, pl))
                return false;
            State b;
            differentiate(m, c, pl, h, b);
            s = {add(s.p, b.p), add(s.v, b.v), add(s.a, b.a)};
        }
        if (body == EPH_EARTH) {
            /* The table holds the Earth-Moon barycentre. */
            State m;
            if (!cached_state(EPH_MOON, jd, true, m)) return false;
            const double k = 1.0 / kEarthMoonRatio1;
            s = {s.p - scale(m.p, k), s.v - scale(m.v, k), s.a - scale(m.a, k)};
        }
        return true;
    }
    Vec3 m, c, pl;
    if (!source(body, jd - h, jd, m) || !source(body, jd, jd, c) ||
        !source(body, jd + h, jd, pl))
        return false;
    differentiate(m, c, pl, h, s);
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

Vec3 cross(const Vec3& a, const Vec3& b) {
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}

/* Longitude and latitude, and their rates, of the direction u (rate ud),
 * both in the IAU-2006 ecliptic of date. */
void direction(const Vec3& u, const Vec3& ud, eph_position* out) {
    double r;
    cartesian_to_spherical(u, out->longitude, out->latitude, r);
    const double rho2 = u.x * u.x + u.y * u.y;
    out->longitude_speed = (u.x * ud.y - u.y * ud.x) / rho2;
    out->latitude_speed = (ud.z * rho2 - u.z * (u.x * ud.x + u.y * ud.y)) /
                          (r * r * std::sqrt(rho2));
    out->distance = 0.0;
    out->distance_speed = 0.0;
}

/*
 * Lunar nodes and apogee ("Black Moon Lilith"), geocentric, ecliptic of
 * date; apparent positions add nutation (no light-time or aberration, as
 * these are not bodies).
 *
 * Mean node and mean perigee: the polynomials of ELP2000-82B as given by
 * Meeus, Astronomical Algorithms (2nd ed.), ch. 47.  The mean apogee lies
 * on the mean lunar orbit (inclination 5.1454 deg), so it has a latitude.
 * True node and osculating apogee: from the osculating (two-body) orbit of
 * the Moon's geometric position r and velocity v, with derivatives from
 * its acceleration: the node is along z x (r x v), the apogee opposite the
 * eccentricity vector (v x h) / mu - r / |r|.
 */
/* The Moon's geometric position alone, from the table or the series. */
bool moon_at(double jd, bool table, Vec3& p) {
    if (table) {
        Vec3 v, a;
        bool has_base;
        if (!compressed_residual(EPH_MOON, jd, p, v, a, has_base)) return false;
        if (!has_base) return true;
        Vec3 b;
        if (!base(EPH_MOON, jd, b)) return false;
        p = add(p, b);
        return true;
    }
    return source(EPH_MOON, jd, jd, p);
}

/* The Moon's state for the osculating elements: the acceleration from
 * +-0.02 day, not the +-0.001 day of state(), where the ~5e-10 day
 * resolution of a double Julian date leaves ~0.2% noise.  The node and
 * apogee rates depend only on the small (~0.6%) perturbing part of the
 * acceleration.  Cached per date and thread. */
struct MoonOsc {
    double jd = -1.0;
    int generation = -1;
    bool table = false;
    State s;
};
thread_local MoonOsc g_moon_osc;

bool moon_osculating_state(double jd, State& s) {
    bool table = eph_compressed_available() && compressed_covers(jd);
    MoonOsc& c = g_moon_osc;
    if (c.jd == jd && c.generation == data_generation() && c.table == table) {
        s = c.s;
        return true;
    }
    const double H = 0.02;
    Vec3 m, pl;
    if (!cached_state(EPH_MOON, jd, table, s) || !moon_at(jd - H, table, m) ||
        !moon_at(jd + H, table, pl)) {
        if (!table) return false;
        table = false;
        if (!cached_state(EPH_MOON, jd, false, s) || !moon_at(jd - H, false, m) ||
            !moon_at(jd + H, false, pl))
            return false;
    }
    s.a = scale(add(pl, m) - scale(s.p, 2.0), 1.0 / (H * H));
    c = {jd, data_generation(), table, s};
    return true;
}

bool lunar_point(eph_body body, double jd, bool apparent, eph_position* out) {
    const double T = (jd - J2000) / DAYS_PER_CENTURY;
    if (body == EPH_MEAN_NODE || body == EPH_MEAN_LILITH) {
        const double node =
            (125.0445479 + T * (-1934.1362891 + T * (0.0020754 +
                                T * (1.0 / 467441.0 - T / 60616000.0)))) * DEG;
        const double node_rate =
            (-1934.1362891 + T * (2 * 0.0020754 + T * (3.0 / 467441.0 -
                                  T * 4.0 / 60616000.0))) * DEG / DAYS_PER_CENTURY;
        if (body == EPH_MEAN_NODE) {
            out->longitude = normalize(node);
            out->latitude = 0.0;
            out->longitude_speed = node_rate;
            out->latitude_speed = 0.0;
            out->distance = out->distance_speed = 0.0;
        } else {
            const double perigee =
                (83.3532465 + T * (4069.0137287 + T * (-0.0103200 +
                                   T * (-1.0 / 80053.0 + T / 18999000.0)))) * DEG;
            const double perigee_rate =
                (4069.0137287 + T * (-2 * 0.0103200 + T * (-3.0 / 80053.0 +
                                     T * 4.0 / 18999000.0))) * DEG / DAYS_PER_CENTURY;
            /* Apogee: argument w from the node along the mean orbit. */
            const double inc = 5.1453964 * DEG;
            const double w = perigee + PI - node, wd = perigee_rate - node_rate;
            const Vec3 q = {std::cos(w), std::sin(w) * std::cos(inc),
                            std::sin(w) * std::sin(inc)};
            const Vec3 qd = {-std::sin(w) * wd, std::cos(w) * std::cos(inc) * wd,
                             std::cos(w) * std::sin(inc) * wd};
            /* Rotate by the node about the ecliptic pole. */
            const double c = std::cos(node), sn = std::sin(node);
            const Vec3 a = {c * q.x - sn * q.y, sn * q.x + c * q.y, q.z};
            const Vec3 ad = {c * qd.x - sn * qd.y - node_rate * a.y,
                             sn * qd.x + c * qd.y + node_rate * a.x, qd.z};
            double r;
            cartesian_to_spherical(a, out->longitude, out->latitude, r);
            out->longitude_speed = (a.x * ad.y - a.y * ad.x) / (a.x * a.x + a.y * a.y);
            out->latitude_speed = ad.z / std::cos(out->latitude);
            out->distance = out->distance_speed = 0.0;
        }
        /* The polynomials refer to the mean ecliptic and equinox of date;
         * the 1976/2006 frame difference is below 0.1" here. */
    } else {
        /* Osculating elements from the Moon's geometric state. */
        State m;
        if (!moon_osculating_state(jd, m)) return false;
        /* The ecliptic of date turns with the precession; osculating
         * elements need the inertial velocity and acceleration (the frame
         * term is ~3e-6 of the velocity, ~20" in the apogee). */
        const Vec3 zv = {-m.v.y, m.v.x, 0.0};
        const Vec3 v = m.v - Vec3{-kPrecessionRate * m.p.y, kPrecessionRate * m.p.x, 0.0};
        const Vec3 acc = m.a - scale(zv, 2.0 * kPrecessionRate);
        m.v = v;
        m.a = acc;
        /* The node is where the orbit crosses the IAU-2006 ecliptic, so
         * the orbit normal goes to that frame first (the 1976 ecliptic is
         * tilted by 0.04", which moves the node by ~0.5"). */
        const Vec3 h = frame_1976_to_2006(jd, cross(m.p, m.v));
        const Vec3 hd = frame_1976_to_2006(jd, cross(m.p, m.a));
        if (body == EPH_TRUE_NODE) {
            direction({-h.y, h.x, 0.0}, {-hd.y, hd.x, 0.0}, out);
            out->latitude = out->latitude_speed = 0.0;
        } else {
            /* G (M_Earth + M_Moon), AU^3/day^2 (DE440). */
            const double mu = 403503.235502 * 86400.0 * 86400.0 /
                              (AU_KM * AU_KM * AU_KM);
            const double r = length(m.p);
            const double rd = dot(m.p, m.v) / r;
            const Vec3 h0 = cross(m.p, m.v), hd0 = cross(m.p, m.a);
            const Vec3 e = add(scale(cross(m.v, h0), 1.0 / mu), scale(m.p, -1.0 / r));
            const Vec3 ed = add(scale(add(cross(m.a, h0), cross(m.v, hd0)), 1.0 / mu),
                                add(scale(m.v, -1.0 / r), scale(m.p, rd / (r * r))));
            direction(frame_1976_to_2006(jd, scale(e, -1.0)),
                      frame_1976_to_2006(jd, scale(ed, -1.0)), out);
        }
        /* Inertial rate -> rate in the precessing ecliptic of date. */
        out->longitude_speed += kPrecessionRate;
    }
    if (apparent) {
        double dpsi, deps, d1, d2;
        nutation(jd, dpsi, deps);
        out->longitude = normalize(out->longitude + dpsi);
        /* Nutation in longitude changes by up to ~0.1"/day. */
        nutation(jd - 0.05, d1, deps);
        nutation(jd + 0.05, d2, deps);
        out->longitude_speed += (d2 - d1) / 0.1;
    }
    return true;
}

bool compute(eph_body body, double jd_tt, Kind kind, eph_position* out) {
    if (body < 0 || (body >= EPH_BODY_COUNT && !is_minor(body)) || body == EPH_EARTH)
        return false;
    if (body >= EPH_MEAN_NODE && body <= EPH_OSCU_LILITH)
        return lunar_point(body, jd_tt, kind == Kind::Apparent, out);
    if (is_uranian(body)) {
        const bool app = kind == Kind::Apparent;
        double l1, b1, r1, l2, b2, r2;
        const double h = 0.5;
        if (!uranian_position(body, jd_tt, app, out->longitude, out->latitude, out->distance) ||
            !uranian_position(body, jd_tt - h, app, l1, b1, r1) ||
            !uranian_position(body, jd_tt + h, app, l2, b2, r2))
            return false;
        out->longitude_speed = angle_diff(l2, l1) / (2 * h);
        out->latitude_speed = (b2 - b1) / (2 * h);
        out->distance_speed = (r2 - r1) / (2 * h);
        return true;
    }
    const bool apparent = (kind == Kind::Apparent);
    /* Light-time applies to everything except the Sun (at the origin). */
    const bool light_time = apparent && body != EPH_SUN;

    /* The compressed table when it covers this date, else the analytical
     * series; bodies the table lacks (asteroids) fall back individually. */
    const bool table = eph_compressed_available() && compressed_covers(jd_tt);
    State e, b;
    if (!cached_state(EPH_EARTH, jd_tt, table, e) &&
        !(table && cached_state(EPH_EARTH, jd_tt, false, e)))
        return false;
    auto body_state = [&](bool t) {
        return body == EPH_MOON ? cached_state(body, jd_tt, t, b)
                                : state(body, jd_tt, t, b);
    };
    if (!body_state(table) && !(table && body_state(false))) return false;

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
    if (body >= EPH_MEAN_NODE && body <= EPH_OSCU_LILITH) { /* directions: no parallax */
        *lon = g.longitude;
        *lat = g.latitude;
        *dist = g.distance;
        return true;
    }
    double ra, dec;
    ecliptic_to_equatorial(g.longitude, g.latitude, jd_tt, &ra, &dec);
    double ra2, dec2;
    topocentric_correct(ra, dec, g.distance, jd_ut, obs, &ra2, &dec2, dist);
    equatorial_to_ecliptic(ra2, dec2, jd_tt, lon, lat);
    return true;
}

}  // namespace

bool earth_heliocentric_state(double jd_tt, Vec3& p, Vec3& v) {
    const bool table = eph_compressed_available() && compressed_covers(jd_tt);
    State e;
    if (!cached_state(EPH_EARTH, jd_tt, table, e) &&
        !(table && cached_state(EPH_EARTH, jd_tt, false, e)))
        return false;
    p = e.p;
    v = e.v;
    return true;
}

void reduce_direction(double jd_tt, bool apparent, const Vec3& g, const Vec3& e,
                      const Vec3& ev, double& lon, double& lat) {
    reduce(EPH_MERCURY, jd_tt, apparent, g, e, ev, lon, lat); /* a generic body */
}

Vec3 ecliptic_j2000_to_date(double jd_tt, const Vec3& v) {
    return j2000_to_date(jd_tt, v);
}

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
    const bool asteroid = (body >= EPH_CERES && body <= EPH_CHIRON) || is_minor(body);
    if (!planet && !asteroid && !is_uranian(body)) return -1;
    Vec3 h;
    if (is_uranian(body) ? !uranian_heliocentric(body, jd_tt, h)
                         : !body_heliocentric(body, jd_tt, h, jd_tt))
        return -1;
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
