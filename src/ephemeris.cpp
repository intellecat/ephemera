#include "internal.h"

#include "../include/ephemeris/ephemeris.h"

namespace eph {
namespace {

Vec3 earth_heliocentric(double jd) {
    double l, b, r;
    vsop87_heliocentric(EPH_EARTH, jd, l, b, r);
    return spherical_to_cartesian(l, b, r);
}

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

Vec3 body_heliocentric(eph_body body, double jd) {
    if (body == EPH_PLUTO) return pluto_heliocentric(jd);
    return planet_heliocentric(static_cast<int>(body), jd);
}

Vec3 normalize_vec(const Vec3& v) {
    const double n = length(v);
    return (n > 0.0) ? Vec3{v.x / n, v.y / n, v.z / n} : v;
}

/* Geometric geocentric vector; optionally retarded for light-time. */
bool geocentric_vector(eph_body body, double jd, bool light_time, Vec3& out,
                       double& dist) {
    if (body == EPH_EARTH) return false;

    /* Optional high-accuracy Chebyshev backend (ecliptic J2000). */
    if (eph_compressed_available()) {
        Vec3 cv;
        double cd;
        if (compressed_geocentric((int)body, jd, light_time, cv, cd)) {
            double l = std::atan2(cv.y, cv.x);
            double b = std::asin(cv.z / cd);
            precess_ecliptic_j2000_to_date(jd, l, b);
            out = spherical_to_cartesian(l, b, cd);
            dist = cd;
            return true;
        }
    }

    const Vec3 earth = earth_heliocentric(jd);

    if (body == EPH_MOON) {
        double l, b, r;
        moon_position(jd, l, b, r);
        out = spherical_to_cartesian(l, b, r);
        dist = r;
        if (light_time) {
            for (int i = 0; i < 2; ++i) {
                const double t_emit = jd - dist / AU_PER_LIGHT_DAY;
                moon_position(t_emit, l, b, r);
                out = spherical_to_cartesian(l, b, r);
                dist = r;
            }
        }
        return true;
    }
    if (body == EPH_SUN) {
        out = Vec3{-earth.x, -earth.y, -earth.z};
        dist = length(out);
        return true;
    }

    Vec3 p = body_heliocentric(body, jd);
    out = p - earth;
    dist = length(out);

    if (light_time) {
        for (int i = 0; i < 3; ++i) {
            const double t_emit = jd - dist / AU_PER_LIGHT_DAY;
            p = body_heliocentric(body, t_emit);
            out = p - earth;
            dist = length(out);
        }
    }
    return true;
}

Vec3 earth_velocity(double jd) {
    const double h = 0.01;
    const Vec3 a = earth_heliocentric(jd + h);
    const Vec3 b = earth_heliocentric(jd - h);
    return {(a.x - b.x) / (2 * h), (a.y - b.y) / (2 * h), (a.z - b.z) / (2 * h)};
}

double angle_diff(double a, double b) { return normalize(a - b + PI) - PI; }

double speed_dt(eph_body body) {
    return (body == EPH_MOON) ? 0.001 : 0.1;
}

enum class Kind { Geometric, Apparent };

bool compute(eph_body body, double jd_tt, Kind kind, eph_position* out) {
    if (body < 0 || body >= EPH_BODY_COUNT) return false;

    const bool apparent = (kind == Kind::Apparent);
    /* Light-time applies to everything except the Sun (at the origin).  The
     * Moon gets light-time but no annual aberration: as a geocentric body
     * its apparent place is dominated by the ~1.3 s retarded position. */
    const bool light_time = apparent && body != EPH_SUN;

    Vec3 v;
    double dist;
    if (!geocentric_vector(body, jd_tt, light_time, v, dist)) return false;

    if (apparent) {
        if (body != EPH_MOON) {
            const Vec3 vel = earth_velocity(jd_tt);
            const Vec3 u = normalize_vec(v);
            v = normalize_vec(
                Vec3{u.x + vel.x / AU_PER_LIGHT_DAY,
                     u.y + vel.y / AU_PER_LIGHT_DAY,
                     u.z + vel.z / AU_PER_LIGHT_DAY});
            v = Vec3{v.x * dist, v.y * dist, v.z * dist};
        }
    }

    double lon, lat;
    cartesian_to_spherical(v, lon, lat, dist);
    if (apparent) {
        double dpsi, deps;
        nutation(jd_tt, dpsi, deps);
        lon = normalize(lon + dpsi);
    }

    out->longitude = lon;
    out->latitude = lat;
    out->distance = dist;

    /* Speeds by central difference of the same quantity. */
    const double dt = speed_dt(body);
    eph_position a, b;
    auto eval = [&](double j, eph_position& p) -> bool {
        Vec3 vv;
        double dd;
        if (!geocentric_vector(body, j, light_time, vv, dd)) return false;
        if (apparent && body != EPH_MOON) {
            const Vec3 vel = earth_velocity(j);
            const Vec3 u = normalize_vec(vv);
            vv = normalize_vec(
                Vec3{u.x + vel.x / AU_PER_LIGHT_DAY,
                     u.y + vel.y / AU_PER_LIGHT_DAY,
                     u.z + vel.z / AU_PER_LIGHT_DAY});
        }
        double lo, la;
        cartesian_to_spherical(vv, lo, la, dd);
        if (apparent) {
            double dp, de;
            nutation(j, dp, de);
            lo = normalize(lo + dp);
        }
        p.longitude = lo;
        p.latitude = la;
        p.distance = dd;
        return true;
    };
    if (!eval(jd_tt - dt, a) || !eval(jd_tt + dt, b)) return false;

    out->longitude_speed = angle_diff(b.longitude, a.longitude) / (2 * dt);
    out->latitude_speed = (b.latitude - a.latitude) / (2 * dt);
    out->distance_speed = (b.distance - a.distance) / (2 * dt);
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
    if (body < EPH_MERCURY || body > EPH_PLUTO) return -1;
    Vec3 v = body_heliocentric(body, jd_tt);
    double r;
    double l, b;
    cartesian_to_spherical(v, l, b, r);
    if (lon) *lon = l;
    if (lat) *lat = b;
    if (distance) *distance = r;
    return 0;
}

}  // extern "C"

}  // namespace eph
