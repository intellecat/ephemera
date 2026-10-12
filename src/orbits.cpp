#include "internal.h"

#include "../include/ephemeris/ephemeris.h"

/* Planetary nodes and apsides: the orbit's elements, mean or osculating,
 * referred to the mean ecliptic and equinox of date, and the ascending and
 * descending nodes, perihelion and aphelion as points on it. */

namespace {

using eph::Vec3;

/* Mean elements referred to the mean equinox of date, polynomials in
 * Julian centuries from J2000 (Simon et al. 1994, A&A 282, 663, as given by
 * Meeus, Astronomical Algorithms, 2nd ed., Table 31.A): a (AU), e, i, node,
 * longitude of perihelion (degrees). */
struct MeanElements {
    double a[4], e[4], i[4], node[4], peri[4];
};

const MeanElements kMean[] = {
    /* Mercury */
    {{0.38709831, 0, 0, 0},
     {0.20563175, 0.000020407, -0.0000000283, -0.00000000018},
     {7.004986, 0.0018215, -0.0000181, 0.000000056},
     {48.330893, 1.1861883, 0.00017542, 0.000000215},
     {77.456119, 1.5564776, 0.00029544, 0.000000009}},
    /* Venus */
    {{0.72332982, 0, 0, 0},
     {0.00677192, -0.000047765, 0.0000000981, 0.00000000046},
     {3.394662, 0.0010037, -0.00000088, -0.000000007},
     {76.67992, 0.9011206, 0.00040618, -0.000000093},
     {131.563703, 1.4022288, -0.00107618, -0.000005678}},
    /* Earth: no node on its own ecliptic */
    {{0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}},
    /* Mars */
    {{1.523679342, 0, 0, 0},
     {0.09340065, 0.000090484, -0.0000000806, -0.00000000025},
     {1.849726, -0.0006011, 0.00001276, -0.000000007},
     {49.558093, 0.7720959, 0.00001557, 0.000002267},
     {336.060234, 1.8410449, 0.00013477, 0.000000536}},
    /* Jupiter */
    {{5.202603209, 0.0000001913, 0, 0},
     {0.04849793, 0.000163225, -0.0000004714, -0.00000000201},
     {1.303267, -0.0054965, 0.00000466, -0.000000002},
     {100.464407, 1.0209774, 0.00040315, 0.000000404},
     {14.331207, 1.6126352, 0.00103042, -0.000004464}},
    /* Saturn */
    {{9.554909192, -0.000002139, 0.000000004, 0},
     {0.05554814, -0.000346641, -0.0000006436, 0.0000000034},
     {2.488879, -0.0037362, -0.00001519, 0.000000087},
     {113.665503, 0.877088, -0.00012176, -0.000002249},
     {93.057237, 1.9637613, 0.00083753, 0.000004928}},
    /* Uranus */
    {{19.218446062, -0.0000000372, 0.00000000098, 0},
     {0.04638122, -0.000027293, 0.0000000789, 0.00000000024},
     {0.773197, 0.0007744, 0.00003749, -0.000000092},
     {74.005957, 0.5211278, 0.00133947, 0.000018484},
     {173.005291, 1.486379, 0.00021406, 0.000000434}},
    /* Neptune (checked against the VSOP87 secular terms: the copy in
     * pymeeus has the first edition's e, node and perihelion) */
    {{30.110386869, -0.0000001663, 0.00000000069, 0},
     {0.00898809, 0.000006408, -0.0000000008, 0},
     {1.769952, -0.0093082, -0.00000708, 0.000000028},
     {131.784057, 1.1022035, 0.00026006, -0.000000636},
     {48.123691, 1.4262677, 0.00037918, -0.000000003}},
};

double poly(const double c[4], double T) { return c[0] + T * (c[1] + T * (c[2] + T * c[3])); }

/* Sun / planet mass ratios (IAU 2009 / DE440) for GM = k^2 (1 + m). */
double planet_mass(eph_body body) {
    switch (body) {
        case EPH_MERCURY: return 1.0 / 6023657.33;
        case EPH_VENUS: return 1.0 / 408523.719;
        case EPH_MARS: return 1.0 / 3098703.59;
        case EPH_JUPITER: return 1.0 / 1047.348644;
        case EPH_SATURN: return 1.0 / 3497.9018;
        case EPH_URANUS: return 1.0 / 22902.944;
        case EPH_NEPTUNE: return 1.0 / 19412.26;
        case EPH_PLUTO: return 1.0 / 136045556.0;
        default: return 0.0;
    }
}

Vec3 cross(const Vec3& a, const Vec3& b) {
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
double dot(const Vec3& a, const Vec3& b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
Vec3 scale(const Vec3& v, double k) { return {v.x * k, v.y * k, v.z * k}; }
Vec3 add(const Vec3& a, const Vec3& b) { return {a.x + b.x, a.y + b.y, a.z + b.z}; }

struct Elements {
    double a, e, i, node, peri;
};

/* Heliocentric position (AU), mean ecliptic of date: the geometric
 * geocentric position plus the Earth's, so that the DE440 corrections and
 * tables (when loaded) apply as they do to the positions. */
bool heliocentric(eph_body body, double jd, Vec3& v) {
    eph_position b, sun;
    if (eph_position_tt(body, jd, &b) != 0 || eph_position_tt(EPH_SUN, jd, &sun) != 0)
        return false;
    const Vec3 g = eph::spherical_to_cartesian(b.longitude, b.latitude, b.distance);
    const Vec3 s = eph::spherical_to_cartesian(sun.longitude, sun.latitude, sun.distance);
    v = {g.x - s.x, g.y - s.y, g.z - s.z};
    return true;
}

/* Rotation from the mean ecliptic of date t1 to that of t0: through J2000
 * (the rotation J2000 -> date is applied to the basis, then transposed). */
Vec3 date_to_date(double t1, double t0, const Vec3& v) {
    auto to_date = [](double t, const Vec3& x) {
        return eph::frame_1976_to_2006(t, eph::ecliptic_j2000_to_date(t, x));
    };
    const Vec3 ex = to_date(t1, {1, 0, 0}), ey = to_date(t1, {0, 1, 0}),
               ez = to_date(t1, {0, 0, 1});
    const Vec3 j2000 = {dot(ex, v), dot(ey, v), dot(ez, v)}; /* transpose */
    return to_date(t0, j2000);
}

/* Osculating elements from the heliocentric state: the velocity from a
 * fourth-order central difference of positions carried to one frame. */
bool osculating(eph_body body, double jd, Elements& el) {
    if (body == EPH_EARTH || body == EPH_SUN || body == EPH_MOON ||
        (body >= EPH_MEAN_NODE && !eph::is_minor(body)))
        return false;
    const double h = 0.1;
    Vec3 p0;
    if (!heliocentric(body, jd, p0)) return false;
    Vec3 p = p0, q[4];
    const double dt[4] = {-2 * h, -h, h, 2 * h};
    for (int k = 0; k < 4; ++k) {
        Vec3 x;
        if (!heliocentric(body, jd + dt[k], x)) return false;
        q[k] = date_to_date(jd + dt[k], jd, x);
    }
    const Vec3 v = scale(add(add(scale(q[0], 1.0), scale(q[1], -8.0)),
                             add(scale(q[2], 8.0), scale(q[3], -1.0))),
                         1.0 / (12 * h));
    const double k = 0.01720209895;
    const double mu = k * k * (1.0 + planet_mass(body));
    const double r = eph::length(p);
    const Vec3 hv = cross(p, v);
    const double hn = eph::length(hv);
    const Vec3 ev = add(scale(cross(v, hv), 1.0 / mu), scale(p, -1.0 / r));
    el.e = eph::length(ev);
    el.a = 1.0 / (2.0 / r - dot(v, v) / mu);
    el.i = std::acos(hv.z / hn);
    el.node = eph::normalize(std::atan2(hv.x, -hv.y));
    const Vec3 n = {std::cos(el.node), std::sin(el.node), 0.0};
    const Vec3 hu = scale(hv, 1.0 / hn);
    const double omega = std::atan2(dot(cross(n, ev), hu), dot(n, ev));
    el.peri = eph::normalize(el.node + omega);
    return true;
}

bool mean(eph_body body, double jd, Elements& el) {
    if (body < EPH_MERCURY || body > EPH_NEPTUNE || body == EPH_EARTH) return false;
    const MeanElements& m = kMean[body];
    const double T = (jd - eph::J2000) / eph::DAYS_PER_CENTURY;
    el.a = poly(m.a, T);
    el.e = poly(m.e, T);
    el.i = poly(m.i, T) * eph::DEG;
    el.node = eph::normalize(poly(m.node, T) * eph::DEG);
    el.peri = eph::normalize(poly(m.peri, T) * eph::DEG);
    return true;
}

/* The four points (heliocentric, AU): at argument of latitude 0, pi,
 * omega and omega + pi, at the orbit's distance there. */
void points(const Elements& el, Vec3 out[4]) {
    const Vec3 n = {std::cos(el.node), std::sin(el.node), 0.0};
    const Vec3 hu = {std::sin(el.i) * std::sin(el.node), -std::sin(el.i) * std::cos(el.node),
                     std::cos(el.i)};
    const Vec3 q = cross(hu, n);
    const double omega = el.peri - el.node, p = el.a * (1.0 - el.e * el.e);
    const double u[4] = {0.0, eph::PI, omega, omega + eph::PI};
    for (int k = 0; k < 4; ++k) {
        double r = p / (1.0 + el.e * std::cos(u[k] - omega));
        if (k == EPH_APHELION && el.e >= 1.0) r = 0.0; /* unbound: no aphelion */
        out[k] = scale(add(scale(n, std::cos(u[k])), scale(q, std::sin(u[k]))), r);
    }
}

/* Elements and points at jd; the geocentric points need the Sun. */
bool evaluate(eph_body body, double jd, eph_orbit_method method, Elements& el,
              Vec3 helio[4], Vec3 geo[4], double& dpsi) {
    const bool ok = method == EPH_ORBIT_MEAN         ? mean(body, jd, el)
                    : method == EPH_ORBIT_OSCULATING ? osculating(body, jd, el)
                                                     : false;
    if (!ok) return false;
    points(el, helio);
    eph_position sun;
    if (eph_position_tt(EPH_SUN, jd, &sun) != 0) return false;
    const Vec3 earth =
        scale(eph::spherical_to_cartesian(sun.longitude, sun.latitude, sun.distance), -1.0);
    for (int k = 0; k < 4; ++k) geo[k] = add(helio[k], scale(earth, -1.0));
    double deps;
    eph::nutation(jd, dpsi, deps);
    return true;
}

void spherical(const Vec3& v, double dpsi, double& lon, double& lat, double& r) {
    eph::cartesian_to_spherical(v, lon, lat, r);
    lon = eph::normalize(lon + dpsi);
}

}  // namespace

extern "C" {

int eph_orbit_tt(eph_body body, double jd_tt, eph_orbit_method method, eph_orbit* out) {
    if (!out) return -1;
    Elements el, e1, e2;
    Vec3 h[4], g[4], h1[4], g1[4], h2[4], g2[4];
    double dpsi, dpsi1, dpsi2;
    const double dt = 0.5; /* speeds: central difference over a day */
    if (!evaluate(body, jd_tt, method, el, h, g, dpsi) ||
        !evaluate(body, jd_tt - dt, method, e1, h1, g1, dpsi1) ||
        !evaluate(body, jd_tt + dt, method, e2, h2, g2, dpsi2))
        return -1;
    out->semi_major_axis = el.a;
    out->eccentricity = el.e;
    out->inclination = el.i;
    out->node = el.node;
    out->perihelion = el.peri;
    for (int k = 0; k < 4; ++k) {
        for (int geo = 0; geo < 2; ++geo) {
            eph_position& p = geo ? out->geocentric[k] : out->heliocentric[k];
            const double n0 = geo ? dpsi : 0.0, n1 = geo ? dpsi1 : 0.0, n2 = geo ? dpsi2 : 0.0;
            double l1, b1, r1, l2, b2, r2;
            spherical(geo ? g[k] : h[k], n0, p.longitude, p.latitude, p.distance);
            spherical(geo ? g1[k] : h1[k], n1, l1, b1, r1);
            spherical(geo ? g2[k] : h2[k], n2, l2, b2, r2);
            p.longitude_speed = (eph::normalize(l2 - l1 + eph::PI) - eph::PI) / (2 * dt);
            p.latitude_speed = (b2 - b1) / (2 * dt);
            p.distance_speed = (r2 - r1) / (2 * dt);
        }
    }
    return 0;
}

int eph_orbit_ut(eph_body body, double jd_ut, double year, eph_orbit_method method,
                 eph_orbit* out) {
    return eph_orbit_tt(body, eph_ut_to_tt(jd_ut, year), method, out);
}

}  // extern "C"
