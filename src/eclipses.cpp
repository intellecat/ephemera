#include "internal.h"

#include "../include/ephemeris/ephemeris.h"

/*
 * Solar and lunar eclipses from the Sun and Moon positions (any data file
 * with the Earth and the Moon), following the Explanatory Supplement to
 * the Astronomical Almanac (1992, ch. 8).
 *
 * Solar eclipses use Besselian elements computed directly from the
 * apparent geocentric positions: the shadow axis runs from the Moon away
 * from the Sun; the fundamental plane passes through the Earth's centre
 * perpendicular to it; (x, y) is where the axis crosses it, l1 and l2 the
 * penumbral and umbral radii there, all in Earth equatorial radii.  The
 * Earth's flattening is handled by scaling y (y1 = y / rho1).
 *
 * Lunar eclipses compare the Moon with the Earth's shadow at the Moon's
 * distance, enlarged for the atmosphere by Danjon's rule (the Earth's
 * radius by 1/85).
 *
 * The searches skip lunations far from a node (Meeus, Astronomical
 * Algorithms, ch. 54), then sample the Sun and Moon once around each
 * candidate and interpolate them (Chebyshev) for the many evaluations the
 * contacts need.
 *
 * All times are TT (Julian days).
 */

namespace {

constexpr double kEarthRadiusKm = 6378.137;
constexpr double kFlattening = 1.0 / 298.257223563;
constexpr double kE2 = kFlattening * (2.0 - kFlattening); /* e^2 */
constexpr double kSunRadiusKm = 696000.0;
/* Moon radius / Earth radius (NASA's eclipse canon): for the penumbra; for the umbra the
 * smaller value through the limb valleys (as NASA's eclipse canon). */
constexpr double kMoonK = 0.272488;
constexpr double kMoonK2 = 0.2722810;
constexpr double kSynodicMonth = 29.530588861;

double year_of(double jd) { return 2000.0 + (jd - eph::J2000) / 365.25; }

/* Apparent geocentric equatorial vector (true equator of date), in Earth
 * equatorial radii. */
bool equatorial(eph_body body, double jd_tt, eph::Vec3& v) {
    eph_position p;
    if (eph_position_apparent_tt(body, jd_tt, &p) != 0) return false;
    double ra, dec;
    eph_ecliptic_to_equatorial(p.longitude, p.latitude, jd_tt, &ra, &dec);
    const double r = p.distance * eph::AU_KM / kEarthRadiusKm;
    v = {r * std::cos(dec) * std::cos(ra), r * std::cos(dec) * std::sin(ra),
         r * std::sin(dec)};
    return true;
}

/* The Sun and Moon vectors over [a, b] as Chebyshev series: degree 8 over
 * a day interpolates the Moon to ~1e-10 Earth radii. */
constexpr int kDegree = 8;
constexpr double kHalfWindow = 0.5;

struct Track {
    double a = 0.0, b = 0.0;
    double c[6][kDegree + 1];
};

bool make_track(double a, double b, Track& tr) {
    const int n = kDegree + 1;
    double f[n][6];
    for (int k = 0; k < n; ++k) {
        const double x = std::cos(eph::PI * (k + 0.5) / n);
        const double t = (a + b) / 2 + x * (b - a) / 2;
        eph::Vec3 s, m;
        if (!equatorial(EPH_SUN, t, s) || !equatorial(EPH_MOON, t, m)) return false;
        const double v[6] = {s.x, s.y, s.z, m.x, m.y, m.z};
        for (int i = 0; i < 6; ++i) f[k][i] = v[i];
    }
    for (int i = 0; i < 6; ++i)
        for (int j = 0; j < n; ++j) {
            double sum = 0.0;
            for (int k = 0; k < n; ++k) sum += f[k][i] * std::cos(eph::PI * j * (k + 0.5) / n);
            tr.c[i][j] = (j == 0 ? 1.0 : 2.0) * sum / n;
        }
    tr.a = a;
    tr.b = b;
    return true;
}

/* The Sun (s) and Moon (m) at t: from the track when it covers t. */
bool vectors(const Track* tr, double t, eph::Vec3& s, eph::Vec3& m) {
    if (!tr || !(t >= tr->a && t <= tr->b))
        return equatorial(EPH_SUN, t, s) && equatorial(EPH_MOON, t, m);
    const double x = (2 * t - tr->a - tr->b) / (tr->b - tr->a);
    double v[6];
    for (int i = 0; i < 6; ++i) { /* Clenshaw */
        double b1 = 0.0, b2 = 0.0;
        for (int j = kDegree; j >= 1; --j) {
            const double b0 = 2 * x * b1 - b2 + tr->c[i][j];
            b2 = b1;
            b1 = b0;
        }
        v[i] = x * b1 - b2 + tr->c[i][0];
    }
    s = {v[0], v[1], v[2]};
    m = {v[3], v[4], v[5]};
    return true;
}

double dot(const eph::Vec3& a, const eph::Vec3& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

/* ---- solar: Besselian elements ------------------------------------------ */

struct Bessel {
    double x, y, z;      /* Moon on the fundamental plane frame       */
    double d, mu;        /* declination, Greenwich hour angle of axis */
    double l1, l2;       /* penumbral, umbral radius on the plane     */
    double tan1, tan2;   /* cone half-angle tangents                  */
    double rho1, d1;     /* flattening: y1 = y / rho1, auxiliary decl. */
};

bool bessel(const Track* tr, double jd_tt, Bessel& b) {
    eph::Vec3 s, m;
    if (!vectors(tr, jd_tt, s, m)) return false;
    const eph::Vec3 g = s - m;
    const double gl = eph::length(g);
    const eph::Vec3 u = {g.x / gl, g.y / gl, g.z / gl};
    b.d = std::asin(u.z);
    const double a = std::atan2(u.y, u.x);
    const eph::Vec3 i = {-std::sin(a), std::cos(a), 0.0};
    const eph::Vec3 j = {-std::sin(b.d) * std::cos(a), -std::sin(b.d) * std::sin(a),
                         std::cos(b.d)};
    b.x = dot(m, i);
    b.y = dot(m, j);
    b.z = dot(m, u);
    const double jd_ut = eph_tt_to_ut(jd_tt, year_of(jd_tt));
    b.mu = eph::normalize(eph_sidereal_time(jd_ut) - a);
    const double rs = kSunRadiusKm / kEarthRadiusKm;
    const double sin1 = (rs + kMoonK) / gl, sin2 = (rs - kMoonK2) / gl;
    b.tan1 = sin1 / std::sqrt(1.0 - sin1 * sin1);
    b.tan2 = sin2 / std::sqrt(1.0 - sin2 * sin2);
    b.l1 = b.z * b.tan1 + kMoonK / std::sqrt(1.0 - sin1 * sin1);
    b.l2 = b.z * b.tan2 - kMoonK2 / std::sqrt(1.0 - sin2 * sin2);
    b.rho1 = std::sqrt(1.0 - kE2 * std::cos(b.d) * std::cos(b.d));
    b.d1 = std::atan2(std::sin(b.d) / b.rho1,
                      std::sqrt(1.0 - kE2) * std::cos(b.d) / b.rho1);
    return true;
}

/* Distance of the shadow axis from the Earth's centre (flattened units). */
double gamma_of(const Bessel& b) { return std::hypot(b.x, b.y / b.rho1); }

/* Signed distance (Earth radii; negative inside) from the shadow axis to
 * the Earth's outline on the fundamental plane, the ellipse x = cos t,
 * y = rho1 sin t; `t` returns the nearest point's parameter. */
double outline_distance(const Bessel& b, double* t_out = nullptr) {
    const double r = b.rho1;
    double t = std::atan2(b.y / r, b.x);
    for (int i = 0; i < 6; ++i) { /* Newton on (P - E(t)) . E'(t) = 0 */
        const double c = std::cos(t), s = std::sin(t);
        const double g = -(b.x - c) * s + (b.y - r * s) * r * c;
        const double dg = -s * s - c * (b.x - c) - r * r * c * c - r * s * (b.y - r * s);
        t -= g / dg;
    }
    if (t_out) *t_out = t;
    const double d = std::hypot(b.x - std::cos(t), b.y - r * std::sin(t));
    return gamma_of(b) < 1.0 ? -d : d;
}

/* The observer in the fundamental frame: xi, eta, zeta (Earth radii). */
void observer_frame(const Bessel& b, const eph_observer& o, double& xi, double& eta,
                    double& zeta) {
    const double phi = o.latitude, h = o.elevation / 1000.0 / kEarthRadiusKm;
    const double u = std::atan(std::sqrt(1.0 - kE2) * std::tan(phi));
    const double rs = std::sqrt(1.0 - kE2) * std::sin(u) + h * std::sin(phi);
    const double rc = std::cos(u) + h * std::cos(phi);
    const double H = b.mu + o.longitude;
    xi = rc * std::sin(H);
    eta = rs * std::cos(b.d) - rc * std::cos(H) * std::sin(b.d);
    zeta = rs * std::sin(b.d) + rc * std::cos(H) * std::cos(b.d);
}

/* ---- root finding and minimisation -------------------------------------- */

/* Zero of f in [a, b] where f(a), f(b) differ in sign (false position,
 * Illinois); f returns false if it cannot be evaluated. */
template <typename F>
bool root(F f, double a, double b, double& t) {
    double fa, fb;
    if (!f(a, fa) || !f(b, fb) || (fa < 0) == (fb < 0)) return false;
    int side = 0;
    for (int i = 0; i < 100 && std::fabs(b - a) > 1e-9; ++i) {
        const double c = (fa * b - fb * a) / (fa - fb);
        double fc;
        if (!f(c, fc)) return false;
        if (fc == 0.0) {
            a = b = c;
            break;
        }
        if ((fc < 0) == (fb < 0)) {
            b = c;
            fb = fc;
            if (side == -1) fa /= 2;
            side = -1;
        } else {
            a = c;
            fa = fc;
            if (side == 1) fb /= 2;
            side = 1;
        }
    }
    t = (a + b) / 2;
    return true;
}

/* Minimum of f in [a, b] (golden section, to ~1e-8 day). */
template <typename F>
bool minimum(F f, double a, double b, double& t, double& value) {
    const double g = (std::sqrt(5.0) - 1.0) / 2.0;
    double c = b - g * (b - a), d = a + g * (b - a), fc, fd;
    if (!f(c, fc) || !f(d, fd)) return false;
    while (b - a > 1e-8) {
        if (fc < fd) {
            b = d;
            d = c;
            fd = fc;
            c = b - g * (b - a);
            if (!f(c, fc)) return false;
        } else {
            a = c;
            c = d;
            fc = fd;
            d = a + g * (b - a);
            if (!f(d, fd)) return false;
        }
    }
    t = (a + b) / 2;
    return f(t, value);
}

/* Meeus (Astronomical Algorithms, ch. 54): for lunation k (k = 0: the
 * new moon of 2000-01-06; full moons at k + 0.5) the approximate time of
 * maximum (minutes) and gamma (~0.01), with u for the shadow radii. */
struct Estimate {
    double jde, gamma, u;
};

Estimate estimate(double k, bool lunar) {
    const double d = eph::DEG, T = k / 1236.85;
    const double E = 1.0 - 0.002516 * T - 0.0000074 * T * T;
    const double M = (2.5534 + 29.10535670 * k) * d;
    const double Mp = (201.5643 + 385.81693528 * k) * d;
    const double F = (160.7108 + 390.67050284 * k) * d;
    const double O = (124.7746 - 1.56375588 * k) * d;
    const double F1 = F - 0.02665 * d * std::sin(O);
    const double A1 = (299.77 + 0.107408 * k) * d;
    double jde = 2451550.09766 + kSynodicMonth * k + 0.00015437 * T * T;
    jde += lunar ? -0.4065 * std::sin(Mp) + 0.1727 * E * std::sin(M)
                 : -0.4075 * std::sin(Mp) + 0.1721 * E * std::sin(M);
    jde += 0.0161 * std::sin(2 * Mp) - 0.0097 * std::sin(2 * F1) +
           0.0073 * E * std::sin(Mp - M) - 0.0050 * E * std::sin(Mp + M) -
           0.0023 * std::sin(Mp - 2 * F1) + 0.0021 * E * std::sin(2 * M) +
           0.0012 * std::sin(Mp + 2 * F1) + 0.0006 * E * std::sin(2 * Mp + M) -
           0.0004 * std::sin(3 * Mp) - 0.0003 * E * std::sin(M + 2 * F1) +
           0.0003 * std::sin(A1) - 0.0002 * E * std::sin(M - 2 * F1) -
           0.0002 * E * std::sin(2 * Mp - M) - 0.0002 * std::sin(O);
    const double P = 0.2070 * E * std::sin(M) + 0.0024 * E * std::sin(2 * M) -
                     0.0392 * std::sin(Mp) + 0.0116 * std::sin(2 * Mp) -
                     0.0073 * E * std::sin(Mp + M) + 0.0067 * E * std::sin(Mp - M) +
                     0.0118 * std::sin(2 * F1);
    const double Q = 5.2207 - 0.0048 * E * std::cos(M) + 0.0020 * E * std::cos(2 * M) -
                     0.3299 * std::cos(Mp) - 0.0060 * E * std::cos(Mp + M) +
                     0.0041 * E * std::cos(Mp - M);
    const double W = std::fabs(std::cos(F1));
    const double gamma = (P * std::cos(F1) + Q * std::sin(F1)) * (1.0 - 0.0048 * W);
    const double u = 0.0059 + 0.0046 * E * std::cos(M) - 0.0182 * std::cos(Mp) +
                     0.0004 * std::cos(2 * Mp) - 0.0005 * std::cos(M + Mp);
    return {jde, gamma, u};
}

/* Whether lunation k may have an eclipse (Meeus' limits with a margin). */
bool possible(const Estimate& e, bool lunar) {
    return std::fabs(e.gamma) < (lunar ? 1.5573 : 1.5433) + e.u + 0.05;
}

/* The lunations to try, from just before (or after) jd. */
double first_lunation(double jd, int direction, bool lunar) {
    const double half = lunar ? 0.5 : 0.0;
    return std::floor((jd - 2451550.09766) / kSynodicMonth - half) + half - direction;
}

/* Saros series: consecutive lunations step 38 series (mod 223). */
int saros(double k, bool lunar) {
    const long n = (long)std::floor(k);
    long s = ((lunar ? 124 : 112) + 38 * n) % 223;
    if (s <= 0) s += 223;
    return (int)s;
}

/* Overlap of two discs (radii r1 = Sun, r2 = Moon, centre distance c) as
 * a fraction of the first disc's area. */
double obscuration(double r1, double r2, double c) {
    if (c >= r1 + r2) return 0.0;
    if (c <= std::fabs(r1 - r2)) return r2 >= r1 ? 1.0 : (r2 * r2) / (r1 * r1);
    const double a1 = std::acos((c * c + r1 * r1 - r2 * r2) / (2 * c * r1));
    const double a2 = std::acos((c * c + r2 * r2 - r1 * r1) / (2 * c * r2));
    const double area = r1 * r1 * (a1 - std::sin(2 * a1) / 2) +
                        r2 * r2 * (a2 - std::sin(2 * a2) / 2);
    return area / (eph::PI * r1 * r1);
}

/* Geodetic latitude / east longitude of the point where the shadow axis
 * meets the Earth at Besselian state b; false if it misses. */
bool axis_point(const Bessel& b, double& lat, double& lon) {
    const double y1 = b.y / b.rho1;
    const double q = 1.0 - b.x * b.x - y1 * y1;
    if (q < 0.0) return false;
    const double z1 = std::sqrt(q);
    const double phi1 = std::asin(y1 * std::cos(b.d1) + z1 * std::sin(b.d1));
    lat = std::atan(std::tan(phi1) / std::sqrt(1.0 - kE2));
    const double H = std::atan2(b.x, z1 * std::cos(b.d1) - y1 * std::sin(b.d1));
    lon = eph::normalize(H - b.mu + eph::PI) - eph::PI;
    return true;
}

/* Local circumstances of a solar eclipse at time t: the distance of the
 * observer from the axis and the shadow radii there (Earth radii), and
 * the Sun's altitude. */
struct LocalState {
    double delta, L1, L2, altitude;
};

/* Fraction of the Sun's diameter covered; inside the umbra or antumbra
 * the ratio of the apparent diameters, Moon / Sun. */
double local_magnitude(const LocalState& s) {
    if (s.delta <= std::fabs(s.L2)) return (s.L1 - s.L2) / (s.L1 + s.L2);
    return std::fmax(0.0, (s.L1 - s.delta) / (s.L1 + s.L2));
}

bool local_state(const Track* tr, double t, const eph_observer& o, LocalState& s) {
    Bessel b;
    if (!bessel(tr, t, b)) return false;
    double xi, eta, zeta;
    observer_frame(b, o, xi, eta, zeta);
    s.delta = std::hypot(b.x - xi, b.y - eta);
    s.L1 = b.l1 - zeta * b.tan1;
    s.L2 = b.l2 - zeta * b.tan2;
    /* Altitude of the Sun (geocentric, near enough for visibility). */
    const double H = b.mu + o.longitude; /* hour angle of the axis ~ Sun */
    s.altitude = std::asin(std::sin(o.latitude) * std::sin(b.d) +
                           std::cos(o.latitude) * std::cos(b.d) * std::cos(H));
    return true;
}

/* ---- lunar ------------------------------------------------------------- */

struct Shadow {
    double sep;        /* Moon centre from the shadow axis (radians) */
    double penumbra;   /* penumbral radius (radians)                 */
    double umbra;      /* umbral radius                              */
    double moon;       /* Moon's angular radius                      */
    double parallax;   /* Moon's equatorial horizontal parallax      */
    double north;      /* sign: Moon north (+1) or south of the axis */
};

bool shadow(const Track* tr, double jd_tt, Shadow& s) {
    eph::Vec3 vs, vm;
    if (!vectors(tr, jd_tt, vs, vm)) return false;
    const double dm = eph::length(vm), ds = eph::length(vs); /* Earth radii */
    const eph::Vec3 um = {vm.x / dm, vm.y / dm, vm.z / dm};
    const eph::Vec3 us = {vs.x / ds, vs.y / ds, vs.z / ds};
    s.sep = std::acos(std::fmax(-1.0, std::fmin(1.0, -dot(um, us))));
    s.parallax = std::asin(1.0 / dm);
    const double pi_s = std::asin(1.0 / ds);
    const double rs = std::asin(kSunRadiusKm / kEarthRadiusKm / ds);
    s.moon = std::asin(kMoonK / dm);
    /* Danjon: the Earth's radius enlarged by 1/85 for the atmosphere;
     * the shadow is also flattened like the Earth (at the Moon's
     * latitude, ~ the mean radius). */
    const double pm = s.parallax * (1.0 + 1.0 / 85.0) * (1.0 - kFlattening / 2.0);
    s.umbra = pm + pi_s - rs;
    s.penumbra = pm + pi_s + rs;
    /* North of the axis: the Moon's ecliptic latitude above the
     * antisolar point's (ecliptic pole of date, to a few arcseconds). */
    const double e = 23.4393 * eph::DEG;
    s.north = (um.y + us.y) * -std::sin(e) + (um.z + us.z) * std::cos(e) >= 0.0 ? 1.0 : -1.0;
    return true;
}

/* The next/previous solar eclipse after jd; `tr` keeps the Sun and Moon
 * around it for the local circumstances. */
bool solar_next(double jd_tt, int direction, eph_solar_eclipse* out, Track& tr) {
    double k = first_lunation(jd_tt, direction, false);
    for (int n = 0; n < 2500; ++n, k += direction) { /* ~200 years */
        const Estimate est = estimate(k, false);
        if (!possible(est, false)) continue;
        const double jde = est.jde;
        if (!make_track(jde - kHalfWindow, jde + kHalfWindow, tr)) return false;
        const Track* t = &tr;
        /* Greatest eclipse: the axis nearest the Earth's centre (as NASA's
         * canon). */
        auto g = [&](double x, double& v) {
            Bessel b;
            if (!bessel(t, x, b)) return false;
            v = std::hypot(b.x, b.y);
            return true;
        };
        double tmax, gmin;
        if (!minimum(g, jde - 0.2, jde + 0.2, tmax, gmin)) return false;
        if ((direction > 0 && tmax <= jd_tt) || (direction < 0 && tmax >= jd_tt))
            continue;
        Bessel b;
        bessel(t, tmax, b);
        const double dmin = outline_distance(b);
        if (dmin > b.l1) continue; /* the penumbra misses the Earth */

        out->maximum = tmax;
        out->gamma = (b.y >= 0 ? 1.0 : -1.0) * std::hypot(b.x, b.y);
        out->saros = saros(k, false);
        /* Contacts: the penumbra (partial), the umbra or antumbra
         * (central) and the axis (center line) touch the Earth. */
        auto contact = [&](int which) {
            return [&, which](double x, double& v) {
                Bessel e;
                if (!bessel(t, x, e)) return false;
                const double r = which == 0 ? e.l1 : which == 1 ? std::fabs(e.l2) : 0.0;
                v = outline_distance(e) - r;
                return true;
            };
        };
        double c;
        out->partial_begin = root(contact(0), tmax - 0.3, tmax, c) ? c : 0.0;
        out->partial_end = root(contact(0), tmax, tmax + 0.3, c) ? c : 0.0;
        const bool central = dmin < 0.0;
        const bool umbral = dmin < std::fabs(b.l2);
        out->central_begin = out->central_end = 0.0;
        out->center_line_begin = out->center_line_end = 0.0;
        if (umbral) {
            if (root(contact(1), tmax - 0.3, tmax, c)) out->central_begin = c;
            if (root(contact(1), tmax, tmax + 0.3, c)) out->central_end = c;
        }
        if (central) {
            if (root(contact(2), tmax - 0.3, tmax, c)) out->center_line_begin = c;
            if (root(contact(2), tmax, tmax + 0.3, c)) out->center_line_end = c;
        }

        /* Where greatest eclipse is seen: the axis point, or for a
         * non-central eclipse the point of the limb nearest the axis. */
        double lat, lon;
        if (!axis_point(b, lat, lon)) {
            double p;
            outline_distance(b, &p);
            Bessel e = b;
            e.x = std::cos(p) * (1.0 - 1e-9);
            e.y = b.rho1 * std::sin(p) * (1.0 - 1e-9);
            axis_point(e, lat, lon);
        }
        const eph_observer place = {lon, lat, 0.0};
        LocalState ls;
        local_state(t, tmax, place, ls);

        /* Type: the umbral radius L2 < 0 there means total; hybrid when its
         * sign differs at the ends of the central line, where zeta = 0. */
        out->type = EPH_ECLIPSE_PARTIAL;
        if (umbral) {
            const double L2 = central ? ls.L2 : b.l2;
            out->type = L2 < 0.0 ? EPH_ECLIPSE_TOTAL : EPH_ECLIPSE_ANNULAR;
            if (central && out->center_line_begin > 0 && out->center_line_end > 0) {
                Bessel e1, e2;
                bessel(t, out->center_line_begin, e1);
                bessel(t, out->center_line_end, e2);
                if ((e1.l2 < 0) != (L2 < 0) || (e2.l2 < 0) != (L2 < 0))
                    out->type = EPH_ECLIPSE_HYBRID;
            }
            if (!central) out->type |= EPH_ECLIPSE_NONCENTRAL;
        }
        out->latitude = lat;
        out->longitude = lon;
        /* NASA's convention: for a non-central eclipse, the fraction of
         * the diameter covered even inside the umbra or antumbra. */
        out->magnitude = central ? local_magnitude(ls)
                                 : (ls.L1 - ls.delta) / (ls.L1 + ls.L2);
        return true;
    }
    return false;
}

void local_eclipse_state(const LocalState& s, eph_local_eclipse_state* out) {
    out->magnitude = local_magnitude(s);
    /* Disc radii in the plane: L1 = Sun + Moon, L2 = Sun - Moon. */
    const double rs = (s.L1 + s.L2) / 2, rm = (s.L1 - s.L2) / 2;
    out->obscuration = obscuration(rs, rm, s.delta);
    out->altitude = s.altitude;
    out->type = s.delta >= s.L1 ? 0
                : s.delta < std::fabs(s.L2)
                    ? (s.L2 < 0 ? EPH_ECLIPSE_TOTAL : EPH_ECLIPSE_ANNULAR)
                    : EPH_ECLIPSE_PARTIAL;
}

}  // namespace

extern "C" {

int eph_solar_eclipse_next(double jd_tt, int direction, eph_solar_eclipse* out) {
    if (!out || (direction != 1 && direction != -1)) return -1;
    Track tr;
    return solar_next(jd_tt, direction, out, tr) ? 0 : -1;
}

int eph_solar_eclipse_where(double jd_tt, double* latitude, double* longitude) {
    Bessel b;
    if (!bessel(nullptr, jd_tt, b)) return -1;
    double lat, lon;
    if (!axis_point(b, lat, lon)) return 1; /* the axis misses the Earth */
    if (latitude) *latitude = lat;
    if (longitude) *longitude = lon;
    return 0;
}

int eph_solar_eclipse_how(double jd_tt, const eph_observer* observer,
                          eph_local_eclipse_state* out) {
    if (!observer || !out) return -1;
    LocalState s;
    if (!local_state(nullptr, jd_tt, *observer, s)) return -1;
    local_eclipse_state(s, out);
    return 0;
}

int eph_solar_eclipse_local(double jd_tt, const eph_observer* observer,
                            int direction, eph_local_solar_eclipse* out) {
    if (!observer || !out || (direction != 1 && direction != -1)) return -1;
    double t = jd_tt;
    Track tr;
    const Track* trp = &tr;
    for (int n = 0; n < 100; ++n) { /* global eclipses: ~40 years */
        eph_solar_eclipse g;
        if (!solar_next(t, direction, &g, tr)) return -1;
        t = g.maximum;
        const eph_observer& o = *observer;
        auto dist = [&](double x, double& v) {
            LocalState s;
            if (!local_state(trp, x, o, s)) return false;
            v = s.delta - s.L1;
            return true;
        };
        double tmax, vmin;
        if (!minimum(dist, g.maximum - 0.2, g.maximum + 0.2, tmax, vmin)) return -1;
        if (vmin >= 0.0) continue; /* not seen from here */
        LocalState s;
        local_state(trp, tmax, o, s);
        out->maximum = tmax;
        double c;
        out->first_contact = root(dist, tmax - 0.25, tmax, c) ? c : 0.0;
        out->fourth_contact = root(dist, tmax, tmax + 0.25, c) ? c : 0.0;
        auto inner = [&](double x, double& v) {
            LocalState z;
            if (!local_state(trp, x, o, z)) return false;
            v = z.delta - std::fabs(z.L2);
            return true;
        };
        out->second_contact = out->third_contact = 0.0;
        if (s.delta < std::fabs(s.L2)) {
            out->second_contact = root(inner, tmax - 0.05, tmax, c) ? c : 0.0;
            out->third_contact = root(inner, tmax, tmax + 0.05, c) ? c : 0.0;
        }
        eph_local_eclipse_state st;
        local_eclipse_state(s, &st);
        out->type = st.type;
        out->magnitude = st.magnitude;
        out->obscuration = st.obscuration;
        out->altitude = st.altitude;
        out->saros = g.saros;
        return 0;
    }
    return -1;
}

int eph_lunar_eclipse_next(double jd_tt, int direction, eph_lunar_eclipse* out) {
    if (!out || (direction != 1 && direction != -1)) return -1;
    double k = first_lunation(jd_tt, direction, true);
    Track tr;
    const Track* t = &tr;
    for (int n = 0; n < 2500; ++n, k += direction) {
        const Estimate est = estimate(k, true);
        if (!possible(est, true)) continue;
        const double jde = est.jde;
        if (!make_track(jde - kHalfWindow, jde + kHalfWindow, tr)) return -1;
        auto sep = [&](double x, double& v) {
            Shadow s;
            if (!shadow(t, x, s)) return false;
            v = s.sep;
            return true;
        };
        double tmax, dmin;
        if (!minimum(sep, jde - 0.2, jde + 0.2, tmax, dmin)) return -1;
        if ((direction > 0 && tmax <= jd_tt) || (direction < 0 && tmax >= jd_tt))
            continue;
        Shadow s;
        shadow(t, tmax, s);
        if (dmin >= s.penumbra + s.moon) continue;

        out->maximum = tmax;
        out->saros = saros(k, true);
        out->gamma = s.north * dmin / s.parallax;
        out->umbral_magnitude = (s.umbra + s.moon - dmin) / (2 * s.moon);
        out->penumbral_magnitude = (s.penumbra + s.moon - dmin) / (2 * s.moon);
        out->type = out->umbral_magnitude >= 1.0   ? EPH_ECLIPSE_TOTAL
                    : out->umbral_magnitude > 0.0 ? EPH_ECLIPSE_PARTIAL
                                                  : EPH_ECLIPSE_PENUMBRAL;
        auto contact = [&](int which) {
            return [&, which](double x, double& v) {
                Shadow z;
                if (!shadow(t, x, z)) return false;
                const double r = which == 0   ? z.penumbra + z.moon
                                 : which == 1 ? z.umbra + z.moon
                                              : z.umbra - z.moon;
                v = z.sep - r;
                return true;
            };
        };
        double c;
        double* begin[3] = {&out->penumbral_begin, &out->partial_begin, &out->total_begin};
        double* end[3] = {&out->penumbral_end, &out->partial_end, &out->total_end};
        for (int i = 0; i < 3; ++i) {
            *begin[i] = root(contact(i), tmax - 0.3, tmax, c) ? c : 0.0;
            *end[i] = root(contact(i), tmax, tmax + 0.3, c) ? c : 0.0;
        }
        return 0;
    }
    return -1;
}

}  // extern "C"
