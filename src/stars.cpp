#include "internal.h"

#include <cctype>
#include <cstdlib>
#include <cstring>

#include "../include/ephemeris/ephemeris.h"
#include "data.h"

/* Fixed stars from the catalogue in a data file (data/eph-stars.bin, made
 * by tools/make_stars.py and pack_series): space motion from the J2000.0
 * epoch, annual parallax, light deflection by the Sun, aberration with the
 * Earth's barycentric velocity, precession and nutation. */

namespace eph {

namespace {

constexpr double kMasPerYearToRadPerDay = ARCSEC / 1000.0 / 365.25;
constexpr double kAuPerParsec = 206264.80624709636;

const StarRecord* star(int index) {
    const Data& d = data();
    return (index >= 0 && index < d.nstars) ? &d.stars[index] : nullptr;
}

/* The Sun's velocity about the solar-system barycentre (AU/day, mean
 * ecliptic J2000) from Keplerian orbits of the four giant planets
 * (Standish's approximate elements): within 0.3 m/s of DE440 over
 * 1850-2140, i.e. 0.0002" of stellar aberration (the Sun moves up to
 * 16 m/s, 0.01").  For the planets this term cancels against their
 * heliocentric light-time; for the stars it does not. */
Vec3 sun_barycentric_velocity(double jd) {
    struct Orbit {
        double a, e, i, l, w, o;       /* AU, -, deg: I, L, varpi, Omega */
        double da, de, di, dl, dw, dO; /* per Julian century            */
        double mass;                   /* planet / Sun                   */
    };
    static const Orbit kOrbits[4] = {
        {5.20288700, 0.04838624, 1.30439695, 34.39644051, 14.72847983,
         100.47390909, -0.00011607, -0.00013253, -0.00183714, 3034.74612775,
         0.21252668, 0.20469106, 1.0 / 1047.348644},
        {9.53667594, 0.05386179, 2.48599187, 49.95424423, 92.59887831,
         113.66242448, -0.00125060, -0.00050991, 0.00193609, 1222.49362201,
         -0.41897216, -0.28867794, 1.0 / 3497.9018},
        {19.18916464, 0.04725744, 0.77263783, 313.23810451, 170.95427630,
         74.01692503, -0.00196176, -0.00004397, -0.00242939, 428.48202785,
         0.40805281, 0.04240589, 1.0 / 22902.944},
        {30.06992276, 0.00859048, 1.77004347, -55.12002969, 44.96476227,
         131.78422574, 0.00026291, 0.00005105, 0.00035372, 218.45945325,
         -0.32241464, -0.00508664, 1.0 / 19412.26},
    };
    const double k = 0.01720209895; /* Gaussian gravitational constant */
    const double T = (jd - J2000) / DAYS_PER_CENTURY;
    Vec3 v = {0.0, 0.0, 0.0};
    double total = 1.0;
    for (const Orbit& o : kOrbits) {
        const double a = o.a + o.da * T, e = o.e + o.de * T;
        const double inc = (o.i + o.di * T) * DEG;
        const double w = (o.w + o.dw * T) * DEG, node = (o.o + o.dO * T) * DEG;
        const double M = std::fmod((o.l + o.dl * T) * DEG - w, TWO_PI);
        double E = M;
        for (int n = 0; n < 6; ++n) E = M + e * std::sin(E);
        const double edot = k / (a * std::sqrt(a)) / (1.0 - e * std::cos(E));
        const double xd = -a * std::sin(E) * edot;
        const double yd = a * std::sqrt(1.0 - e * e) * std::cos(E) * edot;
        const double cw = std::cos(w - node), sw = std::sin(w - node);
        const double cn = std::cos(node), sn = std::sin(node);
        const double ci = std::cos(inc), si = std::sin(inc);
        v.x -= o.mass * ((cw * cn - sw * sn * ci) * xd - (sw * cn + cw * sn * ci) * yd);
        v.y -= o.mass * ((cw * sn + sw * cn * ci) * xd - (sw * sn - cw * cn * ci) * yd);
        v.z -= o.mass * (sw * si * xd + cw * si * yd);
        total += o.mass;
    }
    return {v.x / total, v.y / total, v.z / total};
}

/* ICRS (equatorial) -> mean ecliptic J2000 with the IAU-1980 obliquity, the
 * frame the series and tables use. */
Vec3 icrs_to_ecliptic(const Vec3& v) {
    const double e = 84381.448 * ARCSEC; /* IAU 1980 obliquity at J2000 */
    const double c = std::cos(e), s = std::sin(e);
    return {v.x, c * v.y + s * v.z, -s * v.y + c * v.z};
}

/* Position at jd_tt, ecliptic of date (IAU 2006); `apparent` adds light
 * deflection, aberration and nutation. */
bool position(int index, double jd_tt, bool apparent, double& lon, double& lat,
              double& dist) {
    const StarRecord* s = star(index);
    if (!s) return false;
    Vec3 earth, earth_v;
    if (!earth_heliocentric_state(jd_tt, earth, earth_v)) return false;

    /* Space motion from J2000.0: straight-line motion of the barycentric
     * position (AU), or of the direction for a star without parallax. */
    const double ca = std::cos(s->ra), sa = std::sin(s->ra);
    const double cd = std::cos(s->dec), sd = std::sin(s->dec);
    const Vec3 u = {cd * ca, cd * sa, sd};
    const Vec3 ea = {-sa, ca, 0.0}, ed = {-sd * ca, -sd * sa, cd};
    /* Straight-line motion of the catalogue (as-observed) quantities, as
     * ERFA's pmsafe/starpm: within 0.03" of it for Barnard's Star over
     * 1600-2500. */
    const bool has_plx = s->parallax > 0.0f;
    const double t = jd_tt - J2000;
    const double ma = s->pm_ra * kMasPerYearToRadPerDay * t;
    const double md = s->pm_dec * kMasPerYearToRadPerDay * t;
    const bool far = !has_plx;
    const double d = far ? 1.0 : kAuPerParsec * 1000.0 / s->parallax;
    const double radial = far ? 0.0 : s->rv * 86400.0 / AU_KM * t;
    Vec3 p = {d * (u.x + ma * ea.x + md * ed.x) + radial * u.x,
              d * (u.y + ma * ea.y + md * ed.y) + radial * u.y,
              d * (u.z + ma * ea.z + md * ed.z) + radial * u.z};
    p = ecliptic_j2000_to_date(jd_tt, icrs_to_ecliptic(p));

    /* Annual parallax: seen from the Earth (heliocentric ~ barycentric
     * here: the 0.01 AU difference is ~1e-8" of parallax). */
    const Vec3 g = far ? p : p - earth;
    dist = far ? 0.0 : length(g);

    /* Aberration with the Earth's barycentric velocity. */
    const Vec3 sv = ecliptic_j2000_to_date(jd_tt, sun_barycentric_velocity(jd_tt));
    const Vec3 ev = {earth_v.x + sv.x, earth_v.y + sv.y, earth_v.z + sv.z};
    reduce_direction(jd_tt, apparent, far ? p : g, earth, ev, lon, lat);
    return true;
}

bool fill(int index, double jd_tt, bool apparent, eph_position* out) {
    double l, b, r, l1, b1, r1, l2, b2, r2;
    const double h = 0.5;
    if (!out || !position(index, jd_tt, apparent, l, b, r) ||
        !position(index, jd_tt - h, apparent, l1, b1, r1) ||
        !position(index, jd_tt + h, apparent, l2, b2, r2))
        return false;
    out->longitude = l;
    out->latitude = b;
    out->distance = r;
    out->longitude_speed = (normalize(l2 - l1 + PI) - PI) / (2 * h);
    out->latitude_speed = (b2 - b1) / (2 * h);
    out->distance_speed = (r2 - r1) / (2 * h);
    return true;
}

/* Case-insensitive comparison, treating runs of spaces as one. */
bool same(const char* a, const char* b) {
    while (*a && *b) {
        if (*a == ' ' && *b == ' ') {
            while (*a == ' ') ++a;
            while (*b == ' ') ++b;
            continue;
        }
        if (std::tolower((unsigned char)*a) != std::tolower((unsigned char)*b))
            return false;
        ++a;
        ++b;
    }
    return *a == *b;
}

/* "alpha Tau" -> "alf Tau" (the catalogue's three-letter Greek forms). */
const char* bayer_abbreviation(const char* word, size_t len) {
    static const char* kGreek[][2] = {
        {"alpha", "alf"}, {"beta", "bet"},   {"gamma", "gam"},   {"delta", "del"},
        {"epsilon", "eps"}, {"zeta", "zet"}, {"eta", "eta"},     {"theta", "tet"},
        {"iota", "iot"},  {"kappa", "kap"},  {"lambda", "lam"},  {"mu", "mu"},
        {"nu", "nu"},     {"xi", "ksi"},     {"omicron", "omi"}, {"pi", "pi"},
        {"rho", "rho"},   {"sigma", "sig"},  {"tau", "tau"},     {"upsilon", "ups"},
        {"phi", "phi"},   {"chi", "chi"},    {"psi", "psi"},     {"omega", "ome"}};
    for (auto& g : kGreek)
        if (std::strlen(g[0]) == len) {
            size_t i = 0;
            while (i < len && std::tolower((unsigned char)word[i]) == g[0][i]) ++i;
            if (i == len) return g[1];
        }
    return nullptr;
}

}  // namespace

}  // namespace eph

using namespace eph;

extern "C" {

int eph_star_count(void) { return data().nstars; }

int eph_star_find(const char* name) {
    const Data& d = data();
    if (!name || !d.stars) return -1;
    while (*name == ' ') ++name;
    /* "HIP 21421" */
    if ((name[0] == 'H' || name[0] == 'h') && (name[1] == 'I' || name[1] == 'i') &&
        (name[2] == 'P' || name[2] == 'p')) {
        const long hip = std::strtol(name + 3, nullptr, 10);
        for (int i = 0; i < d.nstars; ++i)
            if (hip > 0 && d.stars[i].hip == (uint32_t)hip) return i;
        return -1;
    }
    /* "alpha Tau" -> "alf Tau" */
    char alt[64] = "";
    const char* space = std::strchr(name, ' ');
    if (space) {
        const char* ab = bayer_abbreviation(name, space - name);
        if (ab && std::strlen(ab) + std::strlen(space) < sizeof alt) {
            std::strcpy(alt, ab);
            std::strcat(alt, space);
        }
    }
    for (int i = 0; i < d.nstars; ++i) {
        const char* n = d.star_strings + d.stars[i].name;
        const char* des = d.star_strings + d.stars[i].designation;
        if (same(name, n) || (*des && (same(name, des) || (*alt && same(alt, des)))))
            return i;
    }
    return -1;
}

const char* eph_star_name(int index) {
    const StarRecord* s = star(index);
    return s ? data().star_strings + s->name : nullptr;
}

const char* eph_star_designation(int index) {
    const StarRecord* s = star(index);
    return s ? data().star_strings + s->designation : nullptr;
}

double eph_star_magnitude(int index) {
    const StarRecord* s = star(index);
    return s ? s->vmag : std::nan("");
}

int eph_star_position_tt(int index, double jd_tt, eph_position* out) {
    return fill(index, jd_tt, false, out) ? 0 : -1;
}

int eph_star_position_apparent_tt(int index, double jd_tt, eph_position* out) {
    return fill(index, jd_tt, true, out) ? 0 : -1;
}

int eph_star_position_apparent_ut(int index, double jd_ut, double year,
                                  eph_position* out) {
    return eph_star_position_apparent_tt(index, eph_ut_to_tt(jd_ut, year), out);
}

}  // extern "C"
