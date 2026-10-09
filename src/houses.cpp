#include "internal.h"

#include "../include/ephemeris/ephemeris.h"

namespace {

double true_obliquity(double jd_tt) {
    double dpsi, deps;
    eph::nutation(jd_tt, dpsi, deps);
    return eph::mean_obliquity(jd_tt) + deps;
}

/* Ecliptic longitude of the point with right ascension ra. */
double ra_to_ecliptic(double ra, double eps) {
    return eph::normalize(std::atan2(std::sin(ra),
                                     std::cos(ra) * std::cos(eps)));
}

/* Declination of the ecliptic point at longitude lon. */
double ecliptic_declination(double lon, double eps) {
    return std::asin(std::sin(eps) * std::sin(lon));
}

/* The ecliptic point on the great circle through the north and south
 * points of the horizon that meets the equator at right ascension ra and
 * has pole height pole (the "oblique ascension" construction behind the
 * Ascendant, Regiomontanus, Campanus and Topocentric cusps).  Of the two
 * (opposite) intersections, the one in the eastern half of that circle:
 * right ascension between ramc and ramc + 180 degrees. */
double oblique_cusp(double ra, double pole, double eps, double ramc) {
    double lon = std::atan2(std::sin(ra), std::cos(ra) * std::cos(eps) -
                                              std::sin(eps) * std::tan(pole));
    const double alpha = std::atan2(std::sin(lon) * std::cos(eps), std::cos(lon));
    if (std::sin(alpha - ramc) < 0.0) lon += eph::PI;
    return eph::normalize(lon);
}

/* The Ascendant when the sidereal time (RAMC) is ramc. */
double ascendant(double ramc, double phi, double eps) {
    return oblique_cusp(ramc + eph::PI / 2, phi, eps, ramc);
}

/* True if a..g run around the ecliptic in one direction, each step less
 * than 180 degrees (sign +1 forwards, -1 backwards). */
bool ordered(const double* c, int n, int sign) {
    for (int i = 0; i + 1 < n; ++i) {
        const double step = eph::normalize(sign * (c[i + 1] - c[i]));
        if (step >= eph::PI) return false;
    }
    return true;
}

/* Diurnal semi-arc (radians of right ascension) of a point of declination
 * dec, or a negative value if it never rises or never sets. */
double semi_arc(double dec, double phi) {
    const double x = -std::tan(phi) * std::tan(dec);
    return (x <= -1.0 || x >= 1.0) ? -1.0 : std::acos(x);
}

/* Fill cusps 4..9 from 10..3 (each opposite) and set 1 and 10. */
void complete(double asc, double mc, double c11, double c12, double c2,
              double c3, double cusps[12]) {
    const double c[6] = {asc, c2, c3, mc, c11, c12}; /* houses 1,2,3,10,11,12 */
    const int idx[6] = {0, 1, 2, 9, 10, 11};
    for (int i = 0; i < 6; ++i) {
        cusps[idx[i]] = eph::normalize(c[i]);
        cusps[(idx[i] + 6) % 12] = eph::normalize(c[i] + eph::PI);
    }
}

/* Placidus: each intermediate cusp is the ecliptic point whose hour angle
 * is the given fraction of its own semi-arc: diurnal, from the meridian
 * (11: 1/3, 12: 2/3), or nocturnal, from the horizon (2: 1/3, 3: 2/3);
 * solved by iteration.  False where a cusp is circumpolar. */
bool placidus_cusp(double ramc, double phi, double eps, double frac,
                   bool above, double& lon) {
    lon = ra_to_ecliptic(ramc + (above ? frac : 1.0 + frac) * eph::PI / 2, eps);
    for (int i = 0; i < 100; ++i) {
        const double dsa = semi_arc(ecliptic_declination(lon, eps), phi);
        if (dsa < 0.0) return false;
        const double ra = above ? ramc + frac * dsa
                                : ramc + eph::PI - (1.0 - frac) * (eph::PI - dsa);
        const double next = ra_to_ecliptic(ra, eps);
        const double step = eph::normalize(next - lon + eph::PI) - eph::PI;
        lon = next;
        if (std::fabs(step) < 1e-12) return true;
    }
    return false;
}

/* Inside the polar circles the cusps of the great-circle systems can run
 * backwards around the ecliptic; the MC is then the opposite meridian
 * point, so that houses 10, 11, 12, 1, 2, 3, 4 stay in order. */
double consistent_mc(double mc, double asc, const double c[4]) {
    double seq[7] = {mc, c[0], c[1], asc, c[2], c[3], mc + eph::PI};
    if (ordered(seq, 7, +1)) return mc;
    seq[0] = mc + eph::PI;
    seq[6] = mc;
    return ordered(seq, 7, -1) ? eph::normalize(mc + eph::PI) : mc;
}

}  // namespace

extern "C" {

int eph_asc_mc(double jd_ut, const eph_observer* obs, double* asc,
               double* mc) {
    if (!obs) return -1;
    const double eps = true_obliquity(jd_ut);
    const double ramc = eph::normalize(eph::gast_degrees(jd_ut) * eph::DEG +
                                       obs->longitude);
    const double ce = std::cos(eps), se = std::sin(eps);
    const double phi = obs->latitude;

    const double lam_mc = ra_to_ecliptic(ramc, eps);
    double lam_asc = eph::normalize(
        std::atan2(std::cos(ramc), -(std::sin(ramc) * ce +
                                     std::tan(phi) * se)));
    /* The ascendant must be east of the MC by 0..180 degrees. */
    if (eph::normalize(lam_asc - lam_mc) > eph::PI)
        lam_asc = eph::normalize(lam_asc + eph::PI);

    if (asc) *asc = lam_asc;
    if (mc) *mc = lam_mc;
    return 0;
}

int eph_houses(double jd_ut, const eph_observer* obs, eph_house_system sys,
               double cusps[12]) {
    if (!obs || !cusps) return -1;
    double asc, mc;
    if (eph_asc_mc(jd_ut, obs, &asc, &mc) != 0) return -1;
    const double eps = true_obliquity(jd_ut);
    const double ramc = eph::normalize(eph::gast_degrees(jd_ut) * eph::DEG +
                                       obs->longitude);
    const double phi = obs->latitude;
    const double third = eph::PI / 6; /* 30 degrees */

    /* Placidus and Koch are undefined inside the polar circles. */
    const bool polar = std::fabs(phi) >= eph::PI / 2 - eps;
    switch (sys) {
        case EPH_HOUSES_PLACIDUS: {
            double c11, c12, c2, c3;
            if (polar || !placidus_cusp(ramc, phi, eps, 1.0 / 3, true, c11) ||
                !placidus_cusp(ramc, phi, eps, 2.0 / 3, true, c12) ||
                !placidus_cusp(ramc, phi, eps, 1.0 / 3, false, c2) ||
                !placidus_cusp(ramc, phi, eps, 2.0 / 3, false, c3)) {
                eph_houses(jd_ut, obs, EPH_HOUSES_PORPHYRY, cusps);
                return 1;
            }
            complete(asc, mc, c11, c12, c2, c3, cusps);
            break;
        }
        case EPH_HOUSES_KOCH: {
            /* Trisect the time the MC degree takes to rise to the
             * meridian: cusps 11, 12, 2, 3 are the Ascendants at those
             * times (and the MC is the Ascendant when it rose). */
            const double dsa = semi_arc(ecliptic_declination(mc, eps), phi);
            if (polar || dsa < 0.0) {
                eph_houses(jd_ut, obs, EPH_HOUSES_PORPHYRY, cusps);
                return 1;
            }
            complete(asc, mc, ascendant(ramc - 2 * dsa / 3, phi, eps),
                     ascendant(ramc - dsa / 3, phi, eps),
                     ascendant(ramc + dsa / 3, phi, eps),
                     ascendant(ramc + 2 * dsa / 3, phi, eps), cusps);
            break;
        }
        case EPH_HOUSES_REGIOMONTANUS:
        case EPH_HOUSES_CAMPANUS: {
            /* Regiomontanus divides the equator, Campanus the prime
             * vertical, in 30 degree steps from the meridian; both cusps
             * lie on great circles through the north and south points. */
            double c[4];
            const int step[4] = {1, 2, 4, 5}; /* houses 11, 12, 2, 3 */
            for (int i = 0; i < 4; ++i) {
                double h = step[i] * third; /* from the meridian, eastwards */
                if (sys == EPH_HOUSES_CAMPANUS)
                    h = std::atan2(std::sin(h) * std::cos(phi), std::cos(h));
                const double pole = std::atan(std::tan(phi) * std::sin(h));
                c[i] = oblique_cusp(ramc + h, pole, eps, ramc);
            }
            complete(asc, consistent_mc(mc, asc, c), c[0], c[1], c[2], c[3], cusps);
            break;
        }
        case EPH_HOUSES_ALCABITIUS: {
            /* Trisect the Ascendant's diurnal and nocturnal semi-arcs in
             * right ascension. */
            const double dsa = semi_arc(ecliptic_declination(asc, eps), phi);
            if (dsa < 0.0) {
                eph_houses(jd_ut, obs, EPH_HOUSES_PORPHYRY, cusps);
                return 1;
            }
            const double nsa = eph::PI - dsa;
            complete(asc, mc, ra_to_ecliptic(ramc + dsa / 3, eps),
                     ra_to_ecliptic(ramc + 2 * dsa / 3, eps),
                     ra_to_ecliptic(ramc + eph::PI - 2 * nsa / 3, eps),
                     ra_to_ecliptic(ramc + eph::PI - nsa / 3, eps), cusps);
            break;
        }
        case EPH_HOUSES_TOPOCENTRIC: {
            /* Polich-Page: 30 degree steps of the equator with pole heights
             * tan P = tan(phi) * k / 3 (k = 1, 2 for houses 11, 12 and 3,
             * 2).  Inside the polar circles, where the system has no
             * accepted form, the cusps follow the same rule as
             * Regiomontanus (other software may differ there). */
            const double t = std::tan(phi);
            const double c[4] = {
                oblique_cusp(ramc + third, std::atan(t / 3), eps, ramc),
                oblique_cusp(ramc + 2 * third, std::atan(2 * t / 3), eps, ramc),
                oblique_cusp(ramc + 4 * third, std::atan(2 * t / 3), eps, ramc),
                oblique_cusp(ramc + 5 * third, std::atan(t / 3), eps, ramc)};
            complete(asc, consistent_mc(mc, asc, c), c[0], c[1], c[2], c[3], cusps);
            break;
        }
        case EPH_HOUSES_VEHLOW:
            for (int i = 0; i < 12; ++i)
                cusps[i] = eph::normalize(asc + (i - 0.5) * third);
            break;
        case EPH_HOUSES_EQUAL_MC:
            for (int i = 0; i < 12; ++i)
                cusps[i] = eph::normalize(mc + (i + 3) * third);
            break;
        case EPH_HOUSES_EQUAL_ARIES:
            for (int i = 0; i < 12; ++i) cusps[i] = i * third;
            break;
        case EPH_HOUSES_MERIDIAN:
            /* Equal 30 degree steps of right ascension from the RAMC. */
            for (int i = 0; i < 12; ++i)
                cusps[i] = ra_to_ecliptic(ramc + (i + 3) * third, eps);
            break;
        case EPH_HOUSES_MORINUS:
            /* The same equator points, projected onto the ecliptic along
             * circles through the ecliptic poles. */
            for (int i = 0; i < 12; ++i) {
                const double ra = ramc + (i + 3) * third;
                cusps[i] = eph::normalize(
                    std::atan2(std::sin(ra) * std::cos(eps), std::cos(ra)));
            }
            break;
        case EPH_HOUSES_WHOLE_SIGN: {
            const double start = std::floor(asc / (30.0 * eph::DEG)) *
                                 30.0 * eph::DEG;
            for (int i = 0; i < 12; ++i)
                cusps[i] = eph::normalize(start + i * 30.0 * eph::DEG);
            break;
        }
        case EPH_HOUSES_EQUAL: {
            for (int i = 0; i < 12; ++i)
                cusps[i] = eph::normalize(asc + i * 30.0 * eph::DEG);
            break;
        }
        case EPH_HOUSES_PORPHYRY: {
            const double ic = eph::normalize(mc + eph::PI);
            const double d1 = eph::normalize(asc - mc);   /* MC -> Asc  */
            const double d2 = eph::normalize(ic - asc);   /* Asc -> IC  */
            cusps[9] = mc;                                /* 10 = MC    */
            cusps[0] = asc;                               /* 1  = Asc   */
            cusps[3] = ic;                                /* 4  = IC    */
            cusps[6] = eph::normalize(asc + eph::PI);     /* 7  = Desc  */
            for (int i = 1; i <= 2; ++i) {
                cusps[9 + i] = eph::normalize(mc + d1 * i / 3.0);   /* 11,12 */
                cusps[i] = eph::normalize(asc + d2 * i / 3.0);      /* 2,3   */
                cusps[(9 + i + 6) % 12] = eph::normalize(cusps[9 + i] + eph::PI);
                cusps[(i + 6) % 12] = eph::normalize(cusps[i] + eph::PI);
            }
            break;
        }
        default:
            return -1;
    }
    return 0;
}

int eph_house_angles(double jd_ut, const eph_observer* obs, eph_angles* out) {
    if (!obs || !out) return -1;
    double asc, mc;
    if (eph_asc_mc(jd_ut, obs, &asc, &mc) != 0) return -1;
    const double eps = true_obliquity(jd_ut);
    const double ramc = eph::normalize(eph::gast_degrees(jd_ut) * eph::DEG +
                                       obs->longitude);
    const double phi = obs->latitude;
    const double ce = std::cos(eps), se = std::sin(eps);
    out->ascendant = asc;
    out->mc = mc;
    out->armc = ramc;

    /* Vertex: the ecliptic meets the prime vertical (whose plane is normal
     * to the north point of the horizon) in the west. */
    const double cd = std::cos(eph::PI / 2 - phi), sd = std::sin(eph::PI / 2 - phi);
    const double n[3] = {cd * std::cos(ramc + eph::PI), cd * std::sin(ramc + eph::PI), sd};
    const double pole[3] = {0.0, -se, ce}; /* ecliptic pole, equatorial */
    double v[3] = {pole[1] * n[2] - pole[2] * n[1], pole[2] * n[0] - pole[0] * n[2],
                   pole[0] * n[1] - pole[1] * n[0]};
    if (v[0] * std::cos(ramc - eph::PI / 2) + v[1] * std::sin(ramc - eph::PI / 2) < 0.0)
        for (double& x : v) x = -x;
    out->vertex = eph::normalize(std::atan2(v[1] * ce + v[2] * se, v[0]));

    out->equatorial_ascendant = ra_to_ecliptic(ramc + eph::PI / 2, eps);
    /* Koch's co-ascendant is the Ascendant for the opposite latitude and
     * the polar ascendant its opposite point; Munkasey's co-ascendant uses
     * the co-latitude (the raw oblique-ascension solution). */
    out->co_ascendant_koch = ascendant(ramc, -phi, eps);
    out->polar_ascendant = eph::normalize(out->co_ascendant_koch + eph::PI);
    const double ra = ramc + eph::PI / 2;
    out->co_ascendant_munkasey = eph::normalize(std::atan2(
        std::sin(ra), std::cos(ra) * ce - se * std::tan(eph::PI / 2 - phi)));
    return 0;
}

double eph_house_position(const double cusps[12], double longitude) {
    for (int i = 0; i < 12; ++i) {
        const double width = eph::normalize(cusps[(i + 1) % 12] - cusps[i]);
        const double into = eph::normalize(longitude - cusps[i]);
        if (into < width) return i + 1 + into / width;
    }
    return 1.0; /* not reached for valid cusps */
}

}  // extern "C"
