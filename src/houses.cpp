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

}  // namespace

extern "C" {

double eph_ayanamsa_lahiri(double jd_tt) {
    /* General precession in longitude since J2000 (IAU 2006, matching the
     * output frame) plus the Lahiri ayanamsa at J2000 (23 deg 51' 11.6").
     * Ayanamsa = tropical - sidereal. */
    const double T = (jd_tt - eph::J2000) / eph::DAYS_PER_CENTURY;
    const double p =
        T * (5028.796195 +
             T * (1.1054348 +
                  T * (0.00007964 + T * (-0.000023857 - T * 0.0000000383)))) *
        eph::ARCSEC;
    return 23.853222 * eph::DEG + p;
}

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

    switch (sys) {
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

}  // extern "C"
