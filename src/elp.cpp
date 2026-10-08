#include "internal.h"

#include "../data/elp_data.h"

namespace eph {

/*
 * ELP2000-82B (Chapront-Touze & Chapront) geocentric lunar position.
 *
 * Evaluates the folded Poisson series and performs the ELP coordinate
 * change to rectangular coordinates in the mean dynamical ecliptic and
 * inertial equinox of J2000.  Output in kilometres.
 */
void elp_moon_j2000(double jd_tt, double& x, double& y, double& z) {
    const double t = (jd_tt - J2000) / DAYS_PER_CENTURY;
    const double t2 = t * t, t3 = t2 * t, t4 = t3 * t;

    const double pscale = TWO_PI / 1073741824.0; /* 30-bit phase -> radians */
    const unsigned int pmask = 0x3FFFFFFFu;

    double r[3] = {0.0, 0.0, 0.0};
    for (int iv = 0; iv < 3; ++iv) {
        const elp::Coord& c = elp::coords[iv];
        double sum = 0.0;
        for (int n = 0; n < c.nmain; ++n) {
            const elp::MainTerm& term = c.main[n];
            const double ph0 = (term.packed & pmask) * pscale;
            const double ph = ph0 + term.ph1 * t + term.ph2 * t2 +
                              term.ph3 * t3 + term.ph4 * t4;
            sum += term.amp * std::sin(ph);
        }
        for (int n = 0; n < c.npert; ++n) {
            const elp::PertTerm& term = c.pert[n];
            const double ph0 = (term.packed & pmask) * pscale;
            const int tpow = term.packed >> 30;
            double a = term.amp;
            if (tpow == 1)
                a *= t;
            else if (tpow == 2)
                a *= t2;
            sum += a * std::sin(ph0 + term.ph1 * t);
        }
        r[iv] = sum;
    }

    const double rad = 648000.0 / PI; /* arcseconds per radian */
    const double a0 = 384747.9806448954;
    const double ath = 384747.9806743165;

    const double w10 = (218.0 + 18.0 / 60.0 + 59.95571 / 3600.0) * DEG;
    const double w11 = 1732559343.73604 / rad;
    const double w12 = -5.8883 / rad;
    const double w13 = 0.6604e-2 / rad;
    const double w14 = -0.3169e-4 / rad;

    const double lon = r[0] / rad + w10 + w11 * t + w12 * t2 + w13 * t3 +
                       w14 * t4;
    const double lat = r[1] / rad;
    const double dist = r[2] * a0 / ath;

    double x1 = dist * std::cos(lat);
    double x2 = x1 * std::sin(lon);
    x1 = x1 * std::cos(lon);
    const double x3 = dist * std::sin(lat);

    /* ELP2000-85 -> J2000 inertial ecliptic rotation. */
    const double p1 = 0.10180391e-4, p2 = 0.47020439e-6,
                 p3 = -0.5417367e-9, p4 = -0.2507948e-11, p5 = 0.463486e-14;
    const double q1 = -0.113469002e-3, q2 = 0.12372674e-6,
                 q3 = 0.1265417e-8, q4 = -0.1371808e-11, q5 = -0.320334e-14;

    double pw = (p1 + p2 * t + p3 * t2 + p4 * t3 + p5 * t4) * t;
    double qw = (q1 + q2 * t + q3 * t2 + q4 * t3 + q5 * t4) * t;
    const double ra = 2.0 * std::sqrt(1.0 - pw * pw - qw * qw);
    const double pwqw = 2.0 * pw * qw;
    const double pw2 = 1.0 - 2.0 * pw * pw;
    const double qw2 = 1.0 - 2.0 * qw * qw;
    pw *= ra;
    qw *= ra;

    x = pw2 * x1 + pwqw * x2 + pw * x3;
    y = pwqw * x1 + qw2 * x2 - qw * x3;
    z = -pw * x1 + qw * x2 + (pw2 + qw2 - 1.0) * x3;
}

}  // namespace eph
