#include "internal.h"

#include "../data/pluto_data.h"

namespace eph {

/*
 * Pluto heliocentric position in the J2000 ecliptic frame.
 *
 * The coefficients in data/pluto_data.h are produced by tools/fit_pluto.py,
 * a least-squares Poisson-series fit to JPL DE440s over 1850-2150 with the
 * argument structure
 *
 *     arg = i*(34.35 + 3034.9057 T) + j*(50.08 + 1222.1138 T)
 *             + k*(238.96 + 144.96 T)        [degrees, T in Julian centuries]
 *
 * Sub-arcsecond heliocentric accuracy across the fitted window.
 */
void pluto_heliocentric_j2000(double jd_tt, double& lon, double& lat, double& r) {
    const double T = (jd_tt - J2000) / DAYS_PER_CENTURY;

    const double J = 34.35 + 3034.9057 * T;
    const double S = 50.08 + 1222.1138 * T;
    const double P = 238.96 + 144.96 * T;

    auto eval = [&](const pluto_fit::Model& m) {
        double v = m.a0 + m.a1 * T;
        for (int n = 0; n < m.count; ++n) {
            const pluto_fit::Term& t = m.terms[n];
            const double arg = (t.i * J + t.j * S + t.k * P) * DEG;
            v += t.s * std::sin(arg) + t.c * std::cos(arg);
        }
        return v;
    };

    lon = normalize(eval(pluto_fit::lon) * DEG);
    lat = eval(pluto_fit::lat) * DEG;
    r = eval(pluto_fit::rad);
}

}  // namespace eph
