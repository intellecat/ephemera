#include "internal.h"

#include "../data/vsop87d_data.h"

namespace eph {

namespace {

/* Sum A*cos(B + C*T) over one series. */
double series_sum(const vsop87d::Series& s, double T) {
    double sum = 0.0;
    const vsop87d::Term* t = s.terms;
    for (int i = 0; i < s.count; ++i, ++t) {
        sum += t->a * std::cos(t->b + t->c * T);
    }
    return sum;
}

/* Evaluate a Poisson series: sum_alpha T^alpha * series_alpha. */
double eval_coordinate(const vsop87d::Series* arr, int n, double T) {
    double val = 0.0;
    for (int alpha = n - 1; alpha >= 0; --alpha) {
        val = val * T + series_sum(arr[alpha], T);
    }
    return val;
}

}  // namespace

void vsop87_heliocentric(int planet, double jd_tt, double& lon, double& lat,
                         double& r) {
    if (planet < 0 || planet > 7) {
        lon = lat = r = 0.0;
        return;
    }
    const double T = (jd_tt - J2000) / DAYS_PER_MILLENNIUM;
    const vsop87d::Body& body = vsop87d::bodies[planet];
    lon = normalize(eval_coordinate(body.l, body.nl, T));
    lat = eval_coordinate(body.b, body.nb, T);
    r = eval_coordinate(body.r, body.nr, T);
}

}  // namespace eph
