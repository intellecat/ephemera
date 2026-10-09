#include "internal.h"

#include "data.h"

namespace eph {

namespace {

/* Sum A*cos(B + C*T) over one series. */
double series_sum(const vsop87d::Series& s, double T) {
    double sum = 0.0;
    const vsop87d::Term* t = s.terms;
    const double bscale = TWO_PI / 4294967296.0; /* uint32 phase -> radians */
    for (int i = 0; i < s.count; ++i, ++t) {
        sum += t->a * std::cos(t->b * bscale + t->c * T);
    }
    return sum;
}

/* Evaluate a Poisson series: sum_alpha T^alpha * series_alpha. */
double eval_coordinate(const vsop87d::Series* arr, int n, double T) {
    double val = 0.0;
    for (int alpha = n - 1; alpha >= 0; --alpha) {
        val = val * T + arr[alpha].poly + series_sum(arr[alpha], T);
    }
    return val;
}

bool eval_body(int role, int planet, double jd_tt, double& lon, double& lat,
               double& r) {
    const SeriesSet& s = data().set[role];
    if (planet < 0 || planet > 7 || !s.vsop_ok[planet]) return false;
    const vsop87d::Body& body = s.vsop[planet];
    const double T = (jd_tt - J2000) / DAYS_PER_MILLENNIUM;
    lon = eval_coordinate(body.l, body.nl, T);
    lat = eval_coordinate(body.b, body.nb, T);
    r = eval_coordinate(body.r, body.nr, T);
    return true;
}

}  // namespace

bool vsop87_heliocentric(int planet, double jd_tt, double& lon, double& lat,
                         double& r) {
    if (!eval_body(kTheory, planet, jd_tt, lon, lat, r)) return false;
    apply_correction(planet, jd_tt, lon, lat, r);
    lon = normalize(lon);
    return true;
}

bool vsop87_base(int planet, double jd_tt, double& lon, double& lat,
                 double& r) {
    return eval_body(kBase, planet, jd_tt, lon, lat, r);
}

}  // namespace eph
