#include "internal.h"

#include <cstdlib>

#include "../data/vsop87d_data.h"

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

void eval_body(const vsop87d::Body& body, double jd_tt, double& lon,
               double& lat, double& r) {
    const double T = (jd_tt - J2000) / DAYS_PER_MILLENNIUM;
    lon = eval_coordinate(body.l, body.nl, T);
    lat = eval_coordinate(body.b, body.nb, T);
    r = eval_coordinate(body.r, body.nr, T);
}

/* Truncated copies of the series (see vsop87_truncate). */
struct Truncated {
    vsop87d::Term* terms[3][6];
    vsop87d::Series series[3][6];
    vsop87d::Body body;
};
Truncated g_truncated[8];

}  // namespace

void vsop87_heliocentric(int planet, double jd_tt, double& lon, double& lat,
                         double& r) {
    if (planet < 0 || planet > 7) {
        lon = lat = r = 0.0;
        return;
    }
    eval_body(vsop87d::bodies[planet], jd_tt, lon, lat, r);
    apply_correction(planet, jd_tt, lon, lat, r);
    lon = normalize(lon);
}

void vsop87_truncate(int planet, double threshold, double tmax_days) {
    /* Mean distances (AU) turn the angular threshold into one for R. */
    static const double kAxis[8] = {0.387, 0.723, 1.0,   1.524,
                                    5.203, 9.537, 19.19, 30.07};
    const double tm = tmax_days / DAYS_PER_MILLENNIUM;
    const vsop87d::Body& full = vsop87d::bodies[planet];
    Truncated& t = g_truncated[planet];
    const vsop87d::Series* src[3] = {full.l, full.b, full.r};
    const int n[3] = {full.nl, full.nb, full.nr};
    for (int c = 0; c < 3; ++c) {
        const double thr = (c == 2) ? threshold * kAxis[planet] : threshold;
        for (int alpha = 0; alpha < n[c]; ++alpha) {
            const vsop87d::Series& s = src[c][alpha];
            vsop87d::Term*& kept = t.terms[c][alpha];
            std::free(kept);
            kept = static_cast<vsop87d::Term*>(
                std::malloc(sizeof(vsop87d::Term) * (s.count + 1)));
            const double f = std::pow(tm, alpha);
            int count = 0;
            for (int i = 0; i < s.count; ++i)
                if (kept && std::fabs(s.terms[i].a) * f >= thr)
                    kept[count++] = s.terms[i];
            t.series[c][alpha] = {s.poly, kept, count};
        }
    }
    t.body = {t.series[0], n[0], t.series[1], n[1], t.series[2], n[2]};
}

void vsop87_truncated(int planet, double jd_tt, double& lon, double& lat,
                      double& r) {
    eval_body(g_truncated[planet].body, jd_tt, lon, lat, r);
}

}  // namespace eph
