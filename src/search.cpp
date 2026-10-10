#include "internal.h"

#include "../include/ephemeris/ephemeris.h"

/* Event searches: longitude crossings (returns, ingresses), node passages
 * and stations.  Each steps through time until a quantity changes sign,
 * then refines the bracket by false position (Illinois variant). */

namespace {

/* Step (days) short enough that no two sign changes fit in one, even
 * around a planet's station. */
double search_step(eph_body body) {
    switch (body) {
        case EPH_MOON: return 0.5;
        case EPH_TRUE_NODE: return 0.25;
        case EPH_OSCU_LILITH: return 0.1;
        case EPH_MEAN_NODE:
        case EPH_MEAN_LILITH: return 5.0;
        default: return 1.0;
    }
}

/* Upper bounds (degrees/day, with margin) on how fast the apparent
 * geocentric longitude can change: far from a target, the search may step
 * as far as the body could possibly move. */
double max_longitude_rate(eph_body body) {
    switch (body) {
        case EPH_SUN: return 1.1;
        case EPH_MOON: return 16.0;
        case EPH_MERCURY: return 2.5;
        case EPH_VENUS: return 1.4;
        case EPH_MARS: return 0.9;
        case EPH_JUPITER: return 0.3;
        case EPH_SATURN: return 0.15;
        case EPH_URANUS: return 0.08;
        case EPH_NEPTUNE:
        case EPH_PLUTO: return 0.05;
        case EPH_CHIRON: return 0.2;
        case EPH_MEAN_NODE: return 0.06;
        case EPH_TRUE_NODE: return 0.6;
        case EPH_MEAN_LILITH: return 0.15;
        case EPH_OSCU_LILITH: return 20.0;
        case EPH_CUPIDO: case EPH_HADES: case EPH_ZEUS: case EPH_KRONOS:
        case EPH_APOLLON: case EPH_ADMETOS: case EPH_VULKANUS: case EPH_POSEIDON:
            return 0.03; /* the Earth's parallax dominates */
        default: return 2.0; /* asteroids */
    }
}

/* The same for the heliocentric longitude. */
double max_heliocentric_rate(eph_body body) {
    switch (body) {
        case EPH_MERCURY: return 6.5;
        case EPH_VENUS: return 1.7;
        case EPH_MARS: return 0.7;
        case EPH_JUPITER: return 0.1;
        case EPH_SATURN: return 0.04;
        case EPH_URANUS: return 0.015;
        case EPH_NEPTUNE:
        case EPH_PLUTO: return 0.01;
        case EPH_CHIRON: return 0.05;
        case EPH_CUPIDO: case EPH_HADES: case EPH_ZEUS: case EPH_KRONOS:
        case EPH_APOLLON: case EPH_ADMETOS: case EPH_VULKANUS: case EPH_POSEIDON:
            return 0.001;
        default: return 0.5; /* asteroids */
    }
}

/* Signed angle a - b in (-pi, pi]. */
double wrap(double a) { return eph::normalize(a + eph::PI) - eph::PI; }

/* First zero of f after (direction +1) or before (-1) jd within max_days.
 * f(jd, value) returns false when it cannot be evaluated.  `angular`
 * values are wrapped angles: a jump of more than pi between two samples is
 * the wrap, not a crossing.  Steps are at least `step` days, and longer
 * while |value| / max_rate (radians/day; 0: unknown) says the value cannot
 * reach zero sooner. */
template <typename F>
bool find_zero(F f, double jd, int direction, double max_days, double step,
               double max_rate, bool angular, double& found) {
    if (direction != 1 && direction != -1) return false;
    if (!(max_days > 0.0)) return false;
    double t0 = jd, v0;
    if (!f(t0, v0)) return false;
    /* A zero exactly at the start does not count. */
    if (v0 == 0.0) {
        t0 = jd + direction * 1e-6;
        if (!f(t0, v0)) return false;
    }
    for (double s = 0.0;;) {
        s += max_rate > 0.0 ? std::fmax(step, std::fabs(v0) / max_rate) : step;
        const double t1 = jd + direction * std::fmin(s, max_days);
        double v1;
        if (!f(t1, v1)) return false;
        const bool sign_change = (v0 < 0.0) != (v1 < 0.0) || v1 == 0.0;
        if (sign_change && !(angular && std::fabs(v1 - v0) > eph::PI)) {
            /* Refine [t0, t1] (false position, Illinois). */
            double a = t0, fa = v0, b = t1, fb = v1;
            int side = 0;
            for (int i = 0; i < 100 && std::fabs(b - a) > 1e-9; ++i) {
                const double c = (fa * b - fb * a) / (fa - fb);
                double fc;
                if (!f(c, fc)) return false;
                if (fc == 0.0) {
                    a = b = c;
                    break;
                }
                if ((fc < 0.0) == (fb < 0.0)) {
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
            found = (a + b) / 2;
            return true;
        }
        if (s >= max_days) return false;
        t0 = t1;
        v0 = v1;
    }
}

double year_of(double jd) {
    int y, m;
    double d;
    eph_calendar_date(jd, &y, &m, &d);
    return y + (jd - eph_julian_date(y, 1, 1.0)) / 365.25;
}

}  // namespace

extern "C" {

int eph_longitude_crossing(eph_body body, double longitude, double jd_tt,
                           int direction, double max_days, double* jd_out) {
    if (!jd_out) return -1;
    auto f = [&](double t, double& v) {
        eph_position p;
        if (eph_position_apparent_tt(body, t, &p) != 0) return false;
        v = wrap(p.longitude - longitude);
        return true;
    };
    return find_zero(f, jd_tt, direction, max_days, search_step(body),
                     max_longitude_rate(body) * eph::DEG, true, *jd_out)
               ? 0
               : -1;
}

int eph_longitude_crossing_ut(eph_body body, double longitude, double jd_ut,
                              int direction, double max_days, double* jd_out) {
    if (!jd_out) return -1;
    double tt;
    if (eph_longitude_crossing(body, longitude, eph_ut_to_tt(jd_ut, year_of(jd_ut)),
                               direction, max_days, &tt) != 0)
        return -1;
    *jd_out = eph_tt_to_ut(tt, year_of(tt));
    return 0;
}

int eph_heliocentric_crossing(eph_body body, double longitude, double jd_tt,
                              int direction, double max_days, double* jd_out) {
    if (!jd_out) return -1;
    auto f = [&](double t, double& v) {
        double lon;
        if (eph_heliocentric_tt(body, t, &lon, nullptr, nullptr) != 0) return false;
        v = wrap(lon - longitude);
        return true;
    };
    return find_zero(f, jd_tt, direction, max_days, 1.0,
                     max_heliocentric_rate(body) * eph::DEG, true, *jd_out)
               ? 0
               : -1;
}

int eph_node_crossing(eph_body body, double jd_tt, int direction,
                      double max_days, double* jd_out, double* longitude) {
    if (!jd_out) return -1;
    auto f = [&](double t, double& v) {
        eph_position p;
        if (eph_position_apparent_tt(body, t, &p) != 0) return false;
        v = p.latitude;
        return true;
    };
    /* Latitude rates: the Moon ~1.2 deg/day, the bodies far less. */
    const double rate = (body == EPH_MOON ? 1.5 : 1.0) * eph::DEG;
    if (!find_zero(f, jd_tt, direction, max_days, search_step(body), rate, false,
                   *jd_out))
        return -1;
    eph_position p;
    if (eph_position_apparent_tt(body, *jd_out, &p) != 0) return -1;
    if (longitude) *longitude = p.longitude;
    /* 1: ascending node (latitude increasing), 2: descending. */
    return p.latitude_speed >= 0.0 ? 1 : 2;
}

int eph_station(eph_body body, double jd_tt, int direction, double max_days,
                double* jd_out) {
    if (!jd_out) return -1;
    auto f = [&](double t, double& v) {
        eph_position p;
        if (eph_position_apparent_tt(body, t, &p) != 0) return false;
        v = p.longitude_speed;
        return true;
    };
    if (!find_zero(f, jd_tt, direction, max_days, search_step(body), 0.0, false,
                   *jd_out))
        return -1;
    /* 1: turns retrograde, 2: turns direct. */
    eph_position after;
    if (eph_position_apparent_tt(body, *jd_out + 0.01, &after) != 0) return -1;
    return after.longitude_speed < 0.0 ? 1 : 2;
}

}  // extern "C"
