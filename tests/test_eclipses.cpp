/* Eclipses against NASA's Five Millennium Canon 1901-2100 (types,
 * greatest eclipse, gamma, magnitudes, lunar phase durations;
 * tests/eclipse_nasa_ref.h) and Swiss Ephemeris 1850-2150 (contact times,
 * local circumstances; tests/eclipse_ref.h). */
#include <cmath>
#include <cstdio>

#include "eclipse_nasa_ref.h"
#include "eclipse_ref.h"
#include "ephemeris/ephemeris.h"
#include "load_data.h"

namespace {

const double kDeg = M_PI / 180.0;
int failures = 0;
double worst[32];
const char* names[32];

/* Track the worst difference of a quantity; fail beyond `limit`. */
void check(double jd, int slot, const char* name, double ours, double ref,
           double limit, double scale = 1.0) {
    if (ref == 0.0 && ours == 0.0) return;
    names[slot] = name;
    const double d = (ref == 0.0 || ours == 0.0) ? 1e9 : std::fabs(ours - ref) * scale;
    if (d > worst[slot]) worst[slot] = d;
    if (d > limit) {
        ++failures;
        if (failures < 20)
            std::printf("FAIL %s at %.5f: %.6f vs %.6f\n", name, jd, ours, ref);
    }
}

}  // namespace

int main(int argc, char** argv) {
    if (!load_data_args(argc, argv)) return 1;
    const double sec = 86400.0;

    for (const NasaSolar& r : kNasaSolar) {
        eph_solar_eclipse e;
        if (eph_solar_eclipse_next(r.max - 5, 1, &e) != 0) {
            std::printf("FAIL: no solar eclipse after %.1f\n", r.max - 5);
            return 1;
        }
        if (e.type != r.type || e.saros != r.saros) {
            ++failures;
            std::printf("FAIL NASA solar at %.5f: type %d saros %d vs %d %d\n", r.max,
                        e.type, e.saros, r.type, r.saros);
        }
        check(r.max, 16, "NASA solar maximum (s)", e.maximum, r.max, 10, sec);
        check(r.max, 17, "NASA solar gamma", e.gamma, r.gamma, 0.0006);
        check(r.max, 18, "NASA solar magnitude", e.magnitude, r.mag, 0.0006);
    }
    for (const NasaLunar& r : kNasaLunar) {
        eph_lunar_eclipse e;
        if (eph_lunar_eclipse_next(r.max - 5, 1, &e) != 0) {
            std::printf("FAIL: no lunar eclipse after %.1f\n", r.max - 5);
            return 1;
        }
        if (e.type != r.type || e.saros != r.saros) {
            ++failures;
            std::printf("FAIL NASA lunar at %.5f: type %d saros %d vs %d %d\n", r.max,
                        e.type, e.saros, r.type, r.saros);
        }
        check(r.max, 19, "NASA lunar maximum (s)", e.maximum, r.max, 30, sec);
        check(r.max, 20, "NASA lunar gamma", e.gamma, r.gamma, 0.001);
        check(r.max, 21, "NASA lunar magnitudes", e.penumbral_magnitude, r.pmag, 0.005);
        check(r.max, 21, "NASA lunar magnitudes", e.umbral_magnitude, r.umag, 0.005);
        const double min = 1440.0;
        check(r.max, 22, "NASA lunar durations (min)",
              e.penumbral_end - e.penumbral_begin, r.pen / min, 1.0, min);
        if (r.par > 0)
            check(r.max, 22, "NASA lunar durations (min)",
                  e.partial_end - e.partial_begin, r.par / min, 1.0, min);
        if (r.total > 0)
            check(r.max, 22, "NASA lunar durations (min)",
                  e.total_end - e.total_begin, r.total / min, 1.0, min);
    }

    /* Backward searches find the same eclipses; the axis point at
     * greatest eclipse is where the eclipse is greatest. */
    for (size_t i = 0; i < sizeof kNasaSolar / sizeof *kNasaSolar; i += 10) {
        const NasaSolar& r = kNasaSolar[i];
        eph_solar_eclipse e;
        eph_lunar_eclipse l;
        if (eph_solar_eclipse_next(r.max + 5, -1, &e) != 0 ||
            eph_lunar_eclipse_next(kNasaLunar[i].max + 5, -1, &l) != 0)
            return 1;
        check(r.max, 24, "backward search (s)", e.maximum, r.max, 10, sec);
        check(r.max, 24, "backward search (s)", l.maximum, kNasaLunar[i].max, 30, sec);
        double lat, lon;
        if (!(e.type & EPH_ECLIPSE_NONCENTRAL) && e.type != EPH_ECLIPSE_PARTIAL) {
            if (eph_solar_eclipse_where(e.maximum, &lat, &lon) != 0) {
                ++failures;
                std::printf("FAIL where at %.5f\n", r.max);
            }
            check(r.max, 25, "where (deg)", lat, e.latitude, 1e-6, 1 / kDeg);
            check(r.max, 25, "where (deg)", lon, e.longitude, 1e-6, 1 / kDeg);
        }
    }

    /* Contact durations from NASA's eclipse figures (Espenak): 2024 Apr 8
     * P1-P4 5h10m06.8s, U1-U4 3h16m44.7s; 2025 Sep 21 P1-P4 4h24m01.5s. */
    {
        eph_solar_eclipse a, b;
        if (eph_solar_eclipse_next(2460400.0, 1, &a) != 0 ||
            eph_solar_eclipse_next(2460930.0, 1, &b) != 0)
            return 1;
        check(a.maximum, 23, "NASA figure durations (s)", a.partial_end - a.partial_begin,
              18606.8 / sec, 3, sec);
        check(a.maximum, 23, "NASA figure durations (s)", a.central_end - a.central_begin,
              11804.7 / sec, 3, sec);
        check(b.maximum, 23, "NASA figure durations (s)", b.partial_end - b.partial_begin,
              15841.5 / sec, 3, sec);
    }

    for (const SolarRef& r : kSolar) {
        eph_solar_eclipse e;
        if (eph_solar_eclipse_next(r.max - 5, 1, &e) != 0) {
            std::printf("FAIL: no solar eclipse after %.1f\n", r.max - 5);
            return 1;
        }
        /* Swiss Ephemeris' global contacts differ from NASA's by up to
         * ~25 s (see the figures below), more where the shadow only grazes
         * the Earth and the contacts are ill-conditioned. */
        if (std::fabs(e.gamma) < 1.45) {
            check(r.max, 1, "solar begin/end (s)", e.partial_begin, r.begin, 60, sec);
            check(r.max, 1, "solar begin/end (s)", e.partial_end, r.end, 60, sec);
        }
        if (std::fabs(e.gamma) < 0.99) {
            check(r.max, 2, "solar central (s)", e.central_begin, r.central_begin, 60, sec);
            check(r.max, 2, "solar central (s)", e.central_end, r.central_end, 60, sec);
            check(r.max, 3, "solar center line (s)", e.center_line_begin, r.center_begin,
                  60, sec);
            check(r.max, 3, "solar center line (s)", e.center_line_end, r.center_end, 60,
                  sec);
        }
    }

    for (const LunarRef& r : kLunar) {
        eph_lunar_eclipse e;
        if (eph_lunar_eclipse_next(r.max - 5, 1, &e) != 0) {
            std::printf("FAIL: no lunar eclipse after %.1f\n", r.max - 5);
            return 1;
        }
        check(r.max, 5, "lunar maximum (s)", e.maximum, r.max, 30, sec);
        check(r.max, 6, "lunar penumbral (s)", e.penumbral_begin, r.pen_begin, 60, sec);
        check(r.max, 6, "lunar penumbral (s)", e.penumbral_end, r.pen_end, 60, sec);
        check(r.max, 7, "lunar partial (s)", e.partial_begin, r.par_begin, 60, sec);
        check(r.max, 7, "lunar partial (s)", e.partial_end, r.par_end, 60, sec);
        check(r.max, 8, "lunar total (s)", e.total_begin, r.tot_begin, 60, sec);
        check(r.max, 8, "lunar total (s)", e.total_end, r.tot_end, 60, sec);
    }

    for (const LocalRef& r : kLocal) {
        const eph_observer o = {kPlaces[r.place][0] * kDeg, kPlaces[r.place][1] * kDeg,
                                kPlaces[r.place][2]};
        eph_local_solar_eclipse e;
        if (eph_solar_eclipse_local(r.max - 5, &o, 1, &e) != 0) {
            std::printf("FAIL: no local eclipse after %.1f\n", r.max - 5);
            return 1;
        }
        /* Swiss Ephemeris gives the maximum while the Sun is up. */
        if (r.altitude > 0.0) {
            check(r.max, 11, "local maximum (s)", e.maximum, r.max, 10, sec);
            check(r.max, 14, "local magnitude", e.magnitude, r.mag, 0.002);
            check(r.max, 15, "local obscuration", e.obscuration, r.obscuration, 0.002);
        }
        check(r.max, 12, "local 1st/4th contact (s)", e.first_contact, r.c1, 10, sec);
        check(r.max, 12, "local 1st/4th contact (s)", e.fourth_contact, r.c4, 10, sec);
        check(r.max, 13, "local 2nd/3rd contact (s)", e.second_contact, r.c2, 10, sec);
        check(r.max, 13, "local 2nd/3rd contact (s)", e.third_contact, r.c3, 10, sec);
    }

    for (int i = 0; i < 32; ++i)
        if (names[i]) std::printf("  %-28s worst %.4f\n", names[i], worst[i]);
    std::printf("%s (%zu + %zu solar, %zu + %zu lunar, %zu local eclipses).\n",
                failures ? "Eclipse checks FAILED" : "All eclipse checks passed", sizeof kNasaSolar / sizeof *kNasaSolar,
                sizeof kSolar / sizeof *kSolar, sizeof kNasaLunar / sizeof *kNasaLunar,
                sizeof kLunar / sizeof *kLunar, sizeof kLocal / sizeof *kLocal);
    return failures ? 1 : 0;
}
