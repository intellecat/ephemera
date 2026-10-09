/* Data files: the core without data, loading from memory and files,
 * rejection of damaged files, and the accuracy of each series level.
 * Usage: test_data <data directory> */
#include <cmath>
#include <cstdio>
#include <algorithm>
#include <cstdlib>
#include <string>
#include <vector>

#include "ephemeris/ephemeris.h"
#include "reference_data.h"

static const double R2D = 180.0 / 3.14159265358979323846;
static int g_fails = 0;

static void check(bool ok, const char* what) {
    if (!ok) {
        std::printf("FAIL: %s\n", what);
        g_fails++;
    }
}

static std::vector<unsigned char> read_file(const std::string& path) {
    std::vector<unsigned char> buf;
    FILE* f = std::fopen(path.c_str(), "rb");
    if (!f) return buf;
    int c;
    while ((c = std::fgetc(f)) != EOF) buf.push_back((unsigned char)c);
    std::fclose(f);
    return buf;
}

static double wrap180(double d) {
    d = std::fmod(d + 180.0, 360.0);
    return (d < 0.0 ? d + 360.0 : d) - 180.0;
}

/* Max apparent error (arcsec) over the reference samples, per body. */
static void max_errors(double err[EPH_BODY_COUNT]) {
    for (int i = 0; i < EPH_BODY_COUNT; ++i) err[i] = 0.0;
    for (const RefSample& r : kRefSamples) {
        eph_position a;
        if (eph_position_apparent_tt((eph_body)r.body, r.tt_jd, &a) != 0) {
            err[r.body] = 1e9;
            continue;
        }
        const double e = std::fmax(
            std::fabs(wrap180(a.longitude * R2D - r.app_lon)),
            std::fabs(a.latitude * R2D - r.app_lat)) * 3600.0;
        err[r.body] = std::fmax(err[r.body], e);
    }
}

int main(int argc, char** argv) {
    if (argc != 2) {
        std::printf("usage: %s <data directory>\n", argv[0]);
        return 1;
    }
    const std::string dir = std::string(argv[1]) + "/";
    eph_position p;

    /* A build with EPHEMERIS_EMBED_DATA works without any load call;
     * eph_load_embedded() reloads it after an unload. */
    const bool embedded = eph_position_tt(EPH_SUN, 2451545.0, &p) == 0;
    std::printf("embedded data: %s\n", embedded ? "yes" : "no");
    eph_unload_data();
    check(eph_load_embedded() == (embedded ? 0 : -1), "eph_load_embedded");
    if (embedded) {
        check(eph_position_apparent_tt(EPH_MOON, 2451545.0, &p) == 0,
              "embedded data reloads");
        eph_unload_data();
    }

    /* The core alone: no positions, but time and houses work. */
    check(eph_position_tt(EPH_SUN, 2451545.0, &p) != 0, "position without data");
    check(std::fabs(eph_julian_date(2000, 1, 1.5) - 2451545.0) < 1e-9,
          "julian date without data");
    double cusps[12];
    eph_observer obs = {0.1, 0.9, 0.0};
    check(eph_houses(2451545.0, &obs, EPH_HOUSES_EQUAL, cusps) == 0,
          "houses without data");

    /* Damaged files are rejected and leave nothing loaded. */
    const std::vector<unsigned char> full = read_file(dir + "eph-full.bin");
    check(!full.empty(), "read eph-full.bin");
    for (size_t n : {size_t(0), size_t(15), size_t(100), full.size() / 2,
                     full.size() - 1}) {
        check(eph_load_data(full.data(), n) != 0, "truncated file rejected");
    }
    std::vector<unsigned char> bad = full;
    std::srand(1);
    for (int k = 0; k < 200 && !bad.empty(); ++k) {
        bad = full;
        for (int j = 0; j < 4; ++j) bad[16 + std::rand() % 400] = (unsigned char)std::rand();
        if (eph_load_data(bad.data(), bad.size()) == 0) eph_unload_data();
    }
    check(eph_position_tt(EPH_SUN, 2451545.0, &p) != 0, "nothing left loaded");

    /* In place: needs 8-byte alignment. */
    std::vector<double> aligned((full.size() + 15) / 8);
    unsigned char* a = reinterpret_cast<unsigned char*>(aligned.data());
    std::copy(full.begin(), full.end(), a);
    check(eph_load_static(a + 1, full.size()) != 0, "misaligned static data");
    check(eph_load_static(a, full.size()) == 0, "static data");
    check(eph_position_tt(EPH_SUN, 2451545.0, &p) == 0, "position from static");
    eph_unload_data();

    /* From memory, then unload. */
    check(eph_load_data(full.data(), full.size()) == 0, "load from memory");
    check(eph_position_tt(EPH_SUN, 2451545.0, &p) == 0, "position after load");
    check(eph_position_tt(EPH_CERES, 2451545.0, &p) != 0,
          "no asteroids without their file");
    eph_unload_data();
    check(eph_position_tt(EPH_SUN, 2451545.0, &p) != 0, "unload");

    /* The residual table alone covers 1850-2140 and nothing else. */
    check(eph_load_file((dir + "ephem.bin").c_str()) == 0, "load ephem.bin");
    check(eph_position_tt(EPH_MOON, 2451545.0, &p) == 0, "table inside range");
    check(eph_position_tt(EPH_MOON, 2300000.0, &p) != 0, "table outside range");
    eph_unload_data();

    /* Each series level within its accuracy (max over the reference,
     * arcsec; the Moon is ~0.45" at the 1850/2140 ends at every level). */
    struct {
        const char* file;
        double planets, moon;
    } levels[] = {{"eph-full.bin", 0.15, 0.5},
                  {"eph-01as.bin", 0.2, 0.6},
                  {"eph-1as.bin", 1.0, 1.0}};
    for (const auto& l : levels) {
        if (eph_load_file((dir + l.file).c_str()) != 0) {
            std::printf("FAIL: cannot load %s\n", l.file);
            g_fails++;
            continue;
        }
        double err[EPH_BODY_COUNT];
        max_errors(err);
        double worst = 0.0;
        for (int b = 0; b < EPH_BODY_COUNT; ++b) {
            if (b == EPH_MOON || b == EPH_EARTH) continue;
            worst = std::fmax(worst, err[b]);
        }
        std::printf("%-14s max %.3f\" (planets, Sun, Pluto)  %.3f\" (Moon)\n",
                    l.file, worst, err[EPH_MOON]);
        check(worst <= l.planets && err[EPH_MOON] <= l.moon, l.file);
        eph_unload_data();
    }

    if (g_fails == 0) std::printf("\nAll data checks passed.\n");
    return g_fails == 0 ? 0 : 1;
}
