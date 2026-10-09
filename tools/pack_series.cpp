/* Pack the series coefficients (the data/ coefficient tables) into a data file.
 *
 *   pack_series series    OUT VSOP ELP   theory: VSOP87D, ELP2000-82B, Pluto
 *                                         and the DE440 corrections
 *   pack_series base      OUT VSOP ELP   base series of a residual table
 *   pack_series asteroids OUT            Ceres, Pallas, Juno, Vesta
 *
 * VSOP is one threshold (radians) or eight, comma-separated, for Mercury..
 * Neptune; ELP is one threshold.  A term is kept when its amplitude, times
 * T^alpha at the 1550/2650 ends of the fitted range, reaches the threshold
 * (scaled by the mean distance for R); 0 keeps everything and a negative
 * threshold leaves the body out.  Built by the Makefile / CMake as a
 * development tool; the format is described in src/data.h.  Assumes a
 * little-endian host (as the loader does). */
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

#include "asteroid_data.h"
#include "corr_data.h"
#include "data.h"
#include "elp_data.h"
#include "pluto_data.h"
#include "vsop87d_data.h"

using namespace eph;

namespace {

const double kTmaxDays = 237400.0; /* 2650 - J2000: the fitted range is 1550-2650 */
const double kPi = 3.14159265358979323846;

struct Section {
    uint32_t type, id;
    std::vector<unsigned char> bytes;

    template <typename T>
    void put(const T& v) {
        const unsigned char* p = reinterpret_cast<const unsigned char*>(&v);
        bytes.insert(bytes.end(), p, p + sizeof(T));
    }
    void put_bytes(const void* p, size_t n) {
        const unsigned char* q = static_cast<const unsigned char*>(p);
        bytes.insert(bytes.end(), q, q + n);
    }
};

std::vector<Section> g_sections;
long g_kept = 0, g_total = 0;

void add_vsop(int planet, double thr, uint32_t role) {
    static const double kAxis[8] = {0.387, 0.723, 1.0,   1.524,
                                    5.203, 9.537, 19.19, 30.07};
    const vsop87d::Body& b = vsop87d::bodies[planet];
    const vsop87d::Series* src[3] = {b.l, b.b, b.r};
    const uint32_t n[3] = {(uint32_t)b.nl, (uint32_t)b.nb, (uint32_t)b.nr};
    std::vector<vsop87d::Term> kept[3][6];
    for (int c = 0; c < 3; ++c) {
        const double t = (c == 2) ? thr * kAxis[planet] : thr;
        for (uint32_t a = 0; a < n[c]; ++a) {
            const double f = std::pow(kTmaxDays / 365250.0, (double)a);
            for (int i = 0; i < src[c][a].count; ++i) {
                ++g_total;
                if (std::fabs(src[c][a].terms[i].a) * f >= t) {
                    kept[c][a].push_back(src[c][a].terms[i]);
                    ++g_kept;
                }
            }
        }
    }
    Section s{kSectionVsop, (uint32_t)planet | role << 8, {}};
    s.put(n[0]);
    s.put(n[1]);
    s.put(n[2]);
    s.put(uint32_t(0));
    for (int c = 0; c < 3; ++c)
        for (uint32_t a = 0; a < n[c]; ++a) {
            s.put(src[c][a].poly);
            s.put((uint32_t)kept[c][a].size());
            s.put(uint32_t(0));
        }
    for (int c = 0; c < 3; ++c)
        for (uint32_t a = 0; a < n[c]; ++a)
            s.put_bytes(kept[c][a].data(), kept[c][a].size() * sizeof(vsop87d::Term));
    g_sections.push_back(std::move(s));
}

void add_elp(double thr, uint32_t role) {
    /* Amplitudes are arcseconds (longitude, latitude) and km (distance). */
    const double arcsec = kPi / 648000.0;
    const double t[3] = {thr / arcsec, thr / arcsec, thr * 385000.0};
    const double tm = kTmaxDays / 36525.0;
    std::vector<elp::MainTerm> main[3];
    std::vector<elp::PertTerm> pert[3];
    for (int iv = 0; iv < 3; ++iv) {
        const elp::Coord& c = elp::coords[iv];
        for (int k = 0; k < c.nmain; ++k, ++g_total)
            if (std::fabs(c.main[k].amp) >= t[iv]) {
                main[iv].push_back(c.main[k]);
                ++g_kept;
            }
        for (int k = 0; k < c.npert; ++k, ++g_total) {
            const double f = std::pow(tm, (double)(c.pert[k].packed >> 30));
            if (std::fabs(c.pert[k].amp) * f >= t[iv]) {
                pert[iv].push_back(c.pert[k]);
                ++g_kept;
            }
        }
    }
    Section s{kSectionElp, role << 8, {}};
    for (int iv = 0; iv < 3; ++iv) {
        s.put((uint32_t)main[iv].size());
        s.put((uint32_t)pert[iv].size());
    }
    for (int iv = 0; iv < 3; ++iv) {
        s.put_bytes(main[iv].data(), main[iv].size() * sizeof(elp::MainTerm));
        s.put_bytes(pert[iv].data(), pert[iv].size() * sizeof(elp::PertTerm));
    }
    g_sections.push_back(std::move(s));
}

template <typename Body>
void add_chebyshev(uint32_t type, uint32_t id, const Body& b) {
    Section s{type, id, {}};
    s.put((uint32_t)b.nseg);
    s.put((uint32_t)b.ncoef);
    s.put_bytes(b.t, (b.nseg + 1) * sizeof(double));
    s.put_bytes(b.c, (size_t)b.nseg * 3 * b.ncoef * sizeof(float));
    g_sections.push_back(std::move(s));
}

void add_correction(uint32_t body, const corr::Series& c) {
    Section s{kSectionCorrection, body, {}};
    s.put(c.tmin);
    s.put(c.tmax);
    s.put((uint32_t)c.n0);
    s.put((uint32_t)c.nh);
    s.put((uint32_t)c.nharm);
    s.put(uint32_t(0));
    s.put_bytes(c.harm, c.nharm * sizeof(corr::Harmonic));
    s.put_bytes(c.coef, 3 * (c.n0 + 2 * c.nharm * c.nh) * sizeof(float));
    g_sections.push_back(std::move(s));
}

void add_asteroids() {
    for (uint32_t i = 0; i < kAsteroids; ++i)
        add_chebyshev(kSectionAsteroid, i, ast::bodies[i]);
}

void add_info(const std::string& text) {
    Section s{kSectionInfo, 0, {}};
    s.put_bytes(text.data(), text.size());
    g_sections.push_back(std::move(s));
}

bool write(const char* path) {
    std::vector<unsigned char> out;
    auto put32 = [&](uint32_t v) {
        for (int i = 0; i < 4; ++i) out.push_back((v >> (8 * i)) & 0xFF);
    };
    auto put64 = [&](uint64_t v) {
        for (int i = 0; i < 8; ++i) out.push_back((v >> (8 * i)) & 0xFF);
    };
    out.insert(out.end(), {'E', 'P', 'H', 'D'});
    put32(2);
    put32((uint32_t)g_sections.size());
    put32(0);
    uint64_t off = 16 + 24 * g_sections.size();
    for (const Section& s : g_sections) {
        off = (off + 7) & ~uint64_t(7);
        put32(s.type);
        put32(s.id);
        put64(off);
        put64(s.bytes.size());
        off += s.bytes.size();
    }
    for (const Section& s : g_sections) {
        while (out.size() % 8) out.push_back(0);
        out.insert(out.end(), s.bytes.begin(), s.bytes.end());
    }
    uint32_t crc = 0xFFFFFFFFu;
    for (size_t i = 16; i < out.size(); ++i) {
        crc ^= out[i];
        for (int k = 0; k < 8; ++k) crc = (crc & 1) ? 0xEDB88320u ^ (crc >> 1) : crc >> 1;
    }
    crc ^= 0xFFFFFFFFu;
    for (int i = 0; i < 4; ++i) out[12 + i] = (crc >> (8 * i)) & 0xFF;
    FILE* f = std::fopen(path, "wb");
    if (!f) return false;
    const bool ok = std::fwrite(out.data(), 1, out.size(), f) == out.size();
    return std::fclose(f) == 0 && ok;
}

/* "t" or "t0,t1,...,t7" -> eight thresholds. */
bool parse_thresholds(const char* arg, double t[8]) {
    std::string s = arg;
    int n = 0;
    size_t pos = 0;
    while (n < 8) {
        const size_t c = s.find(',', pos);
        t[n++] = std::atof(s.substr(pos, c - pos).c_str());
        if (c == std::string::npos) break;
        pos = c + 1;
    }
    if (n != 1 && n != 8) return false;
    for (int i = n; i < 8; ++i) t[i] = t[0];
    return true;
}

}  // namespace

int main(int argc, char** argv) {
    const uint32_t one = 1;
    if (*reinterpret_cast<const unsigned char*>(&one) != 1) {
        std::fprintf(stderr, "pack_series: needs a little-endian host\n");
        return 1;
    }
    const std::string mode = argc > 2 ? argv[1] : "";
    double vsop[8];
    if (mode == "asteroids" && argc == 3) {
        add_info("asteroids: Ceres, Pallas, Juno, Vesta, Chiron (fit to JPL Horizons)");
        add_asteroids();
    } else if ((mode == "series" || mode == "base") && argc == 5 &&
               parse_thresholds(argv[3], vsop)) {
        const bool theory = (mode == "series");
        const double elp = std::atof(argv[4]);
        add_info(std::string(theory ? "series" : "table base") + ": VSOP87D " +
                 argv[3] + ", ELP2000-82B " + argv[4]);
        for (int p = 0; p < 8; ++p)
            if (vsop[p] >= 0.0) add_vsop(p, vsop[p], theory ? kTheory : kBase);
        if (elp >= 0.0) add_elp(elp, theory ? kTheory : kBase);
        if (theory) {
            add_chebyshev(kSectionPluto, 0, pluto_fit::body);
            const corr::Series* c[8] = {&corr::mercury, &corr::venus,
                                        &corr::earth,   &corr::mars,
                                        &corr::jupiter, &corr::saturn,
                                        &corr::uranus,  &corr::neptune};
            for (uint32_t p = 0; p < 8; ++p) add_correction(p, *c[p]);
            add_correction(10, corr::moon);
        }
    } else {
        std::fprintf(stderr,
                     "usage: pack_series series|base OUT VSOP ELP\n"
                     "       pack_series asteroids OUT\n");
        return 1;
    }
    if (!write(argv[2])) {
        std::fprintf(stderr, "pack_series: cannot write %s\n", argv[2]);
        return 1;
    }
    if (g_total)
        std::fprintf(stderr, "%s: %ld of %ld series terms\n", argv[2], g_kept,
                     g_total);
    return 0;
}
