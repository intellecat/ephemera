#include "internal.h"

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "../include/ephemeris/ephemeris.h"

namespace eph {

/*
 * Residual table (data/ephem.bin, written by tools/pack_ephemeris.py).
 *
 * For each body the file stores DE440 minus a base series -- the library's
 * own VSOP87D / ELP2000-82B series with the terms below a threshold
 * dropped, or nothing -- as Chebyshev polynomials over fixed-length
 * segments.  The residual is small and smooth, so its coefficients are
 * quantized with one absolute step per body and packed with a fixed bit
 * width per coefficient.  Positions are in the IAU-1976 ecliptic of date,
 * AU; heliocentric, geocentric for the Moon.  The Earth's entry holds the
 * Earth-Moon barycentre, which is smooth (the Earth is derived from it and
 * the Moon).
 *
 * Layout (little-endian):
 *   0  "EPHP", u32 version (3), u32 body count, f64 jd0, f64 jd1
 *  28  per body, 160 bytes:
 *        u32 body, u32 base (0 none, 1 series), u32 nseg, u32 ncoef,
 *        f64 threshold (rad), f64 start, f64 segment length (days),
 *        f64 quantum (AU), u64 data offset, u32 bits per segment, u32 0,
 *        u8 bits[3][32] (coefficient widths, x/y/z)
 *  data: per segment, x then y then z coefficients, two's complement,
 *        LSB-first bit stream; 8 bytes of padding at the end.
 */

namespace {

constexpr int kMaxCoef = 32;
constexpr long kHeader = 28, kRecord = 160;

struct Entry {
    int body;
    bool base;
    int nseg, ncoef;
    double threshold, start, length, quantum;
    const unsigned char* data;
    int64_t segbits;
    unsigned char bits[3][kMaxCoef];
};

unsigned char* g_data = nullptr;
Entry g_entries[16];
int g_nentries = 0;
double g_jd0 = 0.0, g_jd1 = 0.0;
int g_generation = 1; /* bumped on every load/unload */

template <typename T>
T load_at(const unsigned char* p, long off) {
    T v;
    std::memcpy(&v, p + off, sizeof(T));
    return v;
}

double chebyshev(const double* c, int n, double u) {
    double b1 = 0.0, b2 = 0.0;
    for (int k = n - 1; k > 0; --k) {
        const double b0 = 2.0 * u * b1 - b2 + c[k];
        b2 = b1;
        b1 = b0;
    }
    return u * b1 - b2 + c[0];
}

/* Derivative of a Chebyshev series with respect to its argument u. */
void chebyshev_derivative(const double* c, int n, double* d) {
    d[n - 1] = 0.0;
    if (n > 1) d[n - 2] = 2.0 * (n - 1) * c[n - 1];
    for (int k = n - 3; k >= 0; --k) d[k] = d[k + 2] + 2.0 * (k + 1) * c[k + 1];
    d[0] *= 0.5;
}

/* Signed integer of `width` bits at bit offset `pos`. */
int64_t read_bits(const unsigned char* p, int64_t pos, int width) {
    if (width == 0) return 0;
    const uint64_t w = load_at<uint64_t>(p, pos >> 3) >> (pos & 7);
    const int shift = 64 - width;
    return (int64_t)(w << shift) >> shift;
}

const Entry* find_entry(int body) {
    for (int i = 0; i < g_nentries; ++i)
        if (g_entries[i].body == body) return &g_entries[i];
    return nullptr;
}

bool parse(long size) {
    if (size < kHeader || std::memcmp(g_data, "EPHP", 4) != 0) return false;
    const uint32_t version = load_at<uint32_t>(g_data, 4);
    const uint32_t n = load_at<uint32_t>(g_data, 8);
    if (version != 3 || n > 16 || kHeader + n * kRecord > size) return false;
    g_jd0 = load_at<double>(g_data, 12);
    g_jd1 = load_at<double>(g_data, 20);
    const double tmax = std::fmax(std::fabs(g_jd0 - J2000), std::fabs(g_jd1 - J2000));
    for (uint32_t i = 0; i < n; ++i) {
        const long r = kHeader + i * kRecord;
        Entry& e = g_entries[i];
        e.body = load_at<uint32_t>(g_data, r);
        e.base = load_at<uint32_t>(g_data, r + 4) != 0;
        e.nseg = load_at<uint32_t>(g_data, r + 8);
        e.ncoef = load_at<uint32_t>(g_data, r + 12);
        e.threshold = load_at<double>(g_data, r + 16);
        e.start = load_at<double>(g_data, r + 24);
        e.length = load_at<double>(g_data, r + 32);
        e.quantum = load_at<double>(g_data, r + 40);
        const uint64_t off = load_at<uint64_t>(g_data, r + 48);
        e.segbits = load_at<uint32_t>(g_data, r + 56);
        std::memcpy(e.bits, g_data + r + 64, sizeof e.bits);
        if (e.ncoef < 1 || e.ncoef > kMaxCoef) return false;
        int64_t sum = 0;
        for (int c = 0; c < 3; ++c)
            for (int k = 0; k < e.ncoef; ++k) {
                if (e.bits[c][k] > 56) return false;
                sum += e.bits[c][k];
            }
        if (e.body < 0 || e.body >= EPH_BODY_COUNT || e.nseg < 1 ||
            !(e.length > 0.0) || sum != e.segbits ||
            off + ((uint64_t)e.nseg * e.segbits + 7) / 8 + 8 > (uint64_t)size)
            return false;
        e.data = g_data + off;
        /* Set up the truncated series the residuals were made against. */
        if (e.base) {
            if (e.body == EPH_MOON)
                elp_truncate(e.threshold, tmax);
            else if (e.body <= EPH_NEPTUNE && e.body != EPH_EARTH)
                vsop87_truncate(e.body, e.threshold, tmax);
            else
                return false;
        }
    }
    g_nentries = n;
    return true;
}

}  // namespace

extern "C" int eph_load_compressed(const char* path) {
    eph_unload_compressed();
    FILE* f = std::fopen(path, "rb");
    if (!f) return -1;
    std::fseek(f, 0, SEEK_END);
    const long size = std::ftell(f);
    std::fseek(f, 0, SEEK_SET);
    g_data = static_cast<unsigned char*>(std::malloc(size > 0 ? size : 1));
    const bool ok = g_data && std::fread(g_data, 1, size, f) == (size_t)size;
    std::fclose(f);
    if (!ok || !parse(size)) {
        eph_unload_compressed();
        return -1;
    }
    return 0;
}

extern "C" void eph_unload_compressed(void) {
    std::free(g_data);
    g_data = nullptr;
    g_nentries = 0;
    ++g_generation;
}

extern "C" int eph_compressed_available(void) { return g_data ? 1 : 0; }

bool compressed_residual(int body, double jd, Vec3& p, Vec3& v, Vec3& a,
                         bool& base) {
    if (!g_data) return false;
    const Entry* e = find_entry(body);
    if (!e) return false;
    const double x = (jd - e->start) / e->length;
    if (!(x >= 0.0 && x <= e->nseg)) return false;
    const int i = (x >= e->nseg) ? e->nseg - 1 : (int)x;
    const double u = 2.0 * (x - i) - 1.0;
    const double du = 2.0 / e->length; /* du/dt */

    double out[3][3];
    int64_t pos = i * e->segbits;
    for (int c = 0; c < 3; ++c) {
        double cc[kMaxCoef], d1[kMaxCoef], d2[kMaxCoef];
        const int n = e->ncoef;
        for (int k = 0; k < n; ++k) {
            cc[k] = read_bits(e->data, pos, e->bits[c][k]) * e->quantum;
            pos += e->bits[c][k];
        }
        chebyshev_derivative(cc, n, d1);
        chebyshev_derivative(d1, n, d2);
        out[c][0] = chebyshev(cc, n, u);
        out[c][1] = chebyshev(d1, n, u) * du;
        out[c][2] = chebyshev(d2, n, u) * du * du;
    }
    p = {out[0][0], out[1][0], out[2][0]};
    v = {out[0][1], out[1][1], out[2][1]};
    a = {out[0][2], out[1][2], out[2][2]};
    base = e->base;
    return true;
}

bool compressed_covers(double jd) {
    return g_data && jd >= g_jd0 && jd <= g_jd1;
}

int compressed_generation() { return g_generation; }

}  // namespace eph
