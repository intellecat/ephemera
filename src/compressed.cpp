#include "internal.h"

#include <cstdint>
#include <cstring>

#include "../include/ephemeris/ephemeris.h"
#include "data.h"

namespace eph {

/*
 * Residual table: a section of data/ephem.bin (tools/pack_ephemeris.py).
 *
 * For each body it stores DE440 minus a base series -- truncated VSOP87D /
 * ELP2000-82B series stored in the same file (role kBase), or nothing -- as Chebyshev polynomials over fixed-length
 * segments.  The residual is small and smooth, so its coefficients are
 * quantized with one absolute step per body and packed with a fixed bit
 * width per coefficient.  Positions are in the IAU-1976 ecliptic of date,
 * AU; heliocentric, geocentric for the Moon.  The Earth's entry holds the
 * Earth-Moon barycentre, which is smooth (the Earth is derived from it and
 * the Moon).
 *
 * Layout (little-endian, offsets relative to the section):
 *   0  "EPHP", u32 version (3), u32 body count, f64 jd0, f64 jd1
 *  28  per body, 160 bytes:
 *        u32 body, u32 base (0 none, 1 series), u32 nseg, u32 ncoef,
 *        f64 threshold the base was truncated at (rad, informative),
 *        f64 start, f64 segment length (days),
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

const unsigned char* g_table = nullptr;
Entry g_entries[16];
int g_nentries = 0;
double g_jd0 = 0.0, g_jd1 = 0.0;

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

bool parse(const unsigned char* t, uint64_t size) {
    if (size < kHeader || std::memcmp(t, "EPHP", 4) != 0) return false;
    const uint32_t version = load_at<uint32_t>(t, 4);
    const uint32_t n = load_at<uint32_t>(t, 8);
    if (version != 3 || n > 16 || (uint64_t)(kHeader + n * kRecord) > size)
        return false;
    g_jd0 = load_at<double>(t, 12);
    g_jd1 = load_at<double>(t, 20);
    const SeriesSet& base = data().set[kBase];
    for (uint32_t i = 0; i < n; ++i) {
        const long r = kHeader + i * kRecord;
        Entry& e = g_entries[i];
        e.body = load_at<uint32_t>(t, r);
        e.base = load_at<uint32_t>(t, r + 4) != 0;
        e.nseg = load_at<uint32_t>(t, r + 8);
        e.ncoef = load_at<uint32_t>(t, r + 12);
        e.threshold = load_at<double>(t, r + 16);
        e.start = load_at<double>(t, r + 24);
        e.length = load_at<double>(t, r + 32);
        e.quantum = load_at<double>(t, r + 40);
        const uint64_t off = load_at<uint64_t>(t, r + 48);
        e.segbits = load_at<uint32_t>(t, r + 56);
        std::memcpy(e.bits, t + r + 64, sizeof e.bits);
        if (e.ncoef < 1 || e.ncoef > kMaxCoef) return false;
        int64_t sum = 0;
        for (int c = 0; c < 3; ++c)
            for (int k = 0; k < e.ncoef; ++k) {
                if (e.bits[c][k] > 56) return false;
                sum += e.bits[c][k];
            }
        if (e.body < 0 || e.body >= EPH_BODY_COUNT || e.nseg < 1 ||
            !(e.length > 0.0) || sum != e.segbits ||
            off + ((uint64_t)e.nseg * e.segbits + 7) / 8 + 8 > size)
            return false;
        e.data = t + off;
        /* The base series must be in the loaded data. */
        if (e.base && !(e.body == EPH_MOON ? base.elp_ok
                        : e.body <= EPH_NEPTUNE && e.body != EPH_EARTH
                            ? base.vsop_ok[e.body]
                            : false))
            return false;
    }
    g_nentries = n;
    return true;
}

}  // namespace

bool table_attach(const unsigned char* p, unsigned long size) {
    g_table = nullptr;
    g_nentries = 0;
    if (!p) return true;
    if (!parse(p, size)) {
        g_nentries = 0;
        return false;
    }
    g_table = p;
    return true;
}

bool compressed_residual(int body, double jd, Vec3& p, Vec3& v, Vec3& a,
                         bool& base) {
    if (!g_table) return false;
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
    return g_table && jd >= g_jd0 && jd <= g_jd1;
}

}  // namespace eph
