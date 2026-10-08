#include "internal.h"

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "../include/ephemeris/ephemeris.h"

namespace eph {

namespace {

struct Entry {
    int body;
    int nseg;
    int ncoef;
    const double* tstart;
    const float* scale;
    const unsigned char* coef;  /* 24-bit little-endian */
};

unsigned char* g_data = nullptr;
long g_size = 0;
Entry g_entries[16];
int g_nentries = 0;
double g_jd0 = 0.0, g_jd1 = 0.0;      /* advertised range */
double g_data0 = 0.0, g_data1 = 0.0;  /* actual stored range (with margin) */
int g_cbytes = 3;                     /* coefficient width: 3 (int24) or 4 (int32) */

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

const Entry* find_entry(int body) {
    for (int i = 0; i < g_nentries; ++i)
        if (g_entries[i].body == body) return &g_entries[i];
    return nullptr;
}

}  // namespace

extern "C" int eph_load_compressed(const char* path) {
    eph_unload_compressed();
    FILE* f = std::fopen(path, "rb");
    if (!f) return -1;
    std::fseek(f, 0, SEEK_END);
    g_size = std::ftell(f);
    std::fseek(f, 0, SEEK_SET);
    g_data = static_cast<unsigned char*>(std::malloc(g_size));
    if (!g_data || std::fread(g_data, 1, g_size, f) != (size_t)g_size) {
        std::fclose(f);
        eph_unload_compressed();
        return -1;
    }
    std::fclose(f);

    if (std::memcmp(g_data, "EPHP", 4) != 0) {
        eph_unload_compressed();
        return -1;
    }
    const uint32_t version = load_at<uint32_t>(g_data, 4);
    const uint32_t nindex = load_at<uint32_t>(g_data, 28);
    if (version < 1 || version > 2 || nindex > 16) {
        eph_unload_compressed();
        return -1;
    }
    g_cbytes = (version == 1) ? 4 : 3;
    g_jd0 = load_at<double>(g_data, 12);
    g_jd1 = load_at<double>(g_data, 20);
    const uint32_t idxlen = load_at<uint32_t>(g_data, 32);
    long off = 36;
    for (uint32_t i = 0; i < nindex; ++i) {
        Entry& e = g_entries[i];
        e.body = load_at<uint32_t>(g_data, off + 0);
        e.nseg = load_at<uint32_t>(g_data, off + 4);
        e.ncoef = load_at<uint32_t>(g_data, off + 8);
        const long to = load_at<uint64_t>(g_data, off + 16);
        const long so = load_at<uint64_t>(g_data, off + 24);
        const long co = load_at<uint64_t>(g_data, off + 32);
        e.tstart = reinterpret_cast<const double*>(g_data + to);
        e.scale = reinterpret_cast<const float*>(g_data + so);
        e.coef = g_data + co;
        off += 40;
    }
    (void)idxlen;
    g_nentries = nindex;
    g_data0 = g_jd0;
    g_data1 = g_jd1;
    for (int i = 0; i < g_nentries; ++i) {
        const Entry& e = g_entries[i];
        if (e.tstart[0] < g_data0) g_data0 = e.tstart[0];
        if (e.tstart[e.nseg] > g_data1) g_data1 = e.tstart[e.nseg];
    }
    return 0;
}

extern "C" void eph_unload_compressed(void) {
    std::free(g_data);
    g_data = nullptr;
    g_size = 0;
    g_nentries = 0;
}

extern "C" int eph_compressed_available(void) { return g_data ? 1 : 0; }

namespace {

/* Evaluate the stored position (heliocentric, or geocentric for the Moon)
 * in the mean ecliptic J2000 frame. */
bool eval_entry(const Entry* e, double jd, Vec3& v) {
    int lo = 0, hi = e->nseg - 1;
    while (lo < hi) {
        const int mid = (lo + hi + 1) / 2;
        if (jd < e->tstart[mid])
            hi = mid - 1;
        else
            lo = mid;
    }
    const int i = lo;
    const double t0 = e->tstart[i], t1 = e->tstart[i + 1];
    const double u = (2.0 * jd - (t0 + t1)) / (t1 - t0);

    double out[3];
    for (int c = 0; c < 3; ++c) {
        const unsigned char* q =
            e->coef + (i * 3 + c) * e->ncoef * g_cbytes;
        double cc[24];
        const int n = e->ncoef;
        for (int k = 0; k < n; ++k) {
            int32_t v;
            if (g_cbytes == 4) {
                v = (int32_t)(q[4 * k] | (q[4 * k + 1] << 8) |
                              (q[4 * k + 2] << 16) | (uint32_t)(q[4 * k + 3] << 24));
            } else {
                v = q[3 * k] | (q[3 * k + 1] << 8) | (q[3 * k + 2] << 16);
                if (v & 0x800000) v |= ~0xFFFFFF; /* sign-extend */
            }
            cc[k] = v * (double)e->scale[i * 3 + c];
        }
        out[c] = chebyshev(cc, n, u);
    }
    v = {out[0], out[1], out[2]};
    return true;
}

}  // namespace

bool compressed_raw(int body, double jd, Vec3& v) {
    if (!g_data || jd < g_data0 || jd > g_data1) return false;
    const Entry* e = find_entry(body);
    if (!e) return false;
    return eval_entry(e, jd, v);
}

bool compressed_geocentric(int body, double jd, bool light_time, Vec3& out,
                           double& dist) {
    if (!g_data || jd < g_jd0 || jd > g_jd1) return false;

    if (body == EPH_MOON) {
        Vec3 p;
        if (!compressed_raw(EPH_MOON, jd, p)) return false;
        if (light_time) {
            for (int i = 0; i < 2; ++i) {
                double te = jd - length(p) / AU_PER_LIGHT_DAY;
                if (te < g_data0) te = g_data0;
                if (te > g_data1) te = g_data1;
                Vec3 q;
                if (!compressed_raw(EPH_MOON, te, q)) break;
                p = q;
            }
        }
        out = p;
        dist = length(p);
        return true;
    }

    if (body == EPH_SUN) {
        Vec3 e;
        if (!compressed_raw(EPH_EARTH, jd, e)) return false;
        out = Vec3{-e.x, -e.y, -e.z};
        dist = length(out);
        return true;
    }

    Vec3 e;
    if (!compressed_raw(EPH_EARTH, jd, e)) return false;
    Vec3 p;
    if (!compressed_raw(body, jd, p)) return false;
    out = p - e;
    dist = length(out);
    if (light_time) {
        for (int i = 0; i < 3; ++i) {
            double te = jd - dist / AU_PER_LIGHT_DAY;
            if (te < g_data0) te = g_data0;
            if (te > g_data1) te = g_data1;
            Vec3 q;
            if (!compressed_raw(body, te, q)) break;
            out = q - e;
            dist = length(out);
        }
    }
    return true;
}

}  // namespace eph
