#include "internal.h"

#include <cctype>
#include <cstdlib>
#include <cstring>

#include "../include/ephemeris/ephemeris.h"
#include "data.h"

/* Minor planets from the loaded minor-planet sections (format in data.h,
 * made by tools/make_minor_planets.py). */

namespace eph {

namespace {

/* The record of global index i, with its section. */
const MinorRecord* record(int i, const unsigned char** base, const char** pool) {
    const Data& d = data();
    if (i < 0) return nullptr;
    for (int s = 0; s < d.nminor_sets; ++s) {
        if (i < d.minor[s].count) {
            if (base) *base = d.minor[s].base;
            if (pool) *pool = d.minor[s].pool;
            return &d.minor[s].rec[i];
        }
        i -= d.minor[s].count;
    }
    return nullptr;
}

}  // namespace

bool is_minor(int body) {
    return body >= EPH_MINOR_BASE && body - EPH_MINOR_BASE < data().nminor;
}

bool minor_heliocentric(int body, double jd, Vec3& v) {
    const unsigned char* base;
    const MinorRecord* m = record(body - EPH_MINOR_BASE, &base, nullptr);
    if (!m) return false;
    const double* t = reinterpret_cast<const double*>(base + m->data_off);
    const double* c0 = t + m->nseg + 1;
    const float* c = reinterpret_cast<const float*>(c0 + 3 * m->nseg);
    if (!(jd >= t[0] && jd <= t[m->nseg])) return false;
    int lo = 0, hi = (int)m->nseg - 1;
    while (lo < hi) {
        const int mid = (lo + hi + 1) / 2;
        if (jd < t[mid])
            hi = mid - 1;
        else
            lo = mid;
    }
    const double u = (2.0 * jd - (t[lo] + t[lo + 1])) / (t[lo + 1] - t[lo]);
    const int n = (int)m->ncoef - 1; /* coefficients 1 .. ncoef - 1 */
    double out[3];
    for (int k = 0; k < 3; ++k) {
        const float* ck = c + (lo * 3 + k) * n;
        double b1 = 0.0, b2 = 0.0; /* Clenshaw over T_1 .. T_n */
        for (int j = n; j >= 1; --j) {
            const double b0 = 2.0 * u * b1 - b2 + ck[j - 1];
            b2 = b1;
            b1 = b0;
        }
        out[k] = c0[lo * 3 + k] + u * b1 - b2;
    }
    v = {out[0], out[1], out[2]};
    return true;
}

bool minor_info(int body, double& h, double& g, double& radius, double& q) {
    const MinorRecord* m = record(body - EPH_MINOR_BASE, nullptr, nullptr);
    if (!m) return false;
    h = m->h;
    g = m->g;
    radius = m->radius;
    q = m->q;
    return true;
}

namespace {

bool same(const char* a, const char* b) {
    while (*a && *b) {
        if (std::tolower((unsigned char)*a) != std::tolower((unsigned char)*b)) return false;
        ++a;
        ++b;
    }
    return *a == *b;
}

}  // namespace

}  // namespace eph

using namespace eph;

extern "C" {

int eph_minor_count(void) { return data().nminor; }

int eph_minor_number(int index) {
    const MinorRecord* m = record(index, nullptr, nullptr);
    return m ? (int)m->number : 0;
}

const char* eph_minor_name(int index) {
    const char* pool;
    const MinorRecord* m = record(index, nullptr, &pool);
    return m ? pool + m->name : nullptr;
}

int eph_minor_find(const char* name) {
    if (!name) return -1;
    while (*name == ' ') ++name;
    /* "433", "433 Eros", "(433) Eros" or "Eros". */
    const char* p = name;
    if (*p == '(') ++p;
    char* end;
    const long number = std::strtol(p, &end, 10);
    const bool numbered = end != p;
    if (numbered && *end == ')') ++end;
    while (numbered && *end == ' ') ++end;
    const char* text = numbered ? end : name;
    const int n = data().nminor;
    for (int i = 0; i < n; ++i) {
        const char* pool;
        const MinorRecord* m = record(i, nullptr, &pool);
        if (numbered && (long)m->number != number) continue;
        if (*text && !same(text, pool + m->name)) continue;
        return i;
    }
    return -1;
}

}  // extern "C"
