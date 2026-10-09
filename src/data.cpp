#include "data.h"

#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "../include/ephemeris/ephemeris.h"
#include "internal.h"

namespace eph {

/* Sections are used in place, so the coefficient structs must have the
 * file's layout. */
static_assert(sizeof(vsop87d::Term) == 16 && offsetof(vsop87d::Term, b) == 4 &&
                  offsetof(vsop87d::Term, c) == 8,
              "vsop87d::Term layout");
static_assert(sizeof(elp::MainTerm) == 40 && offsetof(elp::MainTerm, amp) == 32 &&
                  offsetof(elp::MainTerm, packed) == 36,
              "elp::MainTerm layout");
static_assert(sizeof(elp::PertTerm) == 16 && offsetof(elp::PertTerm, amp) == 8 &&
                  offsetof(elp::PertTerm, packed) == 12,
              "elp::PertTerm layout");
static_assert(sizeof(corr::Harmonic) == 16, "corr::Harmonic layout");

namespace {

constexpr int kMaxFiles = 32;

struct File {
    const unsigned char* p;
    uint64_t size;
    bool owned; /* copied by eph_load_data */
    bool has_table;
};

File g_files[kMaxFiles];
int g_nfiles = 0;
Data g_data;
int g_generation = 1;

/* CRC-32 (IEEE), continuing from crc as zlib.crc32(data, crc) does. */
uint32_t crc32(const unsigned char* p, uint64_t n, uint32_t crc = 0) {
    static uint32_t table[256];
    if (!table[1])
        for (uint32_t i = 0; i < 256; ++i) {
            uint32_t c = i;
            for (int k = 0; k < 8; ++k) c = (c & 1) ? 0xEDB88320u ^ (c >> 1) : c >> 1;
            table[i] = c;
        }
    uint32_t c = crc ^ 0xFFFFFFFFu;
    for (uint64_t i = 0; i < n; ++i) c = table[(c ^ p[i]) & 0xFF] ^ (c >> 8);
    return c ^ 0xFFFFFFFFu;
}

template <typename T>
T get(const unsigned char* p, uint64_t off) {
    T v;
    std::memcpy(&v, p + off, sizeof(T));
    return v;
}

/* A bounds-checked view of one section payload. */
struct Reader {
    const unsigned char* p;
    uint64_t size, pos;

    bool take(uint64_t n, const unsigned char** out) {
        if (n > size - pos) return false;
        *out = p + pos;
        pos += n;
        return true;
    }
    template <typename T>
    bool read(T& v) {
        const unsigned char* q;
        if (!take(sizeof(T), &q)) return false;
        std::memcpy(&v, q, sizeof(T));
        return true;
    }
    /* An array of n structs used in place. */
    template <typename T>
    bool array(uint64_t n, const T** out) {
        const unsigned char* q;
        if (n > (size - pos) / sizeof(T) || !take(n * sizeof(T), &q))
            return false;
        *out = reinterpret_cast<const T*>(q);
        return true;
    }
};

bool parse_vsop(Reader r, SeriesSet& s, int planet) {
    uint32_t n[4];
    for (uint32_t& x : n)
        if (!r.read(x)) return false;
    vsop87d::Series(&series)[3][6] = s.vsop_series[planet];
    uint32_t count[3][6];
    for (int c = 0; c < 3; ++c) {
        if (n[c] > 6) return false;
        for (uint32_t a = 0; a < n[c]; ++a) {
            uint32_t pad;
            if (!r.read(series[c][a].poly) || !r.read(count[c][a]) ||
                !r.read(pad))
                return false;
        }
    }
    for (int c = 0; c < 3; ++c)
        for (uint32_t a = 0; a < n[c]; ++a) {
            if (!r.array(count[c][a], &series[c][a].terms)) return false;
            series[c][a].count = (int)count[c][a];
        }
    s.vsop[planet] = {series[0], (int)n[0], series[1], (int)n[1], series[2],
                      (int)n[2]};
    s.vsop_ok[planet] = true;
    return true;
}

bool parse_elp(Reader r, SeriesSet& s) {
    uint32_t n[6];
    for (uint32_t& x : n)
        if (!r.read(x)) return false;
    for (int iv = 0; iv < 3; ++iv) {
        elp::Coord& c = s.elp[iv];
        if (!r.array(n[2 * iv], &c.main) || !r.array(n[2 * iv + 1], &c.pert))
            return false;
        c.nmain = (int)n[2 * iv];
        c.npert = (int)n[2 * iv + 1];
    }
    s.elp_ok = true;
    return true;
}

/* Chebyshev intervals (asteroids, Pluto): u32 nseg, ncoef; f64 t[nseg+1];
 * f32 c[nseg][3][ncoef]. */
template <typename Body>
bool parse_chebyshev(Reader r, Body& b) {
    uint32_t nseg, ncoef;
    if (!r.read(nseg) || !r.read(ncoef) || nseg < 1 || ncoef < 1 ||
        ncoef > 64 || !r.array(nseg + 1ull, &b.t) ||
        !r.array(nseg * 3ull * ncoef, &b.c))
        return false;
    for (uint32_t i = 0; i < nseg; ++i)
        if (!(b.t[i + 1] > b.t[i])) return false;
    b.nseg = (int)nseg;
    b.ncoef = (int)ncoef;
    return true;
}

bool parse_correction(Reader r, Data& d, int body) {
    corr::Series& s = d.corr[body];
    uint32_t n0, nh, nharm, pad;
    if (!r.read(s.tmin) || !r.read(s.tmax) || !r.read(n0) || !r.read(nh) ||
        !r.read(nharm) || !r.read(pad) || n0 < 1 || n0 > 1000 || nh > 100 ||
        nharm > 100 || (nharm > 0 && nh < 1) || !(s.tmax > s.tmin) ||
        !r.array(nharm, &s.harm) ||
        !r.array(3ull * (n0 + 2ull * nharm * nh), &s.coef))
        return false;
    s.n0 = (int)n0;
    s.nh = (int)nh;
    s.nharm = (int)nharm;
    d.corr_ok[body] = true;
    return true;
}

/* Stars: u32 count, u32 string pool size, records, string pool. */
bool parse_stars(Reader r, Data& d) {
    uint32_t n, pool;
    const StarRecord* stars;
    const char* strings;
    if (!r.read(n) || !r.read(pool) || n > 100000 || !r.array(n, &stars) ||
        !r.array(pool, &strings) || pool == 0 || strings[pool - 1] != '\0')
        return false;
    for (uint32_t i = 0; i < n; ++i)
        if (stars[i].name >= pool || stars[i].designation >= pool) return false;
    d.stars = stars;
    d.nstars = (int)n;
    d.star_strings = strings;
    d.star_strings_size = pool;
    return true;
}

/* Apply every section of a file to d (later files override earlier
 * ones).  False if the file is malformed. */
bool apply(const unsigned char* p, uint64_t size, Data& d, bool& has_table) {
    if (size < 16 || std::memcmp(p, "EPHD", 4) != 0 || get<uint32_t>(p, 4) != 2 ||
        get<uint32_t>(p, 12) != crc32(p + 16, size - 16))
        return false;
    const uint32_t n = get<uint32_t>(p, 8);
    if (n > 256 || 16 + 24ull * n > size) return false;
    has_table = false;
    for (uint32_t i = 0; i < n; ++i) {
        const uint64_t h = 16 + 24ull * i;
        const uint32_t type = get<uint32_t>(p, h), id = get<uint32_t>(p, h + 4);
        const uint64_t off = get<uint64_t>(p, h + 8), len = get<uint64_t>(p, h + 16);
        if (off % 8 != 0 || off > size || len > size - off) return false;
        const Reader r = {p + off, len, 0};
        const uint32_t role = id >> 8, sub = id & 0xFF;
        bool ok = true;
        switch (type) {
            case kSectionVsop:
                ok = role <= kBase && sub < 8 && parse_vsop(r, d.set[role], sub);
                if (ok && role == kBase) {
                    d.base_section[sub] = p + off;
                    d.base_size[sub] = len;
                }
                break;
            case kSectionElp:
                ok = role <= kBase && sub == 0 && parse_elp(r, d.set[role]);
                if (ok && role == kBase) {
                    d.base_section[8] = p + off;
                    d.base_size[8] = len;
                }
                break;
            case kSectionPluto:
                ok = d.pluto_ok = parse_chebyshev(r, d.pluto);
                break;
            case kSectionCorrection:
                ok = (id < 8 || id == 10) && parse_correction(r, d, id);
                break;
            case kSectionAsteroid:
                ok = id < kAsteroids &&
                     (d.asteroid_ok[id] = parse_chebyshev(r, d.asteroid[id]));
                break;
            case kSectionStars:
                ok = parse_stars(r, d);
                break;
            case kSectionTable:
                ok = d.ntables < kMaxTables;
                if (ok) d.tables[d.ntables++] = {p + off, len};
                has_table = true;
                break;
            default: /* info and unknown sections are skipped */
                break;
        }
        if (!ok) return false;
    }
    return true;
}

/* Rebuild g_data from the loaded files; false if a residual table cannot
 * be attached. */
bool rebuild() {
    g_data = Data();
    for (int i = 0; i < g_nfiles; ++i) {
        bool t;
        apply(g_files[i].p, g_files[i].size, g_data, t);
    }
    ++g_generation;
    uint32_t crc = 0;
    for (int i = 0; i < 9; ++i)
        if (g_data.base_section[i])
            crc = crc32(g_data.base_section[i], g_data.base_size[i], crc);
    g_data.base_crc = crc;
    return tables_attach();
}

void remove_file(int i) {
    if (g_files[i].owned) std::free(const_cast<unsigned char*>(g_files[i].p));
    for (int k = i + 1; k < g_nfiles; ++k) g_files[k - 1] = g_files[k];
    --g_nfiles;
}

}  // namespace

const Data& data() { return g_data; }

int data_generation() { return g_generation; }

}  // namespace eph

using namespace eph;

namespace {

/* Validate and register a file kept at p (owned: free it on unload). */
int add_file(const unsigned char* p, uint64_t size, bool owned) {
    Data scratch = Data();
    bool has_table;
    if (g_nfiles == kMaxFiles || !apply(p, size, scratch, has_table)) {
        if (owned) std::free(const_cast<unsigned char*>(p));
        return -1;
    }
    g_files[g_nfiles++] = {p, size, owned, has_table};
    if (!rebuild()) {
        remove_file(g_nfiles - 1);
        rebuild();
        return -1;
    }
    return 0;
}

bool little_endian() {
    const uint32_t one = 1;
    unsigned char b;
    std::memcpy(&b, &one, 1);
    return b == 1;
}

}  // namespace

extern "C" int eph_load_data(const void* buffer, size_t size) {
    if (!buffer || !little_endian()) return -1;
    unsigned char* p = static_cast<unsigned char*>(std::malloc(size ? size : 1));
    if (!p) return -1;
    std::memcpy(p, buffer, size);
    return add_file(p, size, true);
}

extern "C" int eph_load_static(const void* buffer, size_t size) {
    /* Used in place: the sections hold doubles. */
    if (!buffer || !little_endian() || reinterpret_cast<uintptr_t>(buffer) % 8)
        return -1;
    return add_file(static_cast<const unsigned char*>(buffer), size, false);
}

extern "C" int eph_load_embedded(void) {
    const unsigned char* p;
    size_t n;
    if (!embedded_data(p, n)) return -1;
    return eph_load_static(p, n);
}

extern "C" int eph_load_file(const char* path) {
#ifdef EPHEMERIS_NO_FILE_IO
    (void)path;
    return -1;
#else
    FILE* f = std::fopen(path, "rb");
    if (!f) return -1;
    int rc = -1;
    if (std::fseek(f, 0, SEEK_END) == 0) {
        const long size = std::ftell(f);
        std::fseek(f, 0, SEEK_SET);
        void* buf = size > 0 ? std::malloc(size) : nullptr;
        if (buf && std::fread(buf, 1, size, f) == (size_t)size)
            rc = eph_load_data(buf, (size_t)size);
        std::free(buf);
    }
    std::fclose(f);
    return rc;
#endif
}

extern "C" void eph_unload_data(void) {
    while (g_nfiles > 0) remove_file(g_nfiles - 1);
    rebuild();
}

/* Compatibility: the table loader of earlier versions. */
extern "C" int eph_load_compressed(const char* path) {
    return eph_load_file(path);
}

extern "C" void eph_unload_compressed(void) {
    for (int i = g_nfiles - 1; i >= 0; --i)
        if (g_files[i].has_table) remove_file(i);
    rebuild();
}

extern "C" int eph_compressed_available(void) {
    return g_data.ntables > 0 ? 1 : 0;
}
