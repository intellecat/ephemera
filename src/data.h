#ifndef EPH_DATA_H
#define EPH_DATA_H

/* Series and tables loaded from data files (src/data.cpp).
 *
 * The coefficient structures are those of the data/ headers; their arrays
 * point into the loaded file buffers.  Each set is null/zero when no
 * loaded file provides it. */

#include "../data/asteroid_data.h"
#include "../data/corr_data.h"
#include "../data/elp_data.h"
#include "../data/pluto_data.h"
#include "../data/vsop87d_data.h"

namespace eph {

/* Container format (little-endian):
 *   0  "EPHD", u32 version (1), u32 section count,
 *      u32 CRC-32 (zlib) of everything from byte 16 on
 *  16  sections: u32 type, u32 id, u64 offset, u64 size  (24 bytes each)
 * Section payloads start at 8-byte aligned offsets. */
enum SectionType {
    kSectionVsop = 1,       /* id: planet | role << 8  */
    kSectionElp = 2,        /* id: role << 8           */
    kSectionPluto = 3,      /* Pluto fit               */
    kSectionCorrection = 4, /* id: planet, or 10 (Moon) */
    kSectionAsteroid = 5,   /* id: 0..3                */
    kSectionTable = 6,      /* DE440 residual table    */
    kSectionInfo = 7,       /* free text               */
};

/* Series role: the theory used for analytical positions, or the truncated
 * base under a residual table. */
enum Role { kTheory = 0, kBase = 1 };

struct SeriesSet {
    bool vsop_ok[8];
    vsop87d::Body vsop[8];
    vsop87d::Series vsop_series[8][3][6];
    bool elp_ok;
    elp::Coord elp[3];
};

struct Data {
    SeriesSet set[2]; /* indexed by Role */
    bool pluto_ok;
    pluto_fit::Model pluto[3]; /* lon, lat, rad */
    bool corr_ok[11];
    corr::Series corr[11];
    bool asteroid_ok[4];
    ast::Body asteroid[4];
    const unsigned char* table; /* residual table section, or null */
    unsigned long table_size;
};

/* The current data (rebuilt on every load/unload). */
const Data& data();

/* Changes on every load/unload. */
int data_generation();

/* Residual table (src/compressed.cpp): parse the table section, or forget
 * it (p == nullptr).  Returns false if it is malformed or its base series
 * are missing. */
bool table_attach(const unsigned char* p, unsigned long size);

}  // namespace eph

#endif
