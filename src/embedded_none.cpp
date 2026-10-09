#include "data.h"

namespace eph {

/* Built without EPHEMERIS_EMBED_DATA: no data file in the library. */
bool embedded_data(const unsigned char*& p, size_t& size) {
    p = nullptr;
    size = 0;
    return false;
}

}  // namespace eph
