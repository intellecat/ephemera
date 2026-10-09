# Compile a data file into the library:
#   cmake -DIN=data/eph-1as.bin -DOUT=embedded_data.cpp -P tools/embed.cmake
# writes a source that defines eph::embedded_data() and loads the file when
# the program starts (used in place, without a copy).
file(READ "${IN}" hex HEX)
string(LENGTH "${hex}" n)
math(EXPR size "${n} / 2")
string(REGEX REPLACE "([0-9a-f][0-9a-f])" "0x\\1," bytes "${hex}")
string(REGEX REPLACE "((0x..,){20})" "\\1\n    " bytes "${bytes}")
get_filename_component(name "${IN}" NAME)
file(WRITE "${OUT}" "/* ${name}, embedded by tools/embed.cmake -- do not edit. */
#include \"data.h\"

#include \"ephemeris/ephemeris.h\"

namespace eph {

namespace {

alignas(8) const unsigned char kData[${size}] = {
    ${bytes}
};

/* Load it at start-up, so the library works without any call. */
struct AutoLoad {
    AutoLoad() { eph_load_embedded(); }
} g_autoload;

}  // namespace

bool embedded_data(const unsigned char*& p, size_t& size) {
    p = kData;
    size = sizeof kData;
    return true;
}

}  // namespace eph
")
