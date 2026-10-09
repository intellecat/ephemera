/* Load the data files named on a test's command line. */
#pragma once
#include <cstdio>

#include "ephemeris/ephemeris.h"

static bool load_data_args(int argc, char** argv) {
    for (int i = 1; i < argc; ++i)
        if (eph_load_file(argv[i]) != 0) {
            std::printf("FAIL: cannot load %s\n", argv[i]);
            return false;
        }
    return true;
}
