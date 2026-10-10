cmake_minimum_required(VERSION 3.16)
project(ephemera C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_C_STANDARD 99)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

option(EPHEMERIS_BUILD_TESTS "Build the test suite" ON)
option(EPHEMERIS_BUILD_EXAMPLES "Build the example programs" ON)
option(EPHEMERIS_BUILD_TOOLS "Build tools/pack_series (data file generator)" OFF)
option(EPHEMERIS_ECLIPSES "Include solar and lunar eclipses" ON)
option(EPHEMERIS_STARS "Include fixed stars (and the star-based ayanamsas)" ON)
set(EPHEMERIS_EMBED_DATA "" CACHE FILEPATH
    "Data file to compile into the library (e.g. data/eph-1as.bin)")

if(EPHEMERIS_EMBED_DATA)
    get_filename_component(EMBED_IN "${EPHEMERIS_EMBED_DATA}" ABSOLUTE)
    set(EMBED_SRC ${CMAKE_CURRENT_BINARY_DIR}/embedded_data.cpp)
    add_custom_command(OUTPUT ${EMBED_SRC}
        COMMAND ${CMAKE_COMMAND} -DIN=${EMBED_IN} -DOUT=${EMBED_SRC}
                -P ${CMAKE_CURRENT_SOURCE_DIR}/tools/embed.cmake
        DEPENDS ${EMBED_IN} ${CMAKE_CURRENT_SOURCE_DIR}/tools/embed.cmake
        COMMENT "Embedding ${EPHEMERIS_EMBED_DATA}")
else()
    set(EMBED_SRC src/embedded_none.cpp)
endif()

add_library(ephemeris STATIC
    src/time.cpp
    src/frames.cpp
    src/vsop87.cpp
    src/moon.cpp
    src/elp.cpp
    src/correction.cpp
    src/compressed.cpp
    src/houses.cpp
    src/asteroids.cpp
    src/pluto.cpp
    src/ephemeris.cpp
    src/data.cpp
    src/util.cpp
    src/ayanamsa.cpp
    src/search.cpp
    src/phenomena.cpp
    src/stars.cpp
    src/eclipses.cpp
    ${EMBED_SRC}
    data/delta_t_data.cpp
)

target_include_directories(ephemeris
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src ${CMAKE_CURRENT_SOURCE_DIR}/data
)

# Left-out features keep their functions as stubs returning -1.
if(NOT EPHEMERIS_ECLIPSES)
    target_compile_definitions(ephemeris PRIVATE EPHEMERIS_NO_ECLIPSES)
endif()
if(NOT EPHEMERIS_STARS)
    target_compile_definitions(ephemeris PRIVATE EPHEMERIS_NO_STARS)
endif()

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(ephemeris PRIVATE -Wall -Wextra -O2)
    # The generated VSOP87 tables are large; keep optimisation but allow time.
endif()

if(EPHEMERIS_BUILD_EXAMPLES)
    add_executable(ephcliapp examples/positions.c)
    target_link_libraries(ephcliapp PRIVATE ephemeris)
endif()

if(EPHEMERIS_BUILD_TOOLS)
    add_executable(pack_series tools/pack_series.cpp
        data/vsop87d_data.cpp data/elp_data.cpp data/pluto_data.cpp
        data/corr_data.cpp data/asteroid_data.cpp data/stars_data.cpp)
    target_include_directories(pack_series PRIVATE src data)
endif()

if(EPHEMERIS_BUILD_TESTS)
    enable_testing()
    set(DATA ${CMAKE_CURRENT_SOURCE_DIR}/data)
    add_executable(test_ephemeris tests/test_ephemeris.cpp)
    target_link_libraries(test_ephemeris PRIVATE ephemeris)
    add_test(NAME ephemeris COMMAND test_ephemeris ${DATA}/eph-full.bin)

    add_executable(test_compressed tests/test_compressed.cpp)
    target_link_libraries(test_compressed PRIVATE ephemeris)
    add_test(NAME compressed COMMAND test_compressed
             ${DATA}/de440/de440-1800.bin ${DATA}/de440/de440-1900.bin
             ${DATA}/de440/de440-2000.bin ${DATA}/de440/de440-2100.bin)

    add_executable(test_topocentric tests/test_topocentric.cpp)
    target_link_libraries(test_topocentric PRIVATE ephemeris)
    add_test(NAME topocentric COMMAND test_topocentric ${DATA}/eph-full.bin)

    add_executable(test_houses tests/test_houses.cpp)
    target_link_libraries(test_houses PRIVATE ephemeris)
    add_test(NAME houses COMMAND test_houses)

    add_executable(test_asteroids tests/test_asteroids.cpp)
    target_link_libraries(test_asteroids PRIVATE ephemeris)
    add_test(NAME asteroids COMMAND test_asteroids ${DATA}/eph-full.bin
             ${DATA}/eph-asteroids.bin)

    add_executable(test_data tests/test_data.cpp)
    target_link_libraries(test_data PRIVATE ephemeris)
    add_test(NAME data COMMAND test_data ${DATA})

    add_executable(test_lunar_points tests/test_lunar_points.cpp)
    target_link_libraries(test_lunar_points PRIVATE ephemeris)
    add_test(NAME lunar_points COMMAND test_lunar_points ${DATA}/eph-full.bin)

    add_executable(test_time tests/test_time.cpp)
    target_link_libraries(test_time PRIVATE ephemeris)
    add_test(NAME time COMMAND test_time ${DATA}/eph-full.bin)

    add_executable(test_search tests/test_search.cpp)
    target_link_libraries(test_search PRIVATE ephemeris)
    add_test(NAME search COMMAND test_search ${DATA}/eph-full.bin)

    add_executable(test_phenomena tests/test_phenomena.cpp)
    target_link_libraries(test_phenomena PRIVATE ephemeris)
    add_test(NAME phenomena COMMAND test_phenomena ${DATA}/eph-full.bin
             ${DATA}/eph-asteroids.bin)

    if(EPHEMERIS_STARS)
        add_executable(test_stars tests/test_stars.cpp)
        target_link_libraries(test_stars PRIVATE ephemeris)
        add_test(NAME stars COMMAND test_stars ${DATA}/eph-full.bin ${DATA}/eph-stars.bin)
    endif()

    if(EPHEMERIS_ECLIPSES)
        add_executable(test_eclipses tests/test_eclipses.cpp)
        target_link_libraries(test_eclipses PRIVATE ephemeris)
        add_test(NAME eclipses COMMAND test_eclipses ${DATA}/eph-full.bin)
    endif()
endif()
