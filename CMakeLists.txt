cmake_minimum_required(VERSION 3.16)
project(ephemeris C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_C_STANDARD 99)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

option(EPHEMERIS_BUILD_TESTS "Build the test suite" ON)
option(EPHEMERIS_BUILD_EXAMPLES "Build the example programs" ON)
option(EPHEMERIS_BUILD_TOOLS "Build tools/pack_series (data file generator)" OFF)
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
    ${EMBED_SRC}
    data/delta_t_data.cpp
)

target_include_directories(ephemeris
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src ${CMAKE_CURRENT_SOURCE_DIR}/data
)

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
        data/corr_data.cpp data/asteroid_data.cpp)
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
endif()
