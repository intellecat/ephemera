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

add_library(ephemeris STATIC
    src/time.cpp
    src/frames.cpp
    src/vsop87.cpp
    src/moon.cpp
    src/pluto.cpp
    src/ephemeris.cpp
    data/vsop87d_data.cpp
    data/pluto_data.cpp
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

if(EPHEMERIS_BUILD_TESTS)
    enable_testing()
    add_executable(test_ephemeris tests/test_ephemeris.cpp)
    target_link_libraries(test_ephemeris PRIVATE ephemeris)
    add_test(NAME ephemeris COMMAND test_ephemeris)
endif()
