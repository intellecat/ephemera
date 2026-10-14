# The Android library's native code: Ephemera from the repository's src/
# with the 1" series and Chiron built in, and the JNI glue.
cmake_minimum_required(VERSION 3.18)
project(ephemera_android C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
get_filename_component(ROOT ${CMAKE_CURRENT_SOURCE_DIR}/../../../../.. ABSOLUTE)

set(EMBED_SRC ${CMAKE_CURRENT_BINARY_DIR}/embedded_data.cpp)
add_custom_command(OUTPUT ${EMBED_SRC}
    COMMAND ${CMAKE_COMMAND} -DIN=${ROOT}/data/eph-1as-chiron.bin -DOUT=${EMBED_SRC}
            -P ${ROOT}/tools/embed.cmake
    DEPENDS ${ROOT}/data/eph-1as-chiron.bin ${ROOT}/tools/embed.cmake)

file(GLOB SOURCES ${ROOT}/src/*.cpp)
list(FILTER SOURCES EXCLUDE REGEX "embedded_none\\.cpp$")

add_library(ephemera SHARED ${SOURCES} ${EMBED_SRC} ${ROOT}/data/delta_t_data.cpp ephemera_jni.c)
target_include_directories(ephemera PRIVATE ${ROOT}/include ${ROOT}/src ${ROOT}/data)
target_compile_options(ephemera PRIVATE -O2)
# 16 KB pages (Android 15 and later).
target_link_options(ephemera PRIVATE -Wl,-z,max-page-size=16384)
target_link_libraries(ephemera m)
