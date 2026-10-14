# The Python package's native library: Ephemera as a shared library with the
# 1" series and Chiron (eph-1as-chiron.bin) built in, loaded with ctypes.
# Built from the sources staged in csrc/ by `make python`.
cmake_minimum_required(VERSION 3.16)
project(pyephemera C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

set(CSRC ${CMAKE_CURRENT_SOURCE_DIR}/csrc)
set(EMBED_SRC ${CMAKE_CURRENT_BINARY_DIR}/embedded_data.cpp)
add_custom_command(OUTPUT ${EMBED_SRC}
    COMMAND ${CMAKE_COMMAND} -DIN=${CSRC}/eph-1as-chiron.bin -DOUT=${EMBED_SRC}
            -P ${CSRC}/tools/embed.cmake
    DEPENDS ${CSRC}/eph-1as-chiron.bin ${CSRC}/tools/embed.cmake)

file(GLOB SOURCES ${CSRC}/src/*.cpp)
list(FILTER SOURCES EXCLUDE REGEX "embedded_none\\.cpp$")

add_library(ephemera_native SHARED ${SOURCES} ${EMBED_SRC} ${CSRC}/data/delta_t_data.cpp)
target_include_directories(ephemera_native PRIVATE ${CSRC}/include ${CSRC}/src ${CSRC}/data)
set_target_properties(ephemera_native PROPERTIES WINDOWS_EXPORT_ALL_SYMBOLS ON)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(ephemera_native PRIVATE -O2)
endif()
if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
    # No dependency on the system's C++ runtime.
    target_link_options(ephemera_native PRIVATE -static-libstdc++ -static-libgcc)
endif()

install(TARGETS ephemera_native LIBRARY DESTINATION ephemera RUNTIME DESTINATION ephemera)
