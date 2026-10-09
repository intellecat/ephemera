# Convenience Makefile.  CMake is the primary build system; see README.
CXX      ?= g++
CC       ?= gcc
CXXFLAGS ?= -O2 -std=c++17 -Wall -Wextra
CFLAGS   ?= -O2 -std=c99 -Wall -Wextra
INCLUDES  = -Iinclude -Isrc -Idata
BUILD     = build

CORE_SRC = src/time.cpp src/frames.cpp src/vsop87.cpp src/moon.cpp \
           src/elp.cpp src/pluto.cpp src/ephemeris.cpp src/correction.cpp \
           src/compressed.cpp src/houses.cpp src/asteroids.cpp src/data.cpp \
           data/delta_t_data.cpp
LIB_SRC = $(CORE_SRC)

# make EMBED=data/eph-1as.bin BUILD=build-1as compiles that data file into
# the library (loaded at start-up; cmake generates the source).  Use a
# separate BUILD directory per EMBED setting.
EMBED ?=
ifeq ($(EMBED),)
LIB_SRC += src/embedded_none.cpp
else
LIB_SRC += $(BUILD)/embedded_data.cpp
endif
LIB_OBJ = $(patsubst %.cpp,$(BUILD)/%.o,$(patsubst $(BUILD)/%,%,$(LIB_SRC)))
LIB     = $(BUILD)/libephemeris.a

.PHONY: all test tools data wasm wasm-test clean

# The DE440 tables covering the test reference (1850-2140).
DE440_1800_2200 = $(foreach y,1800 1900 2000 2100,data/de440/de440-$(y).bin)

TESTS = $(BUILD)/test_ephemeris $(BUILD)/test_compressed $(BUILD)/test_topocentric \
        $(BUILD)/test_houses $(BUILD)/test_asteroids $(BUILD)/test_data \
        $(BUILD)/test_lunar_points

all: $(LIB) $(BUILD)/ephcliapp $(TESTS)

# -MMD -MP: rebuild objects when a header they include changes.
$(BUILD)/%.o: %.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -MMD -MP -c $< -o $@

$(BUILD)/embedded_data.cpp: $(EMBED) tools/embed.cmake
	@mkdir -p $(BUILD)
	cmake -DIN=$(EMBED) -DOUT=$@ -P tools/embed.cmake

$(BUILD)/embedded_data.o: $(BUILD)/embedded_data.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -MMD -MP -c $< -o $@

$(LIB): $(LIB_OBJ)
	rm -f $@
	ar rcs $@ $^

$(BUILD)/ephcliapp: examples/positions.c $(LIB)
	$(CC) $(CFLAGS) -Iinclude $< $(LIB) -lm -o $@

$(BUILD)/test_ephemeris: tests/test_ephemeris.cpp tests/reference_data.h $(LIB)
	$(CXX) $(CXXFLAGS) -Iinclude -Itests $< $(LIB) -lm -o $@

test: $(TESTS)
	$(BUILD)/test_ephemeris data/eph-full.bin
	$(BUILD)/test_compressed $(DE440_1800_2200)
	$(BUILD)/test_topocentric data/eph-full.bin
	$(BUILD)/test_houses
	$(BUILD)/test_asteroids data/eph-full.bin data/eph-asteroids.bin
	$(BUILD)/test_data data
	$(BUILD)/test_lunar_points data/eph-full.bin

$(BUILD)/test_compressed: tests/test_compressed.cpp $(LIB)
	$(CXX) $(CXXFLAGS) -Iinclude -Itests $< $(LIB) -lm -o $@

$(BUILD)/test_topocentric: tests/test_topocentric.cpp tests/topo_ref.h $(LIB)
	$(CXX) $(CXXFLAGS) -Iinclude -Itests $< $(LIB) -lm -o $@

$(BUILD)/test_houses: tests/test_houses.cpp tests/houses_ref.h $(LIB)
	$(CXX) $(CXXFLAGS) -Iinclude -Itests $< $(LIB) -lm -o $@

$(BUILD)/test_asteroids: tests/test_asteroids.cpp tests/asteroid_ref.h $(LIB)
	$(CXX) $(CXXFLAGS) -Iinclude -Itests $< $(LIB) -lm -o $@

$(BUILD)/test_data: tests/test_data.cpp tests/reference_data.h $(LIB)
	$(CXX) $(CXXFLAGS) -Iinclude -Itests $< $(LIB) -lm -o $@

$(BUILD)/test_lunar_points: tests/test_lunar_points.cpp tests/lunar_points_ref.h $(LIB)
	$(CXX) $(CXXFLAGS) -Iinclude -Itests $< $(LIB) -lm -o $@

# Development tools and the series data files (see README).
SERIES_SRC = data/vsop87d_data.cpp data/elp_data.cpp data/pluto_data.cpp \
             data/corr_data.cpp data/asteroid_data.cpp

tools: $(BUILD)/pack_series

$(BUILD)/pack_series: tools/pack_series.cpp $(SERIES_SRC)
	@mkdir -p $(BUILD)
	$(CXX) $(CXXFLAGS) $(INCLUDES) $^ -o $@

# Truncation thresholds (radians) per planet, Mercury..Neptune, and for the
# Moon, chosen so each level stays within its accuracy over 1550-2650.
LEVEL_01AS = 3e-8,1e-8,3e-9,1e-9,1e-8,3e-9,1e-8,1e-8 1e-8
LEVEL_1AS  = 3e-7,1e-7,3e-8,3e-8,1e-7,1e-7,1e-7,3e-7 5e-8

# Each series level with and without Chiron, Chiron alone, and all five
# minor bodies.
data: $(BUILD)/pack_series
	$(BUILD)/pack_series series data/eph-full.bin 0 0
	$(BUILD)/pack_series series data/eph-01as.bin $(LEVEL_01AS)
	$(BUILD)/pack_series series data/eph-1as.bin $(LEVEL_1AS)
	$(BUILD)/pack_series series+chiron data/eph-full-chiron.bin 0 0
	$(BUILD)/pack_series series+chiron data/eph-01as-chiron.bin $(LEVEL_01AS)
	$(BUILD)/pack_series series+chiron data/eph-1as-chiron.bin $(LEVEL_1AS)
	$(BUILD)/pack_series chiron data/eph-chiron.bin
	$(BUILD)/pack_series asteroids data/eph-asteroids.bin

# WebAssembly (wasm32-wasi reactor modules exporting the C API plus
# malloc/free; no file I/O, data comes from JavaScript as bytes):
#   build/wasm/ephemeris.wasm              core only
#   build/wasm/ephemeris-1as.wasm          core with data/eph-1as.bin built in
#   build/wasm/ephemeris-1as-chiron.wasm   ... with data/eph-1as-chiron.bin
# WASM_CXX can also be a wasi-sdk clang, e.g.
#   make wasm WASM_CXX="clang++ --target=wasm32-wasi --sysroot=/opt/wasi-sdk/share/wasi-sysroot"
WASM_CXX   ?= zig c++ -target wasm32-wasi
WASM_FLAGS  = -Os -std=c++17 -fno-exceptions -fno-rtti -DEPHEMERIS_NO_FILE_IO \
              -mexec-model=reactor -Wl,--strip-all $(INCLUDES)
comma  := ,
lparen := (
WASM_API   := $(shell sed -n 's/^[a-z].*[ *]\(eph_[a-z0-9_]*\)$(lparen).*/\1/p' \
                include/ephemeris/ephemeris.h)
WASM_LINK   = $(addprefix -Wl$(comma)--export=,$(WASM_API) malloc free)

wasm: $(BUILD)/wasm/ephemeris.wasm $(BUILD)/wasm/ephemeris-1as.wasm \
      $(BUILD)/wasm/ephemeris-1as-chiron.wasm

$(BUILD)/wasm/ephemeris.wasm: $(CORE_SRC) src/embedded_none.cpp
	@mkdir -p $(dir $@)
	$(WASM_CXX) $(WASM_FLAGS) $^ $(WASM_LINK) -o $@

.PRECIOUS: $(BUILD)/wasm/embedded_%.cpp
$(BUILD)/wasm/embedded_%.cpp: data/eph-%.bin tools/embed.cmake
	@mkdir -p $(dir $@)
	cmake -DIN=$< -DOUT=$@ -P tools/embed.cmake

$(BUILD)/wasm/ephemeris-%.wasm: $(CORE_SRC) $(BUILD)/wasm/embedded_%.cpp
	@mkdir -p $(dir $@)
	$(WASM_CXX) $(WASM_FLAGS) $^ $(WASM_LINK) -o $@

$(BUILD)/test_wasm_ref: tests/test_wasm_ref.c $(LIB)
	$(CC) $(CFLAGS) -Iinclude $< $(LIB) -lm -o $@

wasm-test: wasm $(BUILD)/test_wasm_ref
	node tests/test_wasm.mjs $(BUILD) data

clean:
	rm -rf $(BUILD)

-include $(LIB_OBJ:.o=.d)
