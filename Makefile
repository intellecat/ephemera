# Convenience Makefile.  CMake is the primary build system; see README.
CXX      ?= g++
CC       ?= gcc
CXXFLAGS ?= -O2 -std=c++17 -Wall -Wextra
CFLAGS   ?= -O2 -std=c99 -Wall -Wextra
INCLUDES  = -Iinclude -Isrc -Idata
BUILD     = build

LIB_SRC = src/time.cpp src/frames.cpp src/vsop87.cpp src/moon.cpp \
          src/elp.cpp src/pluto.cpp src/ephemeris.cpp src/correction.cpp \
          src/compressed.cpp src/houses.cpp src/asteroids.cpp src/data.cpp \
          data/delta_t_data.cpp
LIB_OBJ = $(patsubst %.cpp,$(BUILD)/%.o,$(LIB_SRC))
LIB     = $(BUILD)/libephemeris.a

.PHONY: all test tools data clean

TESTS = $(BUILD)/test_ephemeris $(BUILD)/test_compressed $(BUILD)/test_topocentric \
        $(BUILD)/test_houses $(BUILD)/test_asteroids $(BUILD)/test_data

all: $(LIB) $(BUILD)/ephcliapp $(TESTS)

$(BUILD)/%.o: %.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

$(LIB): $(LIB_OBJ)
	rm -f $@
	ar rcs $@ $^

$(BUILD)/ephcliapp: examples/positions.c $(LIB)
	$(CC) $(CFLAGS) -Iinclude $< $(LIB) -lm -o $@

$(BUILD)/test_ephemeris: tests/test_ephemeris.cpp tests/reference_data.h $(LIB)
	$(CXX) $(CXXFLAGS) -Iinclude -Itests $< $(LIB) -lm -o $@

test: $(TESTS)
	$(BUILD)/test_ephemeris data/eph-full.bin
	$(BUILD)/test_compressed data/ephem.bin
	$(BUILD)/test_topocentric data/eph-full.bin
	$(BUILD)/test_houses
	$(BUILD)/test_asteroids data/eph-full.bin data/eph-asteroids.bin
	$(BUILD)/test_data data

$(BUILD)/test_compressed: tests/test_compressed.cpp $(LIB)
	$(CXX) $(CXXFLAGS) -Iinclude -Itests $< $(LIB) -lm -o $@

$(BUILD)/test_topocentric: tests/test_topocentric.cpp tests/topo_ref.h $(LIB)
	$(CXX) $(CXXFLAGS) -Iinclude -Itests $< $(LIB) -lm -o $@

$(BUILD)/test_houses: tests/test_houses.cpp $(LIB)
	$(CXX) $(CXXFLAGS) -Iinclude -Itests $< $(LIB) -lm -o $@

$(BUILD)/test_asteroids: tests/test_asteroids.cpp tests/asteroid_ref.h $(LIB)
	$(CXX) $(CXXFLAGS) -Iinclude -Itests $< $(LIB) -lm -o $@

$(BUILD)/test_data: tests/test_data.cpp tests/reference_data.h $(LIB)
	$(CXX) $(CXXFLAGS) -Iinclude -Itests $< $(LIB) -lm -o $@

# Development tools and the series data files (see README).
SERIES_SRC = data/vsop87d_data.cpp data/elp_data.cpp data/pluto_data.cpp \
             data/corr_data.cpp data/corr_moon.cpp data/asteroid_data.cpp

tools: $(BUILD)/pack_series

$(BUILD)/pack_series: tools/pack_series.cpp $(SERIES_SRC)
	@mkdir -p $(BUILD)
	$(CXX) $(CXXFLAGS) $(INCLUDES) $^ -o $@

# Truncation thresholds (radians) per planet, Mercury..Neptune, and for the
# Moon, chosen so each level stays within its accuracy over 1850-2140.
LEVEL_01AS = 1e-8,1e-8,3e-9,3e-9,1e-8,3e-9,1e-8,1e-8 2e-8
LEVEL_1AS  = 3e-7,1e-7,3e-8,3e-8,3e-7,1e-7,1e-7,3e-7 5e-8

data: $(BUILD)/pack_series
	$(BUILD)/pack_series series data/eph-full.bin 0 0
	$(BUILD)/pack_series series data/eph-01as.bin $(LEVEL_01AS)
	$(BUILD)/pack_series series data/eph-1as.bin $(LEVEL_1AS)
	$(BUILD)/pack_series asteroids data/eph-asteroids.bin

clean:
	rm -rf $(BUILD)
