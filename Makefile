# Convenience Makefile.  CMake is the primary build system; see README.
CXX      ?= g++
CC       ?= gcc
CXXFLAGS ?= -O2 -std=c++17 -Wall -Wextra
CFLAGS   ?= -O2 -std=c99 -Wall -Wextra
INCLUDES  = -Iinclude -Isrc -Idata
BUILD     = build

LIB_SRC = src/time.cpp src/frames.cpp src/vsop87.cpp src/moon.cpp \
          src/elp.cpp src/pluto.cpp src/ephemeris.cpp src/correction.cpp \
          src/compressed.cpp src/houses.cpp src/asteroids.cpp \
          data/vsop87d_data.cpp data/pluto_data.cpp data/elp_data.cpp \
          data/corr_data.cpp data/corr_moon.cpp data/asteroid_data.cpp
LIB_OBJ = $(patsubst %.cpp,$(BUILD)/%.o,$(LIB_SRC))
LIB     = $(BUILD)/libephemeris.a

.PHONY: all test clean

all: $(LIB) $(BUILD)/ephcliapp $(BUILD)/test_ephemeris $(BUILD)/test_compressed $(BUILD)/test_topocentric $(BUILD)/test_houses $(BUILD)/test_asteroids

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

test: $(BUILD)/test_ephemeris $(BUILD)/test_compressed $(BUILD)/test_topocentric \
      $(BUILD)/test_houses $(BUILD)/test_asteroids
	$(BUILD)/test_ephemeris
	$(BUILD)/test_compressed data/ephem.bin
	$(BUILD)/test_topocentric
	$(BUILD)/test_houses
	$(BUILD)/test_asteroids

$(BUILD)/test_compressed: tests/test_compressed.cpp $(LIB)
	$(CXX) $(CXXFLAGS) -Iinclude -Itests $< $(LIB) -lm -o $@

$(BUILD)/test_topocentric: tests/test_topocentric.cpp tests/topo_ref.h $(LIB)
	$(CXX) $(CXXFLAGS) -Iinclude -Itests $< $(LIB) -lm -o $@

$(BUILD)/test_houses: tests/test_houses.cpp $(LIB)
	$(CXX) $(CXXFLAGS) -Iinclude -Itests $< $(LIB) -lm -o $@

$(BUILD)/test_asteroids: tests/test_asteroids.cpp tests/asteroid_ref.h $(LIB)
	$(CXX) $(CXXFLAGS) -Iinclude -Itests $< $(LIB) -lm -o $@

clean:
	rm -rf $(BUILD)
