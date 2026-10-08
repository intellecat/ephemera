# Convenience Makefile.  CMake is the primary build system; see README.
CXX      ?= g++
CC       ?= gcc
CXXFLAGS ?= -O2 -std=c++17 -Wall -Wextra
CFLAGS   ?= -O2 -std=c99 -Wall -Wextra
INCLUDES  = -Iinclude -Isrc -Idata
BUILD     = build

LIB_SRC = src/time.cpp src/frames.cpp src/vsop87.cpp src/moon.cpp \
          src/elp.cpp src/pluto.cpp src/ephemeris.cpp src/correction.cpp \
          data/vsop87d_data.cpp data/pluto_data.cpp data/elp_data.cpp \
          data/corr_data.cpp
LIB_OBJ = $(patsubst %.cpp,$(BUILD)/%.o,$(LIB_SRC))
LIB     = $(BUILD)/libephemeris.a

.PHONY: all test clean

all: $(LIB) $(BUILD)/ephcliapp $(BUILD)/test_ephemeris

$(BUILD)/%.o: %.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

$(LIB): $(LIB_OBJ)
	ar rcs $@ $^

$(BUILD)/ephcliapp: examples/positions.c $(LIB)
	$(CC) $(CFLAGS) -Iinclude $< $(LIB) -lm -o $@

$(BUILD)/test_ephemeris: tests/test_ephemeris.cpp $(LIB)
	$(CXX) $(CXXFLAGS) -Iinclude -Itests $< $(LIB) -lm -o $@

test: $(BUILD)/test_ephemeris
	$(BUILD)/test_ephemeris

clean:
	rm -rf $(BUILD)
