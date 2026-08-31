# Cross-section tools standalone build
# Linux only, standard C++17

include config.mk

CXX      ?= g++
CXXFLAGS  = -std=c++17 -O2 -Wall -Wextra -Wpedantic

TARGETS   = $(BUILD_DIR)/CrossSectionInterpolator \
            $(BUILD_DIR)/CrossSectionScaler \
            $(BUILD_DIR)/SpokeAzimuthalDensity

.PHONY: all clean

all: $(TARGETS) $(BUILD_DIR)/paths.conf

$(BUILD_DIR)/CrossSectionInterpolator: CrossSectionInterpolator.cpp
	@mkdir -p $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -o $@ $<

$(BUILD_DIR)/CrossSectionScaler: CrossSectionScaler.cpp
	@mkdir -p $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -o $@ $<

$(BUILD_DIR)/SpokeAzimuthalDensity: SpokeAzimuthalDensity.cpp
	@mkdir -p $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -o $@ $<

$(BUILD_DIR)/paths.conf: paths.conf
	@mkdir -p $(BUILD_DIR)
	cp $< $@

clean:
	rm -f $(TARGETS) $(BUILD_DIR)/paths.conf
