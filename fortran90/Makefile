FC = gfortran
ARCH := $(shell uname -m)

ifeq ($(ARCH),x86_64)
FFLAGS ?= -O2 -mcmodel=medium
else
FFLAGS ?= -O2
endif

BUILD_DIR := build
TARGET := torelons_4D_adjoint

SOURCES := \
	parameters.f90 \
	lattice.f90 \
	luesher.f90 \
	correlator_construction.f90 \
	torelons_4D_adjoint.f90

.PHONY: all clean test

all: $(TARGET)

$(TARGET): $(SOURCES) | $(BUILD_DIR)
	$(FC) $(FFLAGS) -J$(BUILD_DIR) -I$(BUILD_DIR) $(SOURCES) -o $(TARGET)

test: test_suite

test_suite: parameters.f90 read_field_config.f90 test_suite.f90 | $(BUILD_DIR)
	$(FC) $(FFLAGS) -J$(BUILD_DIR) -I$(BUILD_DIR) parameters.f90 read_field_config.f90 test_suite.f90 -o test_suite

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

clean:
	rm -rf $(BUILD_DIR) $(TARGET) test_suite
