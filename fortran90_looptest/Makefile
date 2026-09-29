FC = gfortran
ARCH := $(shell uname -m)

ifeq ($(ARCH),x86_64)
FFLAGS ?= -O2 -mcmodel=medium -ffree-line-length-none
else
FFLAGS ?= -O2 -ffree-line-length-none
endif

BUILD_DIR := build
TARGET := torelons_4D_adjoint
LOOP_TEST_TARGET := loop_test_simplified

MAIN_SOURCES := \
	parameters.f90 \
	lattice.f90 \
	luesher.f90 \
	correlator_construction.f90 \
	torelons_4D_adjoint.f90

LOOP_TEST_SOURCES := \
	parameters.f90 \
	read_field_config.f90 \
	lattice.f90 \
	thermal_lines_simplified.f90 \
	loop_test_simplified.f90

.PHONY: all main clean test loop-test-simplified run-loop-test

all: loop-test-simplified

main: $(TARGET)

$(TARGET): $(MAIN_SOURCES) | $(BUILD_DIR)
	$(FC) $(FFLAGS) -J$(BUILD_DIR) -I$(BUILD_DIR) $(MAIN_SOURCES) -o $(TARGET)

test: test_suite

test_suite: parameters.f90 read_field_config.f90 test_suite.f90 | $(BUILD_DIR)
	$(FC) $(FFLAGS) -J$(BUILD_DIR) -I$(BUILD_DIR) parameters.f90 read_field_config.f90 test_suite.f90 -o test_suite

loop-test-simplified: $(LOOP_TEST_TARGET)

$(LOOP_TEST_TARGET): $(LOOP_TEST_SOURCES) | $(BUILD_DIR)
	$(FC) $(FFLAGS) -J$(BUILD_DIR) -I$(BUILD_DIR) $(LOOP_TEST_SOURCES) -o $(LOOP_TEST_TARGET)

run-loop-test: $(LOOP_TEST_TARGET)
	./$(LOOP_TEST_TARGET)

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

clean:
	rm -rf $(BUILD_DIR) $(TARGET) test_suite $(LOOP_TEST_TARGET) new_loop_outputs.dat
