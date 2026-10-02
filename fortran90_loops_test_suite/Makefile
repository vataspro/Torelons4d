FC = gfortran
ARCH := $(shell uname -m)

ifeq ($(ARCH),x86_64)
FFLAGS ?= -O2 -mcmodel=medium -ffree-line-length-none
else
FFLAGS ?= -O2 -ffree-line-length-none
endif

BUILD_DIR := build
TARGET := operator_test_suite
OUTPUT := operator_outputs.dat

SOURCES := \
	parameters.f90 \
	read_field_config.f90 \
	lattice.f90 \
	thermal_lines.f90 \
	operator_dispatch.f90 \
	test_suite.f90

.PHONY: all clean run regenerate

all: $(TARGET)

regenerate: operator_dispatch.f90

operator_dispatch.f90: thermal_lines.f90 generate_dispatch.py
	./generate_dispatch.py

$(TARGET): $(SOURCES) | $(BUILD_DIR)
	$(FC) $(FFLAGS) -J$(BUILD_DIR) -I$(BUILD_DIR) $(SOURCES) -o $(TARGET)

run: $(TARGET)
	./$(TARGET)

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

clean:
	rm -rf $(BUILD_DIR) $(TARGET) operator_dispatch.f90 $(OUTPUT)
