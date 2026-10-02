FC = gfortran
ARCH := $(shell uname -m)

ifeq ($(ARCH),x86_64)
FFLAGS ?= -O3 -mcmodel=medium -ffree-line-length-none
else
FFLAGS ?= -O3 -ffree-line-length-none
endif

BUILD_DIR := build
TARGET := torelons_4D_adjoint

SOURCES := \
	parameters.f90 \
	lattice.f90 \
	read_field_config.f90 \
	thermal_lines.f90 \
	here_be_dragons.f90 \
	states_class.f90 \
	main.f90

OBJECTS := $(SOURCES:%.f90=$(BUILD_DIR)/%.o)

.PHONY: all clean test exec

all: $(TARGET)

# ----------------------------------------------------------------------
# Main executable
# ----------------------------------------------------------------------

$(TARGET): $(OBJECTS)
	$(FC) $(FFLAGS) $^ -o $@

# ----------------------------------------------------------------------
# Compile each Fortran source separately.
#
# -J$(BUILD_DIR) puts generated .mod files in build/
# -I$(BUILD_DIR) tells gfortran where to find .mod files
# ----------------------------------------------------------------------

$(BUILD_DIR)/%.o: %.f90 | $(BUILD_DIR)
	$(FC) $(FFLAGS) -J$(BUILD_DIR) -I$(BUILD_DIR) -c $< -o $@

# ----------------------------------------------------------------------
# Fortran module dependencies
# ----------------------------------------------------------------------

$(BUILD_DIR)/lattice.o: $(BUILD_DIR)/parameters.o

$(BUILD_DIR)/thermal_lines.o: \
	$(BUILD_DIR)/parameters.o \
	$(BUILD_DIR)/lattice.o

$(BUILD_DIR)/read_field_config.o: \
	$(BUILD_DIR)/parameters.o \
	$(BUILD_DIR)/lattice.o

$(BUILD_DIR)/here_be_dragons.o: \
	$(BUILD_DIR)/parameters.o \
	$(BUILD_DIR)/lattice.o \
	$(BUILD_DIR)/thermal_lines.o \

$(BUILD_DIR)/states_class.o: \
	$(BUILD_DIR)/parameters.o \
    $(BUILD_DIR)/lattice.o

$(BUILD_DIR)/main.o: \
	$(BUILD_DIR)/parameters.o \
	$(BUILD_DIR)/lattice.o \
	$(BUILD_DIR)/here_be_dragons.o \
	$(BUILD_DIR)/states_class.o

# ----------------------------------------------------------------------
# Extra executable
# ----------------------------------------------------------------------

exec: parameters.f90 lattice.f90 read_field_config.f90 thermal_lines.f90 mock_main.f90 | $(BUILD_DIR)
	$(FC) $(FFLAGS) -J$(BUILD_DIR) -I$(BUILD_DIR) \
		parameters.f90 lattice.f90 read_field_config.f90 \
		thermal_lines.f90 mock_main.f90 -o exec

# ----------------------------------------------------------------------
# Tests
# ----------------------------------------------------------------------

test: test_suite

test_suite: parameters.f90 read_field_config.f90 test_suite.f90 | $(BUILD_DIR)
	$(FC) $(FFLAGS) -J$(BUILD_DIR) -I$(BUILD_DIR) \
		parameters.f90 read_field_config.f90 test_suite.f90 -o test_suite

# ----------------------------------------------------------------------
# Build directory
# ----------------------------------------------------------------------

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

# ----------------------------------------------------------------------
# Clean
# ----------------------------------------------------------------------

clean:
	rm -rf $(BUILD_DIR) $(TARGET) test_suite exec
