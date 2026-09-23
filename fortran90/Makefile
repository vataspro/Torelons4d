FC = gfortran
ARCH := $(shell uname -m)

ifeq ($(ARCH),x86_64)
FFLAGS ?= -O2 -mcmodel=medium -ffree-line-length-none
else
FFLAGS ?= -O2 -ffree-line-length-none
endif

BUILD_DIR := build
TARGET := torelons_4D_adjoint
LATTICE_BUILD := $(BUILD_DIR)/lattice_gfortran.f90

SOURCES := \
	parameters.f90 \
	$(LATTICE_BUILD) \
	luesher.f90 \
	correlator_construction.f90 \
	torelons_4D_adjoint.f90

.PHONY: all clean test

all: $(TARGET)

$(TARGET): $(SOURCES) | $(BUILD_DIR)
	$(FC) $(FFLAGS) -J$(BUILD_DIR) -I$(BUILD_DIR) $(SOURCES) -o $(TARGET)

$(LATTICE_BUILD): lattice.f90 | $(BUILD_DIR)
	cp lattice.f90 $@
	perl -pi -e 's/A_dagger_A\(q,q\)%re/real(A_dagger_A(q,q), kind=real64)/g; s/A_dagger_A\(p,p\)%re/real(A_dagger_A(p,p), kind=real64)/g; s/A_dagger_A\(i,i\)%re/real(A_dagger_A(i,i), kind=real64)/g' $@

test: test_suite

test_suite: parameters.f90 read_field_config.f90 test_suite.f90 | $(BUILD_DIR)
	$(FC) $(FFLAGS) -J$(BUILD_DIR) -I$(BUILD_DIR) parameters.f90 read_field_config.f90 test_suite.f90 -o test_suite

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

clean:
	rm -rf $(BUILD_DIR) $(TARGET) test_suite
