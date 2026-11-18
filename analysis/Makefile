
FC = gfortran
FFLAGS = -O2

BUILDDIR = build

# Prefer standard 32-bit-index LAPACK/BLAS.
# Fall back to 64-bit-index libraries.
LAPACK_LIBS := $(shell \
    printf 'program test\nend program test\n' | \
    $(FC) -x f95 - -o /dev/null -llapack -lblas \
    >/dev/null 2>&1 && \
    echo '-llapack -lblas' || \
    echo '-llapack64 -lblas64')

.PHONY: all clean

all: analysis

# Create build directory
$(BUILDDIR):
	mkdir -p $(BUILDDIR)

# Compile linalg.f90 into build/
$(BUILDDIR)/linalg.o: linalg.f90 | $(BUILDDIR)
	$(FC) $(FFLAGS) -J$(BUILDDIR) -I$(BUILDDIR) -c $< -o $@

# Compile analysis.f90 into build/
# Depends on linalg.o in case it uses a module from linalg.
$(BUILDDIR)/analysis.o: analysis.f90 $(BUILDDIR)/linalg.o
	$(FC) $(FFLAGS) -J$(BUILDDIR) -I$(BUILDDIR) -c $< -o $@

# Link executable in parent directory
analysis: $(BUILDDIR)/analysis.o $(BUILDDIR)/linalg.o
	$(FC) $(FFLAGS) -o $@ $^ $(LAPACK_LIBS)

# Remove generated files
clean:
	rm -rf $(BUILDDIR) analysis
