# Torelon

Fortran code for the calculation of torelon operator correlation
functions on gauge field configurations.

## Compiling

To compile the modern fortran code navigate to the directory 
and use the provided Makefile:

```bash
cd fortran90
make -j2
```

## Usage

The executable must be called from a folder where a file
named `parameters.txt` should be present. An example
parameter file is present in the `fortran90` directory.

```bash
mkdir run
cd run
cp ../parameter_file.txt parameters.txt
../torelons_4D_adjoint
```

## Legacy code

The legacy fortran77 code can also be compiled using the
Makefile in the relevant subdirectory. We advise that the
loading of the configurations is hard-coded in the source
code and suggest that users may preffer to use the modern
code.

```
cd fortran77_scalar
make
```
