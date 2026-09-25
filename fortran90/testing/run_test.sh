#!/bin/bash
gfortran -c ../parameters.f90 -I/../
gfortran -c ../lattice.f90 -I/../
gfortran -c ../read_field_config.f90 -I/../
gfortran -o test_suite test_suite.f90 parameters.o lattice.o read_field_config.o
./test_suite
