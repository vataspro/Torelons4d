#!/bin/bash
gfortran -c parameters.f90 read_field_config.f90
gfortran -o test_suite test_suite.f90 parameters.o read_field_config.o
./test_suite