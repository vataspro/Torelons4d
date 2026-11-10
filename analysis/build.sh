#!/bin/bash
gfortran -O2 -c linalg.f90
gfortran -O2 -o analysis analysis.f90 linalg.o -llapack -lblas
