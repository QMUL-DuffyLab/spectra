PLACEHOLDER README
==================

do `gfortran -g coupling_calc.f08 -llapack -std=f2008 -ffree-form -Wall -fcheck=bounds -o coupling_calc` to make; I haven't bothered with a makefile yet because it's one f08 file. Also requires gfortran and LAPACK obviously. to run it should be sufficient to do `./run.py` which builds an input file, makes output dirs and runs the fortran code.
