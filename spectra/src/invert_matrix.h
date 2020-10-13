#include <stdio.h>
#include <lapacke.h>

#ifndef __INVERT_MATRIX_H__
#define __INVERT_MATRIX_H__

void
print_matrix(char* name, unsigned n, double **matrix);

void
interchange_2d_contiguous(char direction, unsigned n,
              double** array_2d, double* contiguous);

int
invert_matrix_oop(unsigned n, double** in, double** out);

int
invert_matrix_ip(unsigned n, double** inout);

#endif
