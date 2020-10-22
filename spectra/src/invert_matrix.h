#include <stdio.h>
#include <lapacke.h>

#ifndef __INVERT_MATRIX_H__
#define __INVERT_MATRIX_H__

void
print_matrix(FILE *stream, char* name, unsigned n, double** matrix);

void
print_vector(FILE *stream, char* name, unsigned n, double *vector);

void
interchange_2d_contiguous(char direction, unsigned n,
              double** array_2d, double* contiguous);

int
invert_matrix_oop(unsigned n, double** in, double** out);

int
invert_matrix_ip(unsigned n, double** inout);

int
eig_oop(unsigned n, double** in, double** vr, double* eigvals);

int
eig_ip(unsigned n, double** inout, double* eigvals);

void matmul(unsigned n, double **a, double **b, double **ab);
void matvec(unsigned n, double **mat, double *vec, double *res);

#endif
