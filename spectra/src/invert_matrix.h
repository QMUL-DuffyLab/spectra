#include <stdio.h>
#include <lapacke.h>

#ifndef __INVERT_MATRIX_H__
#define __INVERT_MATRIX_H__

int
invert_matrix_oop(unsigned n, double** in, double** out);

int
invert_matrix_ip(unsigned n, double** inout);

#endif
