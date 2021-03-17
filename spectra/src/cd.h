#ifndef __CD_H__
#define __CD_H__
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "input.h"

void cross(double *a, double *b, double *res);
void dot(double *a, double* b, double res);
void cd(unsigned n_chl, unsigned ns, double **chiw, double **mu,
        double **eig, double **com, double *eigvals, double *cd);

#endif
