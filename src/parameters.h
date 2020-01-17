#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>

#ifndef PARAMETERS_H
#define PARAMETERS_H

#define CMS 299792458 /* speed of light _C_ in _M_etres per _S_econd */

typedef struct {
    long unsigned int ns;
    double T;
    char aw_file[200];
    char fw_file[200];
} Protocol;

typedef struct {
    double s0, s1, s2, g0, g1, g2, l0, l1, l2, w1, w2, t, T;
    double (* cw)(double, void *);
    char ligand; /* 1 for chlorophyll, 0 for carotenoid */
} Parameters;

Parameters getParameters(char* filename);
Protocol getProtocol(char* filename);

#endif
