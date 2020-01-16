#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef PARAMETERS_H
#define PARAMETERS_H

#define HBAR 2.90283E-23 /* in cm lol */
#define KB 1.3806503E-23 /* but not this? ask chris */

/* possibly useful for pigment switch? */
typedef enum {
    NONE = 0,
    S100 = 1, 
    S110 = 2, 
    S111 = 3 
} Strain;

typedef struct {
    double s0, s1, s2, g0, g1, g2, l0, l1, l2, w1, w2, t, T;
    double (* cw)(double, void *);
} Parameters;

Parameters getParameters(char* filename);

#endif
