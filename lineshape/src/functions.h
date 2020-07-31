#include <complex.h>
#include <stdio.h>
#include "parameters.h"

#ifndef FUNCTIONS_H
#define FUNCTIONS_H

double cw_chl(double w, void *params);
double cw_car(double w, void *params);
double cw_odo(double w, void *params);
double cw_big(double w, void *params);
double c_n(double w, void *params);
double trig_re(double w, void *params);
double trig_im(double w, void *params);
double reorg_int(double w, void *params);
double complex At(double w, double w0, double re, double im,
    double t, double gamma);
double complex Ft(double w, double w0, double re, double im,
    double reorg, double t, double gamma);

#endif
