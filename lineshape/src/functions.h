#include <complex.h>
#include "parameters.h"

#ifndef FUNCTIONS_H
#define FUNCTIONS_H

double cw_chl(double w, void *params);
double cw_car(double w, void *params);
double cw_odo(double w, void *params);
double trig_re(double w, void *params);
double trig_im(double w, void *params);
double reorg_int(double w, void *params);
double complex At(double w0, double re, double im, double t, double l1, double l2);
double complex Ft(double w0, double re, double im, double reorg, double t);

#endif