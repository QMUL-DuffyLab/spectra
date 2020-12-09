#include <complex.h>
#include <stdio.h>
#include "parameters.h"

#ifndef FUNCTIONS_H
#define FUNCTIONS_H

typedef double (*spec_dens_type)(double, void*);

double (*choose_ansatz(ansatz ans))(double, void *);
double get_offset(Parameters params);
double cw_renger(double w, void *params);
double cw_car(double w, void *params);
double cw_obo(double w, void *params);
double cw_big(double w, void *params);
double c_n(double w, void *params);
double trig_re(double w, void *params);
double trig_im(double w, void *params);
double reorg_int(double w, void *params);
double _Complex At(double w0, double re, double im,
    double t, double gamma);
double _Complex Ft(double w0, double re, double im,
    double reorg, double t, double gamma);

#endif
