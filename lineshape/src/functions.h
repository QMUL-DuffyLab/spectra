#include <complex.h>
#include <stdio.h>
#include "parameters.h"

#ifndef FUNCTIONS_H
#define FUNCTIONS_H

typedef double (*spec_dens_type)(double, void*);

/** Defines which spectral density ansatz to use.
 *
 * OBO    -> simple overdamped Brownian oscillator,
 * Renger -> ansatz from M\"{u}h and Renger derived from FLN data
 * BIG    -> from Novoderzhkin via Mancal - OBO with 48 
 *           high-frequency underdamped modes added on
 * CAR    -> carotenoid ansatz from Chris
 */
typedef enum ansatz {
  OBO = 0,
  RENGER = 1,
  BIG = 2,
  CAR = 3,
} ansatz;

double (*choose_ansatz(ansatz ans))(double, void *);
double get_offset(ansatz ans, Parameters params);
double cw_renger(double w, void *params);
double cw_car(double w, void *params);
double cw_obo(double w, void *params);
double cw_big(double w, void *params);
double c_n(double w, void *params);
double trig_re(double w, void *params);
double trig_im(double w, void *params);
double reorg_int(double w, void *params);
double complex At(double w0, double re, double im,
    double t, double gamma);
double complex Ft(double w0, double re, double im,
    double reorg, double t, double gamma);

#endif
