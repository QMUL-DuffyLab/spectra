#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <fftw3.h>
#include "parameters.h"

#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#ifndef M_PI
	#define M_PI 3.1415926535897932384626433832795L
#endif

#ifdef __cplusplus
  #include <cmath>
  #include <complex>
  /* using namespace std; */
  #define COMPLEX std::complex<double>
  #define REAL(X) std::real(X)
  #define IMAG(X) std::imag(X)
  #define II std::complex<double>(0.0, 1.0)
  #define EXP(X) exp(X)
#else
  #define COMPLEX double _Complex
  #define REAL(X) creal(X)
  #define IMAG(X) cimag(X)
  #define II CMPLX(0.0, 1.0)
  #define EXP(X) cexp(X)
#endif


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
COMPLEX At(double w0, double re, double im,
    double t, double gamma);
COMPLEX Ft(double w0, double re, double im,
    double reorg, double t, double gamma);

#endif
