#include "fluorescence.h"
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_math.h>

#ifndef __STEADY_STATE_H__
#define __STEADY_STATE_H__

typedef enum pulse_type {
  FLAT = 0,
  LORENTZIAN = 1,
  GAUSSIAN = 2,
  DELTA = 3,
} pulse_type;

typedef struct pulse {
  pulse_type type;
  /* centre and width should be in cm^{-1}
   * for LORENTZIAN they're x_0 and \gamma (FWHM)
   * for GAUSSIAN they're \mu and \sigma
   * for DELTA width is meaningless */
  double centre;
  double width;
} pulse;

typedef enum ss_init {
  /* decides what the P_i(0) should look like as the initial input
   * for the steady-state solver (multidim. root finder in GSL) */
    BOLTZ = 0,      /* straight Boltzmann factors */
    BOLTZ_MUSQ = 1, /* Boltzmann weighted by oscillator strength */
    CONST = 2,      /* flat probability distribution */
    MAX = 3,        /* excite only state with highest |μ^2| */
} ss_init;

unsigned short int pop_converge(double *y, double *yprev, 
                                unsigned int N, double thresh);
double* incident(pulse p, unsigned int tau);
gsl_vector* guess(const ss_init p,
                  const double* boltz, const double* musq,
                  unsigned const int max, unsigned const int N);
int pop_steady_f   (const gsl_vector *x, void *params, gsl_vector *f);
int pop_steady_df  (const gsl_vector *x, void *params, gsl_matrix *J);
int pop_steady_fdf (const gsl_vector *x, void *params,
                    gsl_vector *f, gsl_matrix *J);

#endif
