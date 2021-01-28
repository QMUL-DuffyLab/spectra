#include "fluorescence.h"
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_math.h>

#ifndef __STEADY_STATE_H__
#define __STEADY_STATE_H__


/** Determines initial conditions for steady-state solver.
 *
 * For steady state fluorescence we need to solve a coupled set of
 * equations - I use GSL root finding for this. Obviously in order
 * for this to work we need to specify initial conditions - this enum
 * switches between them. options are `BOLTZ` for straight Boltzmann
 * factors of the different exciton states, 
 * `BOLTZ_MUSQ` for Boltzmann factors weighted by the excitons'
 * oscillator strength \f$ {|\mu|}^2 \f$, `FLAT` for all the
 * probabilities equal, `MAX` for only the exciton state with the
 * highest oscillator strength being excited.
 *
 * All of these have the normalisation 
 * \f$ \sum_{i = 1}^N P_i(0) = 1 \f$
 * as expected.
 */
typedef enum ss_init {
  /* decides what the P_i(0) should look like as the initial input
   * for the steady-state solver (multidim. root finder in GSL) */
    BOLTZ = 0,      /* straight Boltzmann factors */
    MUSQ = 1,       /* straight oscillator strength */
    BOLTZ_MUSQ = 2, /* Boltzmann weighted by oscillator strength */
    CONST = 3,      /* flat probability distribution */
    MAX = 4,        /* excite only state with highest |μ^2| */
} ss_init;

unsigned short int pop_converge(double *y, double *yprev, 
                                unsigned int N, double thresh);
double* guess(const ss_init p,
                  const double* boltz, const double* musq,
                  unsigned const int max, unsigned const int N);
int pop_steady_f   (const gsl_vector *x, void *params, gsl_vector *f);
int pop_steady_df  (const gsl_vector *x, void *params, gsl_matrix *J);
int pop_steady_fdf (const gsl_vector *x, void *params,
                    gsl_vector *f, gsl_matrix *J);
double* steady_state_populations(double *guess,
                                 void *params, unsigned n);

#endif
