#include "fluorescence.h"
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_math.h>

#ifndef __STEADY_STATE_H__
#define __STEADY_STATE_H__

/** Defines the shape of the incident light pulse.
 *
 * The population equations for steady-state fluorescence include
 * a source term convolved with the exciton spectra, written as 
 * \f$ \int W(\omega) \chi_i(\omega) d\omega \f$
 * where \f$ \chi_i(\omega) \f$ is the absorption spectrum of the i'th
 * exciton and \f$ W(\omega) \f$ is the incident light.
 * This enum defines the shape of the pulse - either white light
 * (flat), Lorentzian, Gaussian, or narrow laser illumination.
 */
typedef enum pulse_type {
  FLAT = 0,
  LORENTZIAN = 1,
  GAUSSIAN = 2,
  DELTA = 3,
} pulse_type;

/** Additional parameters for incident light term.
 *
 * If the pulse is not `FLAT` then we need to know where its centre is:
 * for `DELTA` this is enough but for `LORENTZIAN` or `GAUSSIAN` we
 * also need the width. We do not need the peak height because I
 * normalise the integral to 1, as will be explained in the 
 * documentation for incident().
 *
 * For `LORENTZIAN` we have the form
 * \f[
      \mathcal{L} = \frac{1}{\pi \gamma} 
      \frac{\gamma^2}{(x - x_0)^2 + \gamma^2}
 * \f]
 * where `centre` \f$ = x_0 \f$ and `width` \f$ = \gamma \f$.
 *
 * For `GAUSSIAN` we have
 * \f[
      \mathcal{N} = \frac{1}{\sigma \sqrt{2 \pi}} 
      \exp{-\frac{1}{2} {(\frac{x - \mu}{\sigma})}^2}
 * \f]
 * where `centre` \f$ = \mu \f$ and `width` \f$ = \sigma \f$.
 *
 * `centre` and `width` should be given in \f$ cm^{-1} \f$.
 *
 */
typedef struct pulse {
  pulse_type type;
  size_t target_state;
  double amplitude;
  double centre;
  double width;
  double t_peak;
  double duration;
} pulse;

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
double* incident(pulse p, unsigned int tau);
double  intensity(double w, double t, pulse p);
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
