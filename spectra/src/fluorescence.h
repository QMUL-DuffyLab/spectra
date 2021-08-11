#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <vector>
#include "input.h"
#include "invert_matrix.h"
#include "../../lineshape/src/parameters.h"

#ifndef __FLUORESCENCE_H__
#define __FLUORESCENCE_H__

typedef struct {
  unsigned int N;
  double **kij;
  double **Tij;
  double *rates;
  double *chiw;
} ode_params;

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


double* incident(pulse p, unsigned int tau);
void Gaussian(unsigned n, double* w,
        double Delta_j0_ba, double Delta_w, double *gau);
double intensity(double w, double t, pulse p);
gsl_matrix* array_to_gsl_matrix(unsigned int n1,
            unsigned int n2, double** mat);
void check_detailed_balance(unsigned n, double t, double thresh,
                            double **kij, double **wij);
double** rate_calc (unsigned int N, double **eig,
                    double** wij, Parameters *p);
double* relaxation_rates (unsigned int N, double* gamma);
double **transfer_matrix (unsigned int N, double* relax, double** kij);
int odefunc (double x, const double *y, double *f, void *params);
int jacobian (double t, const double y[], double *dfdy,
              double dfdt[], void *params);
void bcs (unsigned const int N, const double* eigvals,
             const double T, double *res);
double trapezoid(double *f, double dx, unsigned int n);
void decompose_transfer_matrix(unsigned n, double **Tij,
                               double **Tij_vr, double **Tij_vr_inv,
                               double **Tij_wr);
double mean_excitation_lifetime(unsigned n, double **Tij_vr,
                                double **Tij_vr_inv,
                                double **Tij_wr, double *p0);
void population(unsigned n, double t, double *pt, double **Tij_vr,
                double **Tij_vr_inv, double **Tij_wr, double *p0);
#endif
