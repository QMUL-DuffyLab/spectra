#include <math.h>
#include <complex.h>
#include "functions.h"

#ifndef M_PI
	#define M_PI 3.1415926535897932384626433832795L
#endif

/** Returns a function pointer to the correct ansatz.
 *
 * lmao. hope you like parsing function pointer definitions!
 * this defines a function choose_ansatz, which takes a parameter
 * of type chl_ansatz - this is an enum defined in functions.h.
 * it returns a pointer to the correct spectral density function,
 * which is a function that takes a double argument (the frequency)
 * and a void pointer (the parameters for the spectral density), and
 * returns a double. Cool stuff eh
 *
 */
double (*choose_ansatz(chl_ansatz ansatz))(double, void *)
{
  double (* ansatz_pointer)(double, void *);
  fprintf(stdout, "chl_ansatz: %d\n", ansatz);
  if (ansatz == OBO) {
    ansatz_pointer = &cw_obo;
  } else if (ansatz == RENGER) {
    ansatz_pointer = &cw_renger;
  } else if (ansatz == BIG) {
    ansatz_pointer = &cw_big;
  } else {
    fprintf(stdout, "chl_ansatz is wrong: %d\n", ansatz);
    exit(EXIT_FAILURE);
  }
  fprintf(stdout, "ansatz_pointer address: %p\n", (void *)&ansatz_pointer);
  return ansatz_pointer;
}

/* Chl spectral density */
double
cw_renger(double w, void* params)
{
    /* Pass a void pointer and cast it here for compatibility
     * with gsl_function when we do the quadrature */
    Parameters *p = (Parameters *) params;
    /* double pf = 2 * M_PI * CMS * 100 * 1E-15; */
    fprintf(stdout, "YOU HAVE ARRIVED AT: RENGER ANSATZ\n");

    /* 7! is 5040; this is the Renger form for chlorophyll */
    double c1 = (p->s1 / (5040 * 2 * pow(p->w1, 4.)))
    	      * (exp(-1. * sqrt(fabs(w) / (p->w1))));
    double c2 = (p->s2 / (5040 * 2 * pow(p->w2, 4.)))
    	      * (exp(-1. * sqrt(fabs(w) / (p->w2))));
    return ((M_PI * p->s0 * pow(w, 5.)) / (p->s1 + p->s2)) * (c1 + c2);
}

/* Car spectral density */
double
cw_car(double w, void* params)
{
    Parameters *p = (Parameters *) params;
    /* ansatz from Kieran's paper on carotenoids */
    double c1 = 2. * p->l1 * (w * p->g1 * pow(p->w1, 2.)) /
    	      (pow((pow(w, 2.) - pow(p->w1, 2.)), 2.)
   	    + (pow(w, 2.) * pow(p->g1, 2.)));
    double c2 = 2. * p->l2 * (w * p->g2 * pow(p->w2, 2.)) /
    	      (pow((pow(w, 2.) - pow(p->w2, 2.)), 2.)
    	    + (pow(w, 2.) * pow(p->g2, 2.)));
    return (c1 + c2 + 2 * p->l0 * (w * p->g0)
    	    / (pow(w, 2.) + pow(p->g0, 2.)));
}

/* overdamped brownian oscillator */
double
cw_obo(double w, void* params)
{
    Parameters *p = (Parameters *) params;
    fprintf(stdout, "YOU HAVE ARRIVED AT: OBO ANSATZ\n");
    /* l0 is the reorganisation energy expressed in cm^{-1},
     * gamma0 is the correlation time of fluctuations in cm^{-1} */
    return (2. * p->l0 * p->g0 * w)/(pow(w, 2.) + pow(p->g0, 2.));
}

/* from novoderezhkin/mancal - one OBO +
 * 48 high frequency vibrational modes */
double
cw_big(double w, void* params)
{
    Parameters *p = (Parameters *) params;
    fprintf(stdout, "YOU HAVE ARRIVED AT: BIG ANSATZ\n");
    double c; 
    double sum = 0.;
    double c0 = 2 * p->l0 * (w * p->g0)
    	    / (pow(w, 2.) + pow(p->g0, 2.));
    for (unsigned int i = 0; i < 48; i++) {
      /* p->gsw[j][i] is effectively 3 * 48 doubles
       * with the g's first, s's second, w's last */
      c = 2 * p->gsw[1][i] * p->gsw[2][i] *
        (pow(p->gsw[2][i], 2.) * p->gsw[0][i] * w)
        / (pow((pow(w, 2.) - pow(p->gsw[2][i], 2.)), 2.) 
        + (pow(w, 2.) * pow(p->gsw[0][i], 2.)));
      sum +=c;
    }
    return c0 + sum;
}

double
c_n(double w, void* params)
{
  /* cw_x is C''(w) - spectral density (odd). C_n from
   * the mancal paper is the sum of even and odd parts - 
   * enters in the Redfield rate equations for k_{ij} */
  Parameters *p = (Parameters *) params;
  double hbeta = 1.439 / p->T;
  return (1. + (1. / tanh(0.5 * w * hbeta))) * p->cw(w, p);
}

/* the integral for g(t) includes the spectral density function,
 * which should be switchable, but we need to call a function pointer
 * of a specific form later on. Hence, add a function pointer to the params
 * struct and use that to decide which spectral density function to use. */
double
trig_re(double w, void* params)
{
    Parameters *p = (Parameters *) params;
    /* (H * C)/KB in cm = 1.439 */
    double hbeta = 1.439 / p->T;
    return (creal(p->cw(w, p)) * (1. / (M_PI * pow(w, 2.)))
	   * (1 - cos(w * p->ti))
	   * (1. / tanh(0.5 * w * hbeta)))
           + (cimag(p->cw(w, p)) * (1. / (M_PI * pow(w, 2.)))
    	   * (sin(w * p->ti) - w * p->ti));
}

double
trig_im(double w, void* params)
{
    Parameters *p = (Parameters *) params;
    double hbeta = 1.439 / p->T;
    return (creal(p->cw(w, p)) * (1. / (M_PI * pow(w, 2.)))
    	   * (sin(w * p->ti) - w * p->ti))
           + (cimag(p->cw(w, p)) * (1. / (M_PI * pow(w, 2.)))
	   * (1 - cos(w * p->ti))
	   * (1. / tanh(0.5 * w * hbeta)));
}

double
reorg_int(double w, void* params)
{
    Parameters *p = (Parameters *) params;
    return p->cw(w, p) * (1. / (M_PI * w));
}

double complex
At(double w0, double re, double im, double t,
   double gamma)
{
    /* l1 and l2 from chris's python code - check */
    double complex exponent = -I * ((w0) * t) - (re + (I * im))
      - (0.5 * gamma * t);
    return cexp(exponent);
}

double complex
Ft(double w0, double re, double im, double reorg,
   double t, double gamma)
{
    /* hbar   = 1.054571817E-34 Js
     * 1J     = 6.24150974E18 eV
     * 1ev/hc ~ 8065.5439 cm^-1
     * so 1J  = 5.034117081E22 cm^-1
     * hbar   = (1.054E-34 * 5.034E22) cm^-1 s
     *        = 5.3088E-12 cm^-1 s
     *        = 5.3088 cm^-1 ps
    */
    double complex exponent = -I * ((w0) + (2. * (reorg))) * t 
                            - (re - (I * im) - (0.5 * t * gamma));
    return cexp(exponent);
}
