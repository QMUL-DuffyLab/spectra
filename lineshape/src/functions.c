#include <math.h>
#include <complex.h>
#include "functions.h"

#ifndef M_PI
	#define M_PI 3.1415926535897932384626433832795L
#endif

/* Chl spectral density */
double complex
cw_chl(double w, void* params)
{
    /* Pass a void pointer and cast it here for compatibility
     * with gsl_function when we do the quadrature */
    Parameters *p = (Parameters *) params;
    /* double pf = 2 * M_PI * CMS * 1E-15; */

    /* 7! is 5040; this is the Renger form for chlorophyll */
    double complex c1 = (p->s1 / (5040 * 2 * pow(p->w1, 4.)))
    	      * (cexp(-1. * csqrt((w) / (p->w1))));
    double complex c2 = (p->s2 / (5040 * 2 * pow(p->w2, 4.)))
    	      * (cexp(-1. * csqrt((w) / (p->w2))));
    return ((M_PI * p->s0 * pow(w, 5.)) / (p->s1 + p->s2)) * (c1 + c2);
}

/* Car spectral density */
double complex
cw_car(double w, void* params)
{
    Parameters *p = (Parameters *) params;
    /* ansatz from Kieran's paper on carotenoids */
    double complex c1 = 2. * p->l1 * (w * p->g1 * pow(p->w1, 2.)) /
    	      (pow((pow(w, 2.) - pow(p->w1, 2.)), 2.)
   	    + (pow(w, 2.) * pow(p->g1, 2.)));
    double complex c2 = 2. * p->l2 * (w * p->g2 * pow(p->w2, 2.)) /
    	      (pow((pow(w, 2.) - pow(p->w2, 2.)), 2.)
    	    + (pow(w, 2.) * pow(p->g2, 2.)));
    return (c1 + c2 + 2 * p->l0 * (w * p->g0)
    	    / (pow(w, 2.) + pow(p->g0, 2.)));
}

/* overdamped brownian oscillator */
double complex
cw_odo(double w, void* params)
{
    Parameters *p = (Parameters *) params;
    /* l0 is the reorganisation energy expressed in cm^{-1},
     * gamma0 is the correlation time of fluctuations in cm^{-1} */
    return (2. * p->l0 * p->g0 * w)/(pow(w, 2.) + pow(p->g0, 2.));
}

double complex
cn(double w, void* params)
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
   double l1, double l2, double gamma)
{
    /* l1 and l2 from chris's python code - check */
    double complex exponent = -I * (w0 * t) - (re + (I * im))
      - I * t * (l1 + l2) - (0.5 * gamma * t);
    return cexp(exponent);
}

double complex
Ft(double w0, double re, double im, double reorg, double t)
{
    double complex exponent = -I * ((w0 * t) - (2. * reorg)) - (re - (I * im));
    return cexp(exponent);
}
