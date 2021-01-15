#include "functions.h"

/** Returns a function pointer to the correct ansatz.
 *
 * lmao. hope you like parsing function pointer definitions!
 * this defines a function choose_ansatz, which takes a parameter
 * of type ansatz - this is an enum defined in functions.h.
 * it returns a pointer to the correct spectral density function,
 * which is a function that takes a double argument (the frequency)
 * and a void pointer (the parameters for the spectral density), and
 * returns a double. Cool stuff eh
 *
 */
double (*choose_ansatz(ansatz ans))(double, void *)
{
  double (* ansatz_pointer)(double, void *);
  fprintf(stdout, "chl_ans: %d\n", ans);
  if (ans == OBO) {
    ansatz_pointer = &cw_obo;
  } else if (ans == RENGER) {
    ansatz_pointer = &cw_renger;
  } else if (ans == BIG) {
    ansatz_pointer = &cw_big;
  } else if (ans == CAR) {
    ansatz_pointer = &cw_car;
  } else {
    fprintf(stdout, "ansatz is wrong: %d\n", ans);
    exit(EXIT_FAILURE);
  }
  return ansatz_pointer;
}

/** Return the anomalous offset for this pigment with these parameters.
 *
 * We calculate the reorganisation energy on the fly in this code (see
 * the function reorg_int further down), but the mathematical definition
 * \f[
      \lambda = \frac{1}{\pi} \int_{0}^{\infty} \frac{C''(\omega)}{\omega}
 * \f]
 * doesn't give the correct result (NB: what is the exact reason for this?
 * is it something to do with integrating the spectral density instead of
 * a time-domain correlation function?) because there are extra terms in
 * the ansatzes we use.
 *
 * The problem is that these extra terms, which I call the anomalous offset,
 * are not the same for each ansatz; different things contribute depending
 * on what pigment we're looking at. This function gives the correct
 * result for a given pigment.
 */
double
get_offset(Parameters params)
{
  double offset = 0.;
  if (params.ans == OBO) {
    /* do nothing - the anomalous offset in this case is 0 */
  } else if (params.ans  == RENGER) {
    /* NB: haven't used this ansatz in ages - can't remember what
     * the offset is. check! */
  } else if (params.ans == BIG) {
    for (unsigned int i = 0; i < 48; i++) {
      offset += params.gsw[1][i] * params.gsw[2][i];
    }
  } else if (params.ans == CAR) {
    offset = params.l1 + params.l2;
  }
  return offset;
}

/* Chl spectral density */
double
cw_renger(double w, void* params)
{
    /* Pass a void pointer and cast it here for compatibility
     * with gsl_function when we do the quadrature */
    Parameters *p = (Parameters *) params;

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
    return (REAL(p->cw(w, p)) * (1. / (M_PI * pow(w, 2.)))
	   * (1 - cos(w * p->ti))
	   * (1. / tanh(0.5 * w * hbeta)))
           + (IMAG(p->cw(w, p)) * (1. / (M_PI * pow(w, 2.)))
    	   * (sin(w * p->ti) - w * p->ti));
    /* return (creal(p->cw(w, p)) * (1. / (M_PI * pow(w, 2.))) */
	   /* * (1 - cos(w * p->ti)) */
	   /* * (1. / tanh(0.5 * w * hbeta))) */
    /*        + (cimag(p->cw(w, p)) * (1. / (M_PI * pow(w, 2.))) */
    /* 	   * (sin(w * p->ti) - w * p->ti)); */
}

double
trig_im(double w, void* params)
{
    Parameters *p = (Parameters *) params;
    double hbeta = 1.439 / p->T;
    return (REAL(p->cw(w, p)) * (1. / (M_PI * pow(w, 2.)))
    	   * (sin(w * p->ti) - w * p->ti))
           + (IMAG(p->cw(w, p)) * (1. / (M_PI * pow(w, 2.)))
	   * (1 - cos(w * p->ti))
	   * (1. / tanh(0.5 * w * hbeta)));
    /* return (creal(p->cw(w, p)) * (1. / (M_PI * pow(w, 2.))) */
    /* 	   * (sin(w * p->ti) - w * p->ti)) */
    /*        + (cimag(p->cw(w, p)) * (1. / (M_PI * pow(w, 2.))) */
	   /* * (1 - cos(w * p->ti)) */
	   /* * (1. / tanh(0.5 * w * hbeta))); */
}

double
reorg_int(double w, void* params)
{
    Parameters *p = (Parameters *) params;
    return p->cw(w, p) * (1. / (M_PI * w));
}

COMPLEX
At(double w0, double re, double im, double t,
   double gamma)
{
    /* double cmtofs = (200.0 * M_PI * CMS * 1E-15); */
    COMPLEX exponent = - re - (0.5 * gamma * t) - II * (w0 * t + im);
    return EXP(exponent);
}

COMPLEX
Ft(double w0, double re, double im, double reorg,
   double t, double gamma)
{
    /* double cmtofs = (200.0 * M_PI * CMS * 1E-15); */
    /* NB: check bracket in final two terms here. comes out as + Gamma t/2 */
    COMPLEX exponent = -II * ((w0) - (2. * reorg)) * t 
                       - (re - (II * im) - (0.5 * t * gamma));
    return EXP(exponent);
}
