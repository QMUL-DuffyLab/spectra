#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <gsl/gsl_integration.h>

/* in cm lol */
#define HBAR 2.90283E-23
#define KB 1.3806503E-23

/* contains the relevant parameters for Chl/Car spectral densities */
typedef struct {
    double s0, s1, s2, g0, g1, g2, l0, l1, l2, w1, w2, t, T;
    double (* cw)(double, void *);
} Params;

/* Chl spectral density */
double
cw_chl(double w, void* params)
{
    /* Pass a void pointer and cast it here for compatibility
     * with gsl_function when we do the quadrature */
    Params *p = (Params *) params;
    /* 7! is 5040 */
    double c1 = (p->s1 / (5040 * 2 * pow(p->w1, 4.))) * (exp(-1. * sqrt(w / p->w1)));
    double c2 = (p->s2 / (5040 * 2 * pow(p->w2, 4.))) * (exp(-1. * sqrt(w / p->w2)));
    return ((M_PI * p->s0 * pow(w, 5.)) / (p->s1 + p->s2)) * (c1 + c2);
}

/* Car spectral density */
double
cw_car(double w, void* params)
{
    Params *p = (Params *) params;
    double c1 = 2. * p->l1 * (w * p->g1 * pow(p->w1, 2.)) /
    	      ((pow(w, 2.) - pow(p->w1, 2.)) - (pow(w, 2.) * pow(p->g1, 2.)));
    double c2 = 2. * p->l2 * (w * p->g2 * pow(p->w2, 2.)) /
    	      ((pow(w, 2.) - pow(p->w2, 2.)) - (pow(w, 2.) * pow(p->g2, 2.)));
    return (c1 + c2 + 2 * p->l0 * (w * p->g0) / (pow(w, 2.) + pow(p->g0, 2.)));
}

/* the integral for g(t) includes the spectral density function,
 * which should be switchable, but we need to call a function pointer
 * of a specific form later on. Hence, add a function pointer to the params
 * struct and use that to decide which spectral density function to use. */
double
trig_re(double w, void* params)
{
    Params *p = (Params *) params;
    return p->cw(w, p) * (1. / (M_PI * pow(w, 2.))) * (1 - cos(w * p->t))
	   * (1. / tanh((HBAR * w) / (2. * KB * p->T)));
}

double
trig_im(double w, void* params)
{
    Params *p = (Params *) params;
    return p->cw(w, p) * (1. / (M_PI * pow(w, 2.))) * (sin(w * p->t) - w * p->t);
}

/* double complex *gn(Params p, double t, double T) */
/* { */
/*     double re_res, re_err, im_res, im_err; */
/*     double complex result[2]; */
/*     p.t = t; */
/*     p.T = T; */
/*     gsl_function F; */
/*     /1* possible issue here with the function pointer *1/ */
/*     F.function = &trig_re; */
/*     F.params = &p; */
/*     gsl_integration_workspace * work = gsl_integration_workspace_alloc(1000); */
/*     /1* this will eventually do the integral *1/ */
/*     gsl_integration_qagiu(&F, 0, 1, 1e-7, 1000, work, &re_res, &re_err); */
/*     F.function = &trig_im; */
/*     gsl_integration_qagiu(&F, 0, 1, 1e-7, 1000, work, &im_res, &im_err); */
/*     result[0] = re_res + I * im_res; */
/*     result[1] = re_err + im_err; */

/*     return result; */
/* } */

double
At(double w0, double re, double im, Params p)
{
    /* Params *p = (Params *) params; */
    double complex exponent = -I * (w0 * p.t) - re - (I * im);
    return exp(exponent);
}

int
main(int argc, char** argv)
{
    Params p;
    /* can write stuff to read this in */
    p.s0 = 0.5;
    p.s1 = 0.8;
    p.s2 = 0.5;
    p.w1 = 0.56;
    p.w2 = 1.94;
    p.t = 1.0;
    p.T = 3.0;
    double (*cw)(double, void *);
    cw = &cw_chl;
    p.cw = cw;

    double re_res, re_err, im_res, im_err;
    gsl_integration_workspace * work = gsl_integration_workspace_alloc(1000);
    gsl_function F;

    for (int i = 10; i < 1000; i++) {

	/* SO: 1E-15 means that each step is a femtosecond,
	 * and the 2 Pi c * 100 gives us cm, which is
	 * what we need for the rest of the functions. */
	double cmtime = ((double) i) * 2. * M_PI * 3E8 * 100 * 1E-15;
	p.t = cmtime;
	F.function = &trig_re;
	F.params = &p;
	gsl_integration_qagiu(&F, 0., 1e-4, 1e-7, 1000, work, &re_res, &re_err);
	F.function = &trig_im;
	gsl_integration_qagiu(&F, 0., 1e-4, 1e-7, 1000, work, &im_res, &im_err);

	double w0 = 1.0;
	double Ati = At(w0, re_res, im_res, p);
	fprintf(stdout, "t = %8.5f. result: "
		"(%10.6f + %10.6fi) +- (%10.6f + %10.6fi). At = %10.6f. iterations: %i\n",
		cmtime, re_res, im_res, re_err, im_err, Ati, work->size);
    }
    gsl_integration_workspace_free(work);
    exit(EXIT_SUCCESS);
}
