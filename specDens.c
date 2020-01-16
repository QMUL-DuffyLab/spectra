#include <stdio.h>
#include <math.h>
#include <time.h>
#include <complex.h>
#include <gsl/gsl_integration.h>
#include <fftw3.h>

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
    /* 7! is 5040; this is the Renger form for chlorophyll */
    double c1 = (p->s1 / (5040 * 2 * pow(p->w1, 4.))) * (exp(-1. * sqrt(w / p->w1)));
    double c2 = (p->s2 / (5040 * 2 * pow(p->w2, 4.))) * (exp(-1. * sqrt(w / p->w2)));
    return ((M_PI * p->s0 * pow(w, 5.)) / (p->s1 + p->s2)) * (c1 + c2);
}

/* Car spectral density */
double
cw_car(double w, void* params)
{
    Params *p = (Params *) params;
    /* ansatz from Kieran's paper on carotenoids */
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

double
reorg_int(double w, void* params)
{
    Params *p = (Params *) params;
    return p->cw(w, p) * (1. / (M_PI * w));
}

double
At(double w0, double re, double im, Params p)
{
    /* Params *p = (Params *) params; */
    double complex exponent = -I * (w0 * p.t) - (re + (I * im));
    return exp(exponent);
}

double
Ft(double w0, double re, double im, double reorg, Params p)
{
    double complex exponent = -I * (w0 * p.t) - (re + (I * im)) - (2. * reorg);
    return exp(exponent);
}

int
main(int argc, char** argv)
{
    time_t start_time, end_time;
    time(&start_time);
    Params p;
    double *times, *Atv, *Ftv;
    fftw_complex *in, *out;
    fftw_plan plan;

    /* ALL THIS SHOULD BE READ IN */
    int num_steps = 1000;
    p.s0 = 0.5;
    p.s1 = 0.8;
    p.s2 = 0.5;
    p.w1 = 0.56;
    p.w2 = 1.94;
    p.t = 1.0;
    p.T = 3.0;
    double (*cw)(double, void *);
    double (*re)(double, void *);
    /* note: this'll have to be done via a string 
     * switch in an input file because
     * you can't do &var of a string */
    cw = &cw_chl;
    p.cw = cw;

    double reorg_res, reorg_err, re_res, re_err, im_res, im_err;
    gsl_integration_workspace * work = gsl_integration_workspace_alloc(1000);
    gsl_function F;
    times = malloc(num_steps * sizeof(double));
    Atv   = malloc(num_steps * sizeof(double));
    Ftv   = malloc(num_steps * sizeof(double));

    /* reorganisation energy */
    re = &reorg_int;
    F.function = re;
    F.params = &p;
    gsl_integration_qagiu(&F, 0., 1e-4, 1e-7, 1000,
			  work, &reorg_res, &reorg_err);
    /* need to switch back to the spectral density function:
     * this is a very ugly way of doing this honestly */
    cw = &cw_chl;

    for (int i = 10; i < num_steps; i++) {

	/* SO: 1E-15 means that each step is a femtosecond,
	 * and the 2 Pi c * 100 gives us cm, which is
	 * what we need for the rest of the functions. */
	double cmtime = ((double) i) * 2. * M_PI * 3E8 * 100 * 1E-15;
	p.t = cmtime;
	times[i] = cmtime;
	F.function = &trig_re;
	F.params = &p;
	gsl_integration_qagiu(&F, 0., 1e-4, 1e-7, 1000, work, &re_res, &re_err);
	F.function = &trig_im;
	gsl_integration_qagiu(&F, 0., 1e-4, 1e-7, 1000, work, &im_res, &im_err);

	double w0 = 1.0;
	Atv[i] = At(w0, re_res, im_res, p);
	Ftv[i] = Ft(w0, re_res, im_res, reorg_res, p);
	fprintf(stdout, "t = %8.5f. result: "
		"(%10.6f + %10.6fi) +- (%10.6f + %10.6fi)."
		"At = %10.6f. Ft = %10.6f. iterations: %lu\n",
		cmtime, re_res, im_res, re_err,
		im_err, Atv[i], Ftv[i], work->size);
    }
    time(&end_time);
    /* this is pretty useless, i forgot it only does integer seconds */
    double time_taken = difftime(end_time, start_time);
    fprintf(stdout, "Time taken: %12.8f\n",
	    time_taken);
    gsl_integration_workspace_free(work);
    free(times);
    free(Atv);
    exit(EXIT_SUCCESS);
}
