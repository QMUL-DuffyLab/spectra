#include <time.h>
#include <gsl/gsl_integration.h>
#include <fftw3.h>
#include "functions.h"
#include "parameters.h"

int
main(int argc, char** argv)
{
    time_t start_time, end_time;
    time(&start_time);
    Parameters p;
    double *times;
    double complex *Atv, *Ftv;
    fftw_complex *out, *in;
    fftw_plan plan;
    double (*cw)(double, void *);
    double (*re)(double, void *);
    double reorg_res, reorg_err, re_res, re_err, im_res, im_err;
    gsl_function gsl_reorg, gsl_re, gsl_im;

    /* ALL THIS SHOULD BE READ IN */
    p = getParameters(argv[2]);
    Protocol pr = getProtocol(argv[1]);
    /* this is ugly but need T in parameters
     * might just get rid of protocol file */
    p.T = pr.T;

    times = malloc(pr.ns * sizeof(double));
    Atv   = malloc(pr.ns * sizeof(double complex));
    Ftv   = malloc(pr.ns * sizeof(double complex));

    if (p.ligand == 1) {
	cw = &cw_chl;
    } else if (p.ligand == 0) {
	cw = &cw_car;
    } else {
    	fprintf(stdout, "Ligand selection not working\n");
    	exit(EXIT_FAILURE);
    }

    cw = &cw_chl;
    p.cw = cw;

    gsl_integration_workspace * work = gsl_integration_workspace_alloc(1000);

    /* reorganisation energy */
    re = &reorg_int;
    gsl_reorg.function = re;
    gsl_reorg.params = &p;
    gsl_integration_qagiu(&gsl_reorg, 0., 1e-4, 1e-7, 1000,
			  work, &reorg_res, &reorg_err);
    gsl_re.function = &trig_re;
    gsl_im.function = &trig_im;

    for (int i = 7; i < pr.ns; i++) {

	/* SO: 1E-15 means that each step is a femtosecond,
	 * and the 2 Pi c * 100 gives us cm, which is
	 * what we need for the rest of the functions.
	 * also I start from i = 10 because the integral
	 * seems to be singular as t -> 0; look at this later */
	double cmtime = ((double) i) * 2. * M_PI * 3E8 * 100 * 1E-15;
	times[i] = cmtime;

	p.t = cmtime;
	gsl_re.params = &p;
	gsl_im.params = &p;

	double small_t, small_t_err;
	gsl_integration_qags(&gsl_re, 0., 1e-6, 1e-4, 1e-7, 1000,
			      work, &small_t, &small_t_err);
	gsl_integration_qagiu(&gsl_re, 1e-4, 1e-4, 1e-7, 1000,
			      work, &re_res, &re_err);
	gsl_integration_qagiu(&gsl_im, 0., 1e-4, 1e-7, 1000,
			      work, &im_res, &im_err);

	re_res = re_res + small_t;
	re_err = re_err + small_t_err;
	double w0 = 0.0;
	Atv[i] = At(w0, re_res, im_res, cmtime);
	Ftv[i] = Ft(w0, re_res, im_res, reorg_res, cmtime);
	/* fprintf(stdout, "t = %8.5f. result: " */
	/* 	"(%10.6f + %10.6fi) +- (%10.6f + %10.6fi)." */
	/* 	"At = %10.6f. Ft = %10.6f. iterations: %lu\n", */
	/* 	cmtime, re_res, im_res, re_err, */
		/* im_err, Atv[i], Ftv[i], work->size); */
	/* fprintf(stdout, "t = %8.5f " */
	/* 	"At = %8.5f + %8.5f\n", */
	/* 	cmtime, creal(Atv[i]), cimag(Atv[i])); */
    }

    fprintf(stdout, "Performing FFT.\n");

    /* this assignment is very ugly but trying to copy Chris's code */
    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * pr.ns);
    for (int i = 0; i < pr.ns; i++) {
	in[i][0] = creal(Atv[i]);
	in[i][1] = cimag(Atv[i]);
    }

    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * pr.ns);
    plan = fftw_plan_dft_1d(pr.ns, 
    	 in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(plan); 
    double pf = 2.* M_PI * 3E8 * 100. * 1E-15 * sqrt(pr.ns);
    for (int i = 0; i < pr.ns; i++) {
	fprintf(stdout, "FT = %8.5f\n",
		out[i][0]*pf);
    }
    fprintf(stdout, "FFT performed.\n");

    time(&end_time);
    /* this is pretty useless, i forgot it only does integer seconds */
    double time_taken = difftime(end_time, start_time);
    fprintf(stdout, "Time taken: %12.8f\n",
	    time_taken);
    gsl_integration_workspace_free(work);
    fftw_destroy_plan(plan);
    fftw_free(in); fftw_free(out); free(times); free(Atv); free(Ftv);
    exit(EXIT_SUCCESS);
}
