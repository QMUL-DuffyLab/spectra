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
    double *times, *Atv, *Ftv;
    fftw_complex *out, *in;
    fftw_plan plan;
    double (*cw)(double, void *);
    double (*re)(double, void *);
    double reorg_res, reorg_err, re_res, re_err, im_res, im_err;
    gsl_function gsl_reorg, gsl_re, gsl_im;

    /* ALL THIS SHOULD BE READ IN */
    int num_steps = 1000;
    p = getParameters(argv[1]);

    times = malloc(num_steps * sizeof(double));
    Atv   = malloc(num_steps * sizeof(double));
    Ftv   = malloc(num_steps * sizeof(double));

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

    for (int i = 10; i < num_steps; i++) {

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

	gsl_integration_qagiu(&gsl_re, 0., 1e-4, 1e-7, 1000,
			      work, &re_res, &re_err);
	gsl_integration_qagiu(&gsl_im, 0., 1e-4, 1e-7, 1000,
			      work, &im_res, &im_err);

	double w0 = 0.0;
	Atv[i] = At(w0, re_res, im_res, cmtime);
	Ftv[i] = Ft(w0, re_res, im_res, reorg_res, cmtime);
	fprintf(stdout, "t = %8.5f. result: "
		"(%10.6f + %10.6fi) +- (%10.6f + %10.6fi)."
		"At = %10.6f. Ft = %10.6f. iterations: %lu\n",
		cmtime, re_res, im_res, re_err,
		im_err, Atv[i], Ftv[i], work->size);
    }

    fprintf(stdout, "Performing FFT.\n");

    /* this assignment is very ugly but trying to copy Chris's code */
    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * num_steps);
    for (int i = 0; i < num_steps; i++) {
	in[i][0] = Atv[i];
	in[i][1] = 0.0;
    }

    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * num_steps);
    plan = fftw_plan_dft_1d(num_steps, 
    	 in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(plan); 
    double pf = 2.* M_PI * 3E8 * 100. * 1E-15 * sqrt(num_steps);
    for (int i = 0; i < num_steps; i++) {
	fprintf(stdout, "FFT element = (%8.5f + %8.5f)\n",
		out[i][0]*pf, out[i][1]*pf);
    }
    fprintf(stdout, "FFT performed.\n");

    time(&end_time);
    /* this is pretty useless, i forgot it only does integer seconds */
    double time_taken = difftime(end_time, start_time);
    fprintf(stdout, "Time taken: %12.8f\n",
	    time_taken);
    gsl_integration_workspace_free(work);
    fftw_destroy_plan(plan);
    fftw_free(out);
    free(times);
    free(Atv);
    free(Ftv);
    exit(EXIT_SUCCESS);
}
