#include <time.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
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

    if (argc != 3) {
	fprintf(stdout, "Wrong number of arguments. Try again\n");
	exit(EXIT_FAILURE);
    }

    /* TESTING */
    gsl_set_error_handler_off();

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

    p.cw = cw;

    gsl_integration_workspace * work = gsl_integration_workspace_alloc(1000);

    /* reorganisation energy */
    re = &reorg_int;
    gsl_reorg.function = re;
    gsl_reorg.params = &p;
    gsl_integration_qagiu(&gsl_reorg, 0., 1e-4, 1e-7, 1000,
			  work, &reorg_res, &reorg_err);
    fprintf(stdout, "Reorganisation energy lambda = %12.8f\n",reorg_res);

    gsl_re.function = &trig_re;
    gsl_im.function = &trig_im;

    /* SO: 1E-15 means that each step is a femtosecond,
     * and the 2 Pi C * 100 gives us cm, which is
     * what we need for the rest of the functions. */
    double pf = 2.* M_PI * CMS * 100. * 1E-15;
    double pf_norm = pf * (1. / sqrt(pr.ns));
    double w0 = 0.0;

    for (unsigned long i = 0; i < pr.ns; i++) {

	double cmtime = ((double) i) * pf;
	times[i] = cmtime;

	p.t = cmtime;
	gsl_re.params = &p;
	gsl_im.params = &p;

	double small_t, small_t_err;
	gsl_integration_qags(&gsl_re, 0., 1e-6, 1e-4, 1e-7, 1000,
			      work, &small_t, &small_t_err);
	gsl_integration_qagiu(&gsl_re, 1e-5, 1e-8, 1e-8, 1000,
			      work, &re_res, &re_err);
	gsl_integration_qagiu(&gsl_im, 0., 1e-8, 1e-8, 1000,
			      work, &im_res, &im_err);

	re_res = re_res + small_t;
	re_err = re_err + small_t_err;
	Atv[i] = At(w0, re_res, im_res, cmtime);
	Ftv[i] = Ft(w0, re_res, im_res, reorg_res, cmtime);
    }

    fprintf(stdout, "Performing FFT.\n");

    /* this assignment is very ugly but trying to copy Chris's code */
    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * pr.ns);

    for (unsigned long i = 0; i < pr.ns; i++) {
	in[i][0] = creal(Atv[i]);
	in[i][1] = cimag(Atv[i]);
    }

    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * pr.ns);
    plan = fftw_plan_dft_1d(pr.ns, 
    	 in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(plan); 

    FILE *fp = fopen(pr.aw_file, "w");

    for (unsigned long i = 0; i < pr.ns; i++) {
    	double k = i * 2. * M_PI / (pr.ns);
    	double freq = fmod(k, M_PI) - (M_PI * floor(k / M_PI));
    	/* the 6.4 here is from an N in the python code */
	fprintf(fp, " %12.8f %12.8f\n",
		freq / pf, out[i][0]* pf_norm * 6.4);
    }

    fclose(fp);
    fprintf(stdout, "FFT on A(t) performed.\n");

    for (unsigned long i = 0; i < pr.ns; i++) {
	in[i][0] = creal(Ftv[i]);
	in[i][1] = cimag(Ftv[i]);
    }
    fftw_execute(plan); 

    fp = fopen(pr.fw_file, "w");

    for (unsigned long i = 0; i < pr.ns; i++) {
    	double k = i * 2. * M_PI / (pr.ns);
    	double freq = fmod(k, M_PI) - (M_PI * floor(k / M_PI));
    	/* the 6.4 here is from an N in the python code */
	fprintf(fp, " %12.8f %12.8f\n",
		freq / pf, out[i][0]* pf_norm * 6.4);
    }

    fclose(fp);
    fprintf(stdout, "FFT on F(t) performed.\n");

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
