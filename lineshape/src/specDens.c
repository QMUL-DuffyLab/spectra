#include <time.h>
#include <stdio.h>
#include <math.h>
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
    double (* re)(double, void *);
    double reorg_res, reorg_err, re_res, re_err, im_res, im_err;
    gsl_function gsl_reorg, gsl_re, gsl_im;

    if (argc != 3) {
	fprintf(stdout, "Wrong number of arguments. Try again\n");
	exit(EXIT_FAILURE);
    }

    /* TESTING */
    gsl_set_error_handler_off();

    Protocol pr = get_protocol(argv[1]);
    p = get_parameters(argv[2], pr.chl_ansatz);
    /* this is ugly but need T in parameters */
    p.T = pr.T;

    times = malloc(pr.ns * sizeof(double));
    Atv   = malloc(pr.ns * sizeof(double complex));
    Ftv   = malloc(pr.ns * sizeof(double complex));

    /* choose spectral density based on the ansatz */
    p.cw = choose_ansatz(p.ans);
    p.cn = &c_n;

    gsl_integration_workspace * work = gsl_integration_workspace_alloc(1000);

    /* reorganisation energy */
    /* w0 is the offset to put the 0-0 line at 0;
     * it's set to zero in get_parameters so 
     * we're fine to just add to it in the loop */
    p.offset = get_offset(p);
    fprintf(stdout, "Anomalous phase shift = %12.8e\n", p.offset);
    FILE *fp = fopen(p.offset_file, "w");
    fprintf(fp, "%18.10f", p.offset);
    fclose(fp);

    re = &reorg_int;
    gsl_reorg.function = re;
    gsl_reorg.params = &p;
    gsl_integration_qagiu(&gsl_reorg, 0., 1e-4, 1e-7, 1000,
			  work, &reorg_res, &reorg_err);
    fprintf(stdout, "Reorganisation energy lambda = %12.8f\n",reorg_res);
    reorg_res -= p.offset;
    fp = fopen(p.lambda_file, "w");
    /* p.nu comes from mancal */
    fprintf(fp, "%18.10f %18.10f\n", reorg_res * p.nu, reorg_err);
    fclose(fp);

    gsl_re.function = &trig_re;
    gsl_im.function = &trig_im;

    /* SO: 1E-15 means that each step is a femtosecond,
     * and the 2 Pi C * 100 gives us cm, which is
     * what we need for the rest of the functions. */
    double pf = 2.* M_PI * CMS * 100. * 1E-15;
    double pf_norm = pf * (1. / sqrt(pr.ns));

    fp = fopen(p.gt_file, "w");

    for (unsigned long i = 0; i < pr.ns; i++) {

	double cmtime = ((double) i) * pf;
	times[i] = cmtime;

        p.ti = cmtime;
	gsl_re.params = &p;
	gsl_im.params = &p;

	double small_t, small_t_err;
	gsl_integration_qags(&gsl_re, 0., 1e-6, 1e-4, 1e-7, 1000,
			      work, &small_t, &small_t_err);
	gsl_integration_qagiu(&gsl_re, 1e-5, 1e-8, 1e-8, 1000,
			      work, &re_res, &re_err);
	gsl_integration_qagiu(&gsl_im, 0., 1e-8, 1e-8, 1000,
			      work, &im_res, &im_err);

	re_res = p.nu * (re_res + small_t);
	re_err = p.nu * (re_err + small_t_err);
	im_res = p.nu * (im_res + (p.offset * p.ti));

	fprintf(fp, "%18.10f %18.10f %18.10f\n",
		(float) i, re_res, im_res);

	Atv[i] = At(0.0, re_res, im_res, cmtime, 0.0);
	Ftv[i] = Ft(0.0, re_res, im_res, reorg_res,
	            cmtime, 0.);
    }
    fclose(fp);

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

    fp = fopen(p.aw_file, "w");

    for (unsigned long i = 0; i < pr.ns; i++) {
    	/* unpack the ordering used by FFTW */
    	double k = i * 2. * M_PI / (pr.ns);
    	double freq = fmod(k, M_PI) - (M_PI * floor(k / M_PI));
	/* fprintf(fp, " %18.10f %18.10f\n", */
	/* 	freq / pf, out[i][0] * pf_norm * 6.4); */
	fprintf(fp, " %18.10f %18.10f\n",
		freq / pf, out[i][0]);
    }

    fclose(fp);
    fprintf(stdout, "FFT on A(t) performed.\n");

    for (unsigned long i = 0; i < pr.ns; i++) {
	in[i][0] = creal(Ftv[i]);
	in[i][1] = cimag(Ftv[i]);
    }
    fftw_execute(plan); 

    fp = fopen(p.fw_file, "w");
    /* FILE *gp = fopen("car_cw.dat", "w"); */

    for (unsigned long i = 0; i < pr.ns; i++) {
    	double k = i * 2. * M_PI / (pr.ns);
    	double freq = fmod(k, M_PI) - (M_PI * floor(k / M_PI));
    	/* the 6.4 here is from an N in the python code */
	fprintf(fp, " %18.10f %18.10f\n",
		freq / pf, out[i][0]* pf_norm * 6.4);
	/* fprintf(gp, " %18.10f %18.10f\n", */
	/* 	freq / pf, p.cw(freq / pf, &p)); */
    }

    fclose(fp);
    /* fclose(gp); */
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
