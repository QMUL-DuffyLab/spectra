#include "input.h"
#include "steady_state.h"
#include <fftw3.h>
#include <stdio.h>

unsigned short int
pop_converge(double *y, double *yprev, unsigned int N, double thresh)
{
  unsigned short int result = 0;
  double *diff = calloc(N, sizeof(double));
  unsigned int conv = 0;
  for (unsigned int i = 0; i < N; i++) {
    diff[i] = fabs(yprev[i] - y[i]);
    if (diff[i] < thresh) {
      conv++;
    }
  }
  if (conv == N) {
    result = 1;
  }
  return result;
}

int
main(int argc, char** argv)
{
  FILE *fp;
  unsigned int tau, i, j;
  char *line, **lineshape_files;
  double kd;
  double complex *ex, **gi_array;
  double *eigvals, *gamma, *rates, *musq, *lambda, *integral,
         *chiw_ints, **wij, **kij, **Jij, **mu, **eig, **chiw;
  Parameters *line_params;
  fftw_complex *out, *in;
  fftw_plan plan;

  fprintf(stdout, "== Spectra calculation code ==\n");

  if (argc != 3) {
    fprintf(stdout, "Wrong number of arguments - should be 2: "
        "Input file and list of lineshape defs\n");
    exit(EXIT_FAILURE);
  }

  tau = 2000; /* again probably shouldn't hardcode this but oh well */
  /* double (*chi_p)[tau]; /1* pointer to chi(w)[i] for each exciton *1/ */
  double *chi_p;

  Input *p = read_input_file(argv[1]);

  /* malloc 1d stuff, read them in */
  eigvals = calloc(p->N, sizeof(double));
  gamma = calloc(p->N, sizeof(double));
  rates = calloc(p->N, sizeof(double));
  musq = calloc(p->N, sizeof(double));
  lambda = calloc(p->N, sizeof(double));
  line_params = malloc(p->N * sizeof(Parameters));
  chiw_ints = calloc(p->N, sizeof(double));
  gamma = read(p->gamma_file, p->N);
  lambda = read(p->lambda_file, p->N);
  eigvals = read(p->eigvals_file, p->N);
  line = malloc(200 * sizeof(char));

  /* malloc 2d stuff */
  lineshape_files = malloc(p->N * sizeof(char*));
  mu = calloc(p->N, sizeof(double*));
  gi_array = calloc(p->N, sizeof(double complex*));
  eig = calloc(p->N, sizeof(double*));
  wij = calloc(p->N, sizeof(double*));
  kij = calloc(p->N, sizeof(double*));
  Jij = calloc(p->N, sizeof(double*));
  chiw = calloc(p->N, sizeof(double*));
  for (i = 0; i < p->N; i++) {
    lineshape_files[i] = malloc(200 * sizeof(char));
    gi_array[i] = calloc(tau, sizeof(double complex));
    mu[i] = calloc(3, sizeof(double));
    eig[i] = calloc(p->N, sizeof(double));
    wij[i] = calloc(p->N, sizeof(double));
    kij[i] = calloc(p->N, sizeof(double));
    Jij[i] = calloc(p->N, sizeof(double));
    chiw[i] = calloc(tau, sizeof(double));

  }

  gi_array = read_gi(p->gi_files, p->N, tau);
  eig = read_eigvecs(p->eigvecs_file, p->N);
  mu = read_mu(p->mu_file, p->N);

  fp = fopen(argv[2], "r"); /* read in list of lineshape files here */
  for (i = 0; i < p->N; i++) {
    /* now load in the parameters for each ligand */
    fgets(line, 200, fp);
    line[strcspn(line, "\n")] = 0;
    strcpy(lineshape_files[i], line);
    line_params[i] = get_parameters(lineshape_files[i]);
    if (line_params[i].ligand == 0) {
      line_params[i].cw = &cw_car;
    } else if (line_params[i].ligand == 1) {
      /* line_params[i].cw = &cw_chl; */
      line_params[i].cw = &cw_odo;
    } else if (line_params[i].ligand == 2) {
      line_params[i].cw = &cw_odo;
    } else {
      fprintf(stdout, "Ligand number of params %d is %d."
          "No idea what's happened here.\n",
          i, line_params[i].ligand);
    }

    for (j = 0; j < p->N; j++) {
      /* this might just be eigval[i] - eigval[j] */
      /* the equations aren't balanced! need to include rates away from i etc. */
      wij[i][j] = ((eigvals[i] - lambda[i]) - (eigvals[j] - lambda[j]));
    }
  }
  int cl = fclose(fp);
  if (cl != 0) {
      fprintf(stdout, "Failed to close list of lineshape files %d.\n", cl);
      exit(EXIT_FAILURE);
  }

  kij = rate_calc(p->N, eig, eigvals, wij, line_params);
  rates = relaxation_rates(p->N, gamma, kij);

  integral = calloc(tau, sizeof(double));

  in =  (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * tau);
  out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * tau);
  plan = fftw_plan_dft_1d(tau, 
  	 in, out, FFTW_BACKWARD, FFTW_ESTIMATE);

  ex = calloc(tau, sizeof(double complex));
  for (i = 0; i < p->N; i++) {
    musq[i] = pow(mu[i][0], 2.) + pow(mu[i][1], 2.) + pow(mu[i][2], 2.);

    for (unsigned int j = 0; j < tau; j++) {
      in[j] = At(eigvals[i], creal(gi_array[i][j]), cimag(gi_array[i][j]),
                 (double)j * TOFS, line_params[i].l1, line_params[i].l2,
                 1. / gamma[i]);
    }

    fftw_execute(plan); 
    for (unsigned int j = 0; j < tau; j++) {
      chiw[i][j] = creal(out[j]) * musq[i] * 2.0;
      integral[j] += creal(out[j]) * musq[i] * 2.0;
    }

    chi_p = chiw[i];
    chiw_ints[i] = trapezoid(chi_p, tau);

  }

  fprintf(stdout, "\nWriting A(w) file\n");
  fp = fopen(p->aw_file, "w");
  for (i = 0; i < tau; i++) {
    /* unpack the ordering used by FFTW */
    kd = i * 2. * M_PI / (tau);
    fprintf(fp, "%18.10f %18.10f\n", kd / TOFS, 
        creal(integral[i]) * TOFS * (1./ sqrt(tau)) * 6.4);
  }
  cl = fclose(fp);
  if (cl != 0) {
      fprintf(stdout, "Failed to close A(w) output file %d.\n", cl);
      exit(EXIT_FAILURE);
  }


  /* one with the highest oscillator strength gets excited? */
  unsigned int max = 0;
  double musq_max = 0.0;
  fprintf(stdout, "\n----------------------------------\n"
      "Osc. strengths and chi(w) integral\n"
      "----------------------------------\n\n");
  fprintf(stdout, "Pigment        |μ^2|      ∫χ_i(w)\n");
  for (i = 0; i < p->N; i++) {
    fprintf(stdout, "%7d %10.6e %10.6e\n", i + 601, musq[i], chiw_ints[i]);
    if (musq[i] > musq_max) {
      max = i;
      musq_max = musq[i];
    }
  }

  ode_params odep;
  odep.N = p->N;
  odep.kij = kij;
  odep.rates = rates;
  odep.chiw = chiw_ints;
  double *f = calloc(p->N, sizeof(double));
  double *y = calloc(p->N, sizeof(double));
  /* check convergence */
  double *yprev = calloc(p->N, sizeof(double));
  double *boltz = calloc(p->N, sizeof(double));

  boltz = bcs(p->N, eigvals);
  fprintf(stdout, "\n-----------------\nBOLTZMANN WEIGHTS\n"
                  "-----------------\n\n");
  fprintf(stdout, "Pigment          p_i    |μ^2|*p_i\n");
  for (i = 0; i < p->N; i++) {
    fprintf(stdout, "%7d %10.6e %10.6e\n", i + 601, boltz[i],
        boltz[i] * musq[i]);
    /* possible initial values for transient absorption? */
    /* if (i == max) { */
    /*   y[i] = 1.0; */
    /* } else { */
    /*   y[i] = 0.0; */
    /* } */
  }
  fprintf(stdout, "-----------------\n");
  double xtest = 0.0;
  void *params = &odep;

  fprintf(stdout, "\n-------------------------\nSTEADY STATE FLUORESCENCE\n"
                  "-------------------------\n\n");

  const gsl_multiroot_fdfsolver_type *T;
  gsl_multiroot_fdfsolver *s;
  gsl_vector *x = gsl_vector_alloc(p->N);
  /* initial guesses for populations */
  fprintf(stdout, "Initial population guess:\n");
  for (i = 0; i < p->N; i++){
    gsl_vector_set(x, i, boltz[i] * musq[i]);
    fprintf(stdout, "%d %8.6f ", i, boltz[i] * musq[i]);
  }
  fprintf(stdout, "\n");

  gsl_multiroot_function_fdf FDF;
  FDF.f = &pop_steady_f;
  FDF.df = &pop_steady_df;
  FDF.fdf = &pop_steady_fdf;
  FDF.n = p->N;
  FDF.params = params;

  T = gsl_multiroot_fdfsolver_hybridsj;
  s = gsl_multiroot_fdfsolver_alloc(T, p->N);
  gsl_multiroot_fdfsolver_set(s, &FDF, x);

  int status;
  unsigned int iter = 0;

  do
    {
      iter++;
      status = gsl_multiroot_fdfsolver_iterate (s);

      fprintf(stdout, "iter %d: ", iter);
      for (i = 0; i < p->N; i++) {
        fprintf(stdout, "x(%2d) = %8.6f ", i, gsl_vector_get(s->x, i));
      }
      fprintf(stdout, "\n");

      if (status)   /* check if solver is stuck */
        break;

      status =
        gsl_multiroot_test_residual (s->f, 1e-7);
    }
  while (status == GSL_CONTINUE && iter < 1000);

  printf ("status = %s\n", gsl_strerror (status));

  gsl_multiroot_fdfsolver_free (s);
  gsl_vector_free (x);
  /* return 0; */

  Jij = jacmat(odep);

  int ode_success = odefunc(xtest, y, f, params);
  if (ode_success != GSL_SUCCESS) {
    fprintf(stdout, "ode_success failed!\n");
  }

  /* do ODE solving
   * this doesn't work yet */
  fprintf(stdout, "\n----------------------------------\n"
      "ODE solver for exciton populations\n"
      "----------------------------------\n\n");

  gsl_odeiv2_system sys = {odefunc, jacobian, p->N, params};

  double thresh = 1e-6;
  gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new(
      &sys, 
      gsl_odeiv2_step_bsimp,
      /* gsl_odeiv2_step_rkf45, */
      thresh,
      thresh,
      0.0);

  double t1 = 0.0; double ti = 0.0; double dt = 1.0;
  status = 0;
  for (i = 0; i < tau; i++) {
    ti = (i * dt);
    for (j = 0; j < p->N; j++) {
      yprev[j] = y[j];
    }

    /* fprintf(stdout, "iteration %d: ", i); */
    status = gsl_odeiv2_driver_apply(d, &t1, ti, y);
    fprintf(stdout, "ti = %6.3f, ", ti);
    for (j = 0; j < p->N; j++) {
      fprintf(stdout, "%8.6f ", y[j]);
    }
    fprintf(stdout, "\n");
    if (status != GSL_SUCCESS) {
      fprintf(stdout, "error: return value = %d\n", status);
      exit(EXIT_FAILURE);
    }
    if (i > 0) {
      unsigned short int converged = pop_converge(y, yprev, p->N, thresh);
      if (converged) {
        fprintf(stdout, "all populations converged to within %f; "
            "ending integration at interation %d\n", thresh, i);
        break;
      }
    }

  }
  gsl_odeiv2_driver_free(d);

  double *ynorm = calloc(p->N, sizeof(double));
  double ysum = 0.0;
  fprintf(stdout, "\n----------------------\nNormalised populations\n"
      "----------------------\n\n");
  for (i = 0; i < p->N; i++) {
    ysum += y[i];
  }
  for (i = 0; i < p->N; i++) {
    ynorm[i] = y[i] / ysum;
    fprintf(stdout, "%8.6f ", ynorm[i]);
  }
  fprintf(stdout, "\n----------------------\n");

  free(line); free(lineshape_files); free(ex); free(integral);
  free(gi_array); free(eigvals); free(gamma); free(lambda); free(mu);
  free(eig); free(wij); free(kij); free(p); free(line_params);
  free(in); free(out); free(y); free(f); free(boltz); free(yprev);
  exit(EXIT_SUCCESS);
}
