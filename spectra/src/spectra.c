#include "input.h"
#include "steady_state.h"
#include <fftw3.h>
#include <stdio.h>
#include <gsl/gsl_eigen.h>

int
main(int argc, char** argv)
{
  FILE *fp;
  unsigned int tau, i, j;
  char *line, **lineshape_files;
  double kd;
  double complex **gi_array;
  double *eigvals, *gamma, *rates, *musq, *chi_p,
         *lambda, *integral, *chiw_ints, *ww,
         **wij, **kij, **Jij, **mu, **eig, **chiw, **pump;
  Parameters *line_params;
  fftw_complex *out, *in;
  fftw_plan plan;

  fprintf(stdout, "== Spectra calculation code ==\n");

  if (argc != 3) {
    fprintf(stdout, "Wrong number of arguments - should be 2: "
        "Input file and list of lineshape defs\n");
    exit(EXIT_FAILURE);
  }

  tau = 2048; /* again probably shouldn't hardcode this but oh well */

  /* specifies form of incident light for source term in P_i eqns */
  pulse pump_properties = { .type=0, .centre=15000., .width=300. };
  /* form of P_i(0) - check steady_state.h for details */
  ss_init population_guess = MAX;

  Input *p = read_input_file(argv[1]);

  /* malloc 1d stuff, read them in */
  eigvals = calloc(p->N, sizeof(double));
  gamma = calloc(p->N, sizeof(double));
  rates = calloc(p->N, sizeof(double));
  musq = calloc(p->N, sizeof(double));
  lambda = calloc(p->N, sizeof(double));
  line_params = malloc(p->N * sizeof(Parameters));
  chiw_ints = calloc(p->N, sizeof(double));
  ww = calloc(tau, sizeof(double));
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
  pump = calloc(p->N, sizeof(double*));
  for (i = 0; i < p->N; i++) {
    lineshape_files[i] = malloc(200 * sizeof(char));
    gi_array[i] = calloc(tau, sizeof(double complex));
    mu[i] = calloc(3, sizeof(double));
    eig[i] = calloc(p->N, sizeof(double));
    wij[i] = calloc(p->N, sizeof(double));
    kij[i] = calloc(p->N, sizeof(double));
    Jij[i] = calloc(p->N, sizeof(double));
    chiw[i] = calloc(tau, sizeof(double));
    pump[i] = calloc(tau, sizeof(double));
  }

  gi_array = read_gi(p->gi_files, p->N, tau);
  eig = read_eigvecs(p->eigvecs_file, p->N);
  mu = read_mu(p->mu_file, p->N);
  ww = incident(pump_properties, tau);

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
    line_params[i].cn = &c_n;

    for (j = 0; j < p->N; j++) {
      wij[i][j] = ((eigvals[i] - lambda[i]) - (eigvals[j] - lambda[j]));
    }
  }
  int cl = fclose(fp);
  if (cl != 0) {
      fprintf(stdout, "Failed to close list of lineshape files %d.\n", cl);
      exit(EXIT_FAILURE);
  }

  kij = rate_calc(p->N, eig, wij, line_params);
  rates = relaxation_rates(p->N, gamma);

  integral = calloc(tau, sizeof(double));

  in =  (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * tau);
  out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * tau);
  plan = fftw_plan_dft_1d(tau, 
  	 in, out, FFTW_BACKWARD, FFTW_ESTIMATE);

  for (i = 0; i < p->N; i++) {
    musq[i] = pow(mu[i][0], 2.) + pow(mu[i][1], 2.) + pow(mu[i][2], 2.);

    for (unsigned int j = 0; j < tau; j++) {
      in[j] = At(eigvals[i], creal(gi_array[i][j]), cimag(gi_array[i][j]),
                 (double)j * TOFS,
                 1. / (1000 * gamma[i]));
    }

    fftw_execute(plan); 
    for (unsigned int j = 0; j < tau; j++) {
      chiw[i][j] = creal(out[j]) * musq[i];
      pump[i][j] = chiw[i][j] * ww[j];
      integral[j] += creal(out[j]) * musq[i];
    }

    chi_p = pump[i];
    chiw_ints[i] = trapezoid(chi_p, tau);

  }

  fprintf(stdout, "\nWriting A(w) file\n");
  fp = fopen(p->aw_file, "w");
  for (i = 0; i < tau; i++) {
    /* unpack the ordering used by FFTW */
    kd = i * 2. * M_PI / (tau);
    fprintf(fp, "%18.10f %18.10f\n", kd / TOFS, 
        kd * creal(integral[i]) * (1./ sqrt(tau)) * 6.4);
  }
  cl = fclose(fp);
  if (cl != 0) {
      fprintf(stdout, "Failed to close A(w) output file %d.\n", cl);
      exit(EXIT_FAILURE);
  }

  fprintf(stdout, "\nWriting χ_i(w) files\n");
  /* replace g_ with chi_ in gi filenames so we can print out the 
   * individual parts of the spectra. strings in C are so annoying.
   * NB: assumes 200 is enough bytes for gi filename + a few chars! */
  char fn[200]; char *pch;
  for (i = 0; i < p->N; i++) {
    strcpy(fn, p->gi_files[i]);
    if((pch = strstr(fn, "g_")) == NULL) {
      fprintf(stdout, "Cannot find string 'g_' in "
          "g_i filename no. %d, needed to print out chi_i. "
          "Something got renamed accidentally?\n", i);
      break;
    }
    /* + 4 bc strlen("chi_") = 4; + 2 bc strlen("g_") = 2.
     * strlen(pch + 2) is the length of the string left after "g_",
     * so we move the end of the string (4 - 2) bytes along */
    memmove(pch + 4, pch + 2, strlen(pch + 2) + 1);
    memcpy(pch, "chi_", 4);
    fp = fopen(fn, "w");
    for (j = 0; j < tau; j++) {
      /* unpack the ordering used by FFTW */
      kd = j * 2. * M_PI / (tau);
      fprintf(fp, "%18.10f %18.10f\n", kd / TOFS, 
          kd * creal(chiw[i][j]) * (1./ sqrt(tau)) * 6.4);
    }
    cl = fclose(fp);
    if (cl != 0) {
        fprintf(stdout, "Failed to close χ_i(w) "
            "output file no. %d, error no. %d.\n", i, cl);
        exit(EXIT_FAILURE);
    }
  }

  /* one with the highest oscillator strength gets excited? */
  unsigned int max = 0;
  double musq_max = musq[0];
  double musq_sum = 0.;
  fprintf(stdout, "\n----------------------------------\n"
      "OSC. STRENGTHS AND CHI(W) INTEGRAL\n"
      "----------------------------------\n\n");
  fprintf(stdout, "Pigment        |μ^2|      ∫χ_i(w)\n");
  for (i = 0; i < p->N; i++) {
    fprintf(stdout, "%7d %10.6e %10.6e\n", i + 601, musq[i], chiw_ints[i]);
    if (musq[i] > musq_max) {
      max = i;
      musq_max = musq[i];
      musq_sum += musq[i];
    }
  }

  ode_params odep;
  odep.N = p->N;
  odep.kij = kij;
  odep.rates = rates;
  odep.chiw = chiw_ints;
  odep.Tij = transfer_matrix(p->N, rates, kij);
  double *f = calloc(p->N, sizeof(double));
  double *y = calloc(p->N, sizeof(double));
  /* check convergence */
  double *yprev = calloc(p->N, sizeof(double));
  double *boltz = calloc(p->N, sizeof(double));

  boltz = bcs(p->N, eigvals);
  fprintf(stdout, "\n-----------------\nBOLTZMANN WEIGHTS\n"
                  "-----------------\n\n");
  fprintf(stdout, "Pigment    p_i   |μ^2|*p_i\n");
  for (i = 0; i < p->N; i++) {
    fprintf(stdout, "%7d %8.6f %8.6f\n", i + 601, boltz[i],
        boltz[i] * musq[i]);
    /* possible initial values for transient absorption? */
    if (i == max) {
      y[i] = 1.0;
    } else {
      y[i] = 0.0;
    }
  }
  fprintf(stdout, "\n");
  double xtest = 0.0;
  void *params = &odep;

  fprintf(stdout, "\n-------------------------\nSTEADY STATE FLUORESCENCE\n"
                  "-------------------------\n\n");

  const gsl_multiroot_fdfsolver_type *T;
  gsl_multiroot_fdfsolver *s;
  gsl_vector *x = guess(population_guess, boltz, musq, max, p->N);
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

  double *p_i_equib = calloc(p->N, sizeof(double));
  double p_i_sum = 0.0;
  double boltz_sum = 0.0;
  for (i = 0; i < p->N; i++) {
    /* amazingly there does not seem to be a GSL function for this */
    p_i_sum += gsl_vector_get(s->x, i);
    boltz_sum += boltz[i] * musq[i];
  }
  if (p_i_sum <= 1E-10) {
    fprintf(stdout, "Sum of steady-state populations is zero!!!\n");
    /* this stops it from e.g. normalising a vector (0, 1e-23)
     * to (0, 1) and making it look normal */
    p_i_sum = 1.;
  } else {
    fprintf(stdout, "Σ_i p_i = %8.6f\n", p_i_sum);
  }

  /* write out the normalised steady-state populations */
  strcpy(fn, p->pop_file);
  pch = strstr(fn, "populations");
  if((pch = strstr(fn, "populations")) == NULL) {
    fprintf(stdout, "Cannot find string 'populations' in "
        "pop_file needed to print out steady-state populations. "
        "Something got renamed accidentally?\n");
    exit(EXIT_FAILURE);
  }
  memmove(pch + 14, pch + 11, strlen(pch + 11) + 1);
  memcpy(pch, "ss_populations", 14);
  fp = fopen(fn, "w");

  fprintf(stdout, "i\t p_i^eq(raw)\t p_i^eq(norm)\t boltz*|μ^2|\n");
  fprintf(fp, "# i\t p_i^eq(raw)\t p_i^eq(norm)\t boltz*|μ^2|\n");
  for (i = 0; i < p->N; i++) {
    p_i_equib[i] = gsl_vector_get(s->x, i) / p_i_sum;
    fprintf(stdout, "%2d\t%+12.8e\t%+12.8e\t%+12.8e\n", i,
        gsl_vector_get(s->x, i), p_i_equib[i],
        (boltz[i] * musq[i]) / boltz_sum);
    fprintf(fp, "%3d\t%+12.8e\t%+12.8e\t%+12.8e\n", i,
        gsl_vector_get(s->x, i), p_i_equib[i],
        (boltz[i] * musq[i]) / boltz_sum);
  }
  cl = fclose(fp);
  if (cl != 0) {
      fprintf(stdout, "Failed to close steady state"
          " population file, error no. %d.\n", cl);
      exit(EXIT_FAILURE);
  }

  /* fluorescence spectrum */
  /* need to zero the running integral before FFT! */
  free(integral);
  integral = calloc(tau, sizeof(double));
  for (i = 0; i < p->N; i++) {
    for (unsigned int j = 0; j < tau; j++) {
      in[j] = p_i_equib[i] * Ft(eigvals[i],
              creal(gi_array[i][j]), cimag(gi_array[i][j]),
              lambda[i], (double)j * TOFS,
              1. / (1000 * gamma[i]));
    }

    fftw_execute(plan); 
    for (unsigned int j = 0; j < tau; j++) {
      chiw[i][j] = creal(out[j]) * musq[i];
      integral[j] += creal(out[j]) * musq[i];
    }

  }

  fprintf(stdout, "\nWriting F(w) file\n");
  fp = fopen(p->fw_file, "w");
  for (i = 0; i < tau; i++) {
    /* unpack the ordering used by FFTW */
    kd = i * 2. * M_PI / (tau);
    fprintf(fp, "%18.10f %18.10f\n", kd / TOFS, 
        (pow(kd, 3.) * creal(integral[i]) * (1./ sqrt(tau))) * 6.4);
  }
  cl = fclose(fp);
  if (cl != 0) {
      fprintf(stdout, "Failed to close F(w) output file %d.\n", cl);
      exit(EXIT_FAILURE);
  }

  fprintf(stdout, "\nWriting χ^{bar}_i(w) files\n");
  /* replace with chi_bar_ this time */
  for (i = 0; i < p->N; i++) {
    strcpy(fn, p->gi_files[i]);
    pch = strstr(fn, "g_");
    if((pch = strstr(fn, "g_")) == NULL) {
      fprintf(stdout, "Cannot find string 'g_' in "
          "g_i filename no. %d, needed to print out chi_i. "
          "Something got renamed accidentally?\n", i);
      break;
    }
    memmove(pch + 8, pch + 2, strlen(pch + 2) + 1);
    memcpy(pch, "chi_bar_", 8);
    fp = fopen(fn, "w");
    for (j = 0; j < tau; j++) {
      /* unpack the ordering used by FFTW */
      kd = j * 2. * M_PI / (tau);
      fprintf(fp, "%18.10f %18.10f\n", kd / TOFS, 
          kd * creal(chiw[i][j]) * (1./ sqrt(tau)) * 6.4);
    }
    cl = fclose(fp);
    if (cl != 0) {
        fprintf(stdout, "Failed to close χ_i(w) "
            "output file no. %d, error no. %d.\n", i, cl);
        exit(EXIT_FAILURE);
    }
  }

  gsl_multiroot_fdfsolver_free(s);
  gsl_vector_free(x);

  /* do ODE solving
   * this doesn't work yet */
  int ode_success = odefunc(xtest, y, f, params);
  if (ode_success != GSL_SUCCESS) {
    fprintf(stdout, "ode_success failed!\n");
  }

  fprintf(stdout, "\n----------------------------------\n"
      "ODE SOLVER FOR EXCITON POPULATIONS\n"
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
  unsigned int MAX_ITER = 20000;
  unsigned int print_pop = 0;
  status = 0;

  fprintf(stdout, "\nWriting populations\n");
  fp = fopen(p->pop_file, "w");

  for (i = 0; i < MAX_ITER; i++) {
    ti = (i * dt);
    for (j = 0; j < p->N; j++) {
      yprev[j] = y[j];
    }

    /* fprintf(stdout, "iteration %d: ", i); */
    status = gsl_odeiv2_driver_apply(d, &t1, ti, y);

    fprintf(fp, "%6.3f ", ti);
    if (print_pop) {
      fprintf(stdout, "ti = %6.3f ", ti);
    }
    for (j = 0; j < p->N; j++) {
      fprintf(fp, "%+12.8e ", y[j]);
      if (print_pop) {
        fprintf(stdout, "%+12.8e ", y[j]);
      }
    }
    fprintf(fp, "\n");
    if(print_pop) {
      fprintf(stdout, "\n");
    }

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
  fclose(fp);
  gsl_odeiv2_driver_free(d);

  fprintf(stdout, "\n-------------------\nEXCITATION LIFETIME\n"
      "-------------------\n\n");
  gsl_complex one, item;
  GSL_SET_COMPLEX(&one, 1., 0.);
  /* first assign the transfer matrix to a GSL matrix,
   * then call the relevant eigensystem function; we can
   * use the eigendecomposed bits to calculate <\tau> */
  gsl_eigen_nonsymmv_workspace *w = gsl_eigen_nonsymmv_alloc(p->N);
  gsl_matrix *Tij_gsl = array_to_gsl_matrix(p->N, p->N, odep.Tij);
  /* these are real actually */
  gsl_vector_complex *eval = gsl_vector_complex_alloc(p->N);
  gsl_matrix_complex *evec = gsl_matrix_complex_alloc(p->N, p->N);

  status = gsl_eigen_nonsymmv(Tij_gsl, eval, evec, w);
  for (i = 0; i < p->N; i++) {
    for (j = 0; j < p->N; j++) {
      fprintf(stdout, "(%+10.6e %+10.6e) ",
          GSL_REAL(gsl_matrix_complex_get(evec, i, j)),
          GSL_IMAG(gsl_matrix_complex_get(evec, i, j)));
    }
    fprintf(stdout, "\n");
  }

  /* now we need to invert the eigenvector matrix as well */
  /* there must be a more efficient way of doing this! */
  /* the matrix C^-1, the inverse of the eigenvector matrix,
   * is complex; so are the corresponding eigenvalues */
  gsl_vector_complex *eval_temp = gsl_vector_complex_alloc(p->N);
  gsl_matrix_complex *evec_inv  = gsl_matrix_complex_alloc(p->N, p->N);
  gsl_matrix_complex *intermed  = gsl_matrix_complex_alloc(p->N, p->N);
  gsl_matrix_complex *eval_m1   = gsl_matrix_complex_calloc(p->N, p->N);
  gsl_vector_complex *p_vector  = gsl_vector_complex_alloc(p->N);
  gsl_vector_complex *final_vec = gsl_vector_complex_alloc(p->N);

  for (i = 0; i < p->N; i++) {
    /* this is so dumb. idk if we can assume everything's real or not */
    GSL_SET_COMPLEX(&item,
        1. / GSL_REAL(gsl_vector_complex_get(eval, i)), 0.);
    /* set diagonals to 1/eigenvalues */
    gsl_matrix_complex_set(eval_m1, i, i, item);
    GSL_SET_COMPLEX(&item, p_i_equib[i], 0.);
    /* GSL_SET_COMPLEX(&item, musq[i] / musq_sum, 0.); */
    /* GSL_SET_COMPLEX(&item, 1./p->N, 0.); */
    gsl_vector_complex_set(p_vector, i, item);
  }

  status = gsl_eigen_nonsymmv(evec, eval_temp, evec_inv, w);

  /* now we do evec * eval^-1 * evec^-1 * P(0) */
  gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, one,
      evec, eval_m1, one, intermed);
  gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, one,
      intermed, evec_inv, one, evec);
  gsl_blas_zgemv(CblasNoTrans, one,
      evec, p_vector, one, final_vec);
  double excite = 0.;
  for (i = 0; i < p->N; i++) {
    excite += GSL_REAL(gsl_vector_complex_get(final_vec, i));
  }
  fprintf(stdout, "mean excitation lifetime = %12.8e\n", excite);

  gsl_matrix_free(Tij_gsl);
  gsl_vector_complex_free(eval);
  gsl_vector_complex_free(eval_temp);
  gsl_vector_complex_free(p_vector);
  gsl_vector_complex_free(final_vec);
  gsl_matrix_complex_free(evec);
  gsl_matrix_complex_free(evec_inv);
  gsl_matrix_complex_free(intermed);
  gsl_matrix_complex_free(eval_m1);

  double *ynorm = calloc(p->N, sizeof(double));
  double ysum = 0.0;
  fprintf(stdout, "\n----------------------\nNORMALISED POPULATIONS\n"
      "----------------------\n\n");
  for (i = 0; i < p->N; i++) {
    ysum += y[i];
  }
  for (i = 0; i < p->N; i++) {
    ynorm[i] = y[i] / ysum;
    fprintf(stdout, "%8.6f ", ynorm[i]);
  }
  fprintf(stdout, "\n----------------------\n");

  free(line); free(lineshape_files); free(integral);
  free(gi_array); free(eigvals); free(gamma); free(lambda); free(mu);
  free(eig); free(wij); free(kij); free(p); free(line_params);
  free(in); free(out); free(y); free(f); free(boltz); free(yprev);
  free(ynorm); free(ww);
  exit(EXIT_SUCCESS);
}
