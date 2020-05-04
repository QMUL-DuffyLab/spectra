#include "input.h"
#include "fluorescence.h"
#include <fftw3.h>
#include <stdio.h>

int
main(int argc, char** argv)
{
  FILE *fp;
  unsigned int tau, i, j;
  char *line, **lineshape_files;
  double musq, kd;
  double complex *ex, **gi_array;
  double *eigvals, *gamma, *lambda, *integral,
         **wij, **kij, **Jij, **mu, **eig;
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

  Input *p = read_input_file(argv[1]);

  /* malloc 1d stuff, read them in */
  eigvals = calloc(p->N, sizeof(double));
  gamma = calloc(p->N, sizeof(double));
  lambda = calloc(p->N, sizeof(double));
  line_params = malloc(p->N * sizeof(Parameters));
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
  for (i = 0; i < p->N; i++) {
    lineshape_files[i] = malloc(200 * sizeof(char));
    gi_array[i] = calloc(tau, sizeof(double complex));
    mu[i] = calloc(3, sizeof(double));
    eig[i] = calloc(p->N, sizeof(double));
    wij[i] = calloc(p->N, sizeof(double));
    kij[i] = calloc(p->N, sizeof(double));
    Jij[i] = calloc(p->N, sizeof(double));

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
      wij[i][j] = ((eigvals[i] - lambda[i]) - (eigvals[j] - lambda[j]));
    }
  }
  int cl = fclose(fp);
  if (cl != 0) {
      fprintf(stdout, "Failed to close list of lineshape files %d.\n", cl);
      exit(EXIT_FAILURE);
  }

  kij = rate_calc(p->N, eig, wij, line_params);
  Jij = jacobian_calc(p->N, kij, gamma);
  fprintf(stdout, "Jacobian matrix:\n\n");
  for (i = 0; i < p->N; i++) {
    for (j = 0; j < p->N; j++) {
      fprintf(stdout, "%12.6e ", Jij[i][j]);
    }
    fprintf(stdout, "\n");
  }

  integral = calloc(tau, sizeof(double));

  in =  (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * tau);
  out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * tau);
  plan = fftw_plan_dft_1d(tau, 
  	 in, out, FFTW_BACKWARD, FFTW_ESTIMATE);

  ex = calloc(tau, sizeof(double complex));
  for (i = 0; i < p->N; i++) {
    musq = pow(mu[i][0], 2.) + pow(mu[i][2], 2.) + pow(mu[i][2], 2.);

    for (unsigned int j = 0; j < tau; j++) {
      in[j] = At(eigvals[i], creal(gi_array[i][j]), cimag(gi_array[i][j]),
                 (double)j * TOFS, line_params[i].l1, line_params[i].l2,
                 1. / gamma[i]);
    }

    fftw_execute(plan); 
    for (unsigned int j = 0; j < tau; j++) {
      integral[j] += creal(out[j]) * musq * 2.0;
    }

  }

  ode_params odep;
  odep.N = p->N;
  odep.kij = kij;
  odep.gamma = gamma;
  double *f = calloc(p->N, sizeof(double));
  double *y = calloc(p->N, sizeof(double));
  /* test boundary condition */
  for (i = 0; i < p->N; i++) {
    if (i == 0) {
      y[i] = 1.0;
    } else {
      y[i] = 0.0;
    }
  }
  double xtest = 0.0;
  void *params = &odep;

  int ode_success = odefunc(xtest, y, f, params);
  if (ode_success != GSL_SUCCESS) {
    fprintf(stdout, "ode_success failed!\n");
  }

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

  free(line); free(lineshape_files); free(ex); free(integral);
  free(gi_array); free(eigvals); free(gamma); free(lambda); free(mu);
  free(eig); free(wij); free(kij); free(p); free(line_params);
  free(in); free(out);
  exit(EXIT_SUCCESS);
}
