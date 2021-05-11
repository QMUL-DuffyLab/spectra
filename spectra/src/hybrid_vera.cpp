#include "hybrid_vera.h"

std::vector<double>
k_i_xa_hybrid(std::vector<VERA> x, unsigned n_chl, unsigned n_car, unsigned tau,
               double **eig, double *eigvals,
               double **Jij, double **normed_ai, double **normed_fi,
               pulse v_abs, double beta, double **car_decays)
{
  std::vector<size_t> vib_totals(n_car, 0.);
  for (unsigned i = 0; i < n_car; i++) {
    vib_totals[i] = pow(x[i].n_vib + 1, x[i].n_normal);
  }
  double ji = 0., ji_work = 0.;
  double *abs = (double *)calloc(tau, sizeof(double));
  double *fi_ad  = (double *)calloc(tau, sizeof(double));
  double *ai_fd  = (double *)calloc(tau, sizeof(double));
  std::vector<double> k_i_xa;
  unsigned short print_ji = 0;
  unsigned short print_delta_fc = 0;
  std::vector<std::vector<size_t>> pop_extents;
  for (unsigned i = 0; i < n_car; i++) {
    pop_extents.push_back(x[i].get_pop_extents()); 
  }
  std::vector<std::vector<double>> car_rates;
  for (unsigned i = 0; i < n_car; i++) {
    car_rates.push_back(x[i].intra_rates()); 
  }
  std::vector<size_t> n_s_totals(n_car, 0);
  for (unsigned i = 0; i < n_car; i++) {
    n_s_totals[i] = vib_totals[i] * x[i].n_elec;
  }
  bool print_decay_details = false;
  bool print_details = false;
  bool output_lineshapes = false;

  for (unsigned chl_index = 0; chl_index < n_chl; chl_index++) {
    for (unsigned carotenoid = n_chl;
        carotenoid < n_chl + n_car; carotenoid++) {
      /* assumes the carotenoids are at the end of the arrays;
       * fortran code puts them there so shouldn't be an issue */

      ji = 0.;
      /* for (unsigned n = 0; n < n_chl; n++) { */
        for (unsigned m = 0; m < n_chl; m++) {
          ji += eig[m][chl_index]
              /* * pow(eig[n][chl_index], 2.) */
              /* * Jij[m][carotenoid] * Jij[n][carotenoid]; */
              * Jij[m][carotenoid];
        }
      /* } */

      if (print_ji) {
        fprintf(stdout, "%5u %12.8e\n", chl_index, ji);
      }

      /* outer loop runs from vib_total to 2 * vib_total because
       * we're looping over the first excited electronic state (S1).
       * the upward rate is just the rate from the true g/s to this
       * state; then the inner loop is over the whole set of ground
       * states, to get the effective decay rate. */
      size_t car_index = carotenoid - n_chl;
      for (unsigned ii = vib_totals[car_index];
          ii < 2 * vib_totals[car_index]; ii++) {
        std::vector<size_t> a = ind2sub(ii, pop_extents[car_index]);
        /* we're only considering the transition from the true g/s
         * up to (1, {a}) */
        double e_xa = x[car_index].get_w_elec(a[0]);
        double chl_car = 0.;
        double car_chl = 0.;
        double fc = 1.;

        for (unsigned alpha = 0; alpha < x[car_index].n_normal; alpha++) {
          e_xa += x[car_index].get_w_normal(alpha) * a[alpha + 1];
          /* fc_sq *= pow(x.get_fc({0, 1, alpha, */
          /*       a[alpha + 1], 0}), 2.); */
          fc *= x[car_index].get_fc({0, 1, alpha, a[alpha + 1], 0});
        }
        ji_work = ji * fc;

        /* e_xa is also delta_xy_ba in this case since y_b = {0} */
        /* first element is the S0-S1 width which is all we need */
        v_abs.width = x[car_index].get_widths(0);
        /* now we have to calculate the rate between every delta_xy_ba */
        v_abs.centre = e_xa;
        abs = incident(v_abs, tau);

        for (unsigned step = 0; step < tau; step++) {
          fi_ad[step] = normed_fi[chl_index][step] * abs[step];
          ai_fd[step] = normed_ai[chl_index][step] * abs[step];
        }
        
        chl_car += pow(ji_work, 2.)
                * trapezoid(fi_ad, 2. * M_PI / (TOFS * tau), tau);
        car_chl += pow(ji_work, 2.)
                * trapezoid(ai_fd, 2. * M_PI / (TOFS * tau), tau);

        if (print_delta_fc) {
          fprintf(stdout, "%2u %1u %2u %10.6e %10.6e %10.6e "
              "%10.6e %10.6e %10.6e %10.6e\n",
              chl_index, carotenoid, ii, e_xa, fc, ji_work, 
              trapezoid(fi_ad, 2. * M_PI / (TOFS * tau), tau),
              trapezoid(ai_fd, 2. * M_PI / (TOFS * tau), tau),
              chl_car,
              car_chl
              );
        }


        /* only calculate this once for each carotenoid lol */
        if (chl_index == 0) {
          for (unsigned kk = 0; kk < vib_totals[car_index]; kk++) {

            /* NB: indexing!!! */
            car_decays[car_index][ii - vib_totals[car_index]] += 
              car_rates[car_index][sub2ind({kk, ii},
                  {n_s_totals[car_index], n_s_totals[car_index]})];
            if (print_decay_details) {
              fprintf(stdout, "%2u %2u %10.6e\n",
                  ii, kk, car_rates[car_index][sub2ind({kk, ii},
                  {n_s_totals[car_index], n_s_totals[car_index]})]);
            }

          } // kk

        }

        if (output_lineshapes) {
          /* both carotenoids have same lineshapes */
          /* this is a nightmare lol strings in C!!!
           * snprintf will add a null char which we don't want,
           * so define another char[2] and memcpy the formatted 
           * number to that without the null char. then gen filename.
           * set the number back to 01 each time by redefining fn,
           * otherwise we'd have to keep track of the previous number
           * and search for that. i guess alternatively we could
           * hardcode the index of the underscore and overwrite
           * two characters after that. but who cares. it's ugly either way
           * */
          char snum[3], car_num[2], num[2];
          if (carotenoid == n_chl) {
            char fn[200] = "out/NLLZ/7_PROD_2/2/1000/car_01.dat\0";
            int status = snprintf(snum, 3, "%02lu", ii - vib_totals[car_index]);
            memcpy(car_num, snum, 2);
            status = generate_filename(sizeof(fn), fn, "01", snum);
            if (status != 0) {
              fprintf(stdout, "filename generation for outputting"
                  "car lineshape %d didn't work\n", ii);
              break;
            } else {
              FILE *fp = fopen(fn, "w");
              fprintf(fp, "# chl_car = %10.6e\n",
                      chl_car * CM_PER_PS * PI * 2);
              fprintf(fp, "# car_chl = %10.6e\n",
                      car_chl * CM_PER_PS * PI * 2);

              for (unsigned step = 0; step < tau; step++) {
                fprintf(fp, "%10.6e\n", abs[step]);
              }
              int cl = fclose(fp);
              if (cl != 0) {
                fprintf(stdout, "Failed to close car lineshape"
                    "output file no. %d, error no. %d.\n", ii, cl);
                exit(EXIT_FAILURE);
              }
            }
          }
          if (ii == vib_totals[car_index]) {
            char fn[200] = "out/NLLZ/7_PROD_2/2/1000/chl_01.dat\0";
            int status = snprintf(snum, 3, "%02u", chl_index);
            memcpy(num, snum, 2);
            /* fprintf(stdout, "num = %s", num); */
            status = generate_filename(sizeof(fn), fn, "01", snum);
            /* fprintf(stdout, "%s", fn); */
            if (status != 0) {
              fprintf(stdout, "filename generation for outputting"
                  "chl lineshape %d didn't work\n", chl_index);
              break;
            } else {
              FILE *fp = fopen(fn, "w");
              for (unsigned step = 0; step < tau; step++) {
                fprintf(fp, "%10.6e %10.6e\n",
                    normed_ai[chl_index][step],
                    normed_fi[chl_index][step]);
              }
              int cl = fclose(fp);
              if (cl != 0) {
                fprintf(stdout, "Failed to close chl lineshape"
                    "output file no. %d, error no. %d.\n", chl_index, cl);
                exit(EXIT_FAILURE);
              }
            }
          }
        } // output_lineshapes

        chl_car *= CM_PER_PS * 2. * PI;
        car_chl *= CM_PER_PS * 2. * PI;
        if (eigvals[chl_index] > e_xa) {
          car_chl *= exp(-beta * (eigvals[chl_index] - e_xa));
        } else {
          chl_car *= exp(-beta * (e_xa - eigvals[chl_index]));
        }
        k_i_xa.push_back(chl_car);
        k_i_xa.push_back(car_chl);

        if (print_details) {
          fprintf(stdout, "chl = %2u, car = %2u,"
          " state = (%2lu, %2lu, %2lu): "
          /* "e_chl = %10.6e, e_car = %10.6e, " */
          "FC = %10.6e, chl->car = %10.6e,"
          " car->chl = %10.6e\n", chl_index, carotenoid,
          a[0], a[1], a[2],
          /* eigvals[chl_index], e_xa, */
          fc, chl_car, car_chl);
        }

      } // ii
    }
  } // chl_index
  free(abs);
  free(fi_ad);
  free(ai_fd);

  return k_i_xa;
}

double**
car_transfer(VERA x, double *decays)
{
  std::vector<size_t> extents = x.get_pop_extents();
  std::vector<size_t> subscripts(extents.size(), 0.),
                      j_subs(extents.size(), 0.);

  unsigned s1_total = pow(x.n_vib + 1, x.n_normal) + 1;
  unsigned n_total = x.n_total;
  std::vector<double> rates = x.intra_rates();

  double **t = (double **)calloc(n_total, sizeof(double*));
  for (unsigned i = 0; i < n_total; i++) {
    t[i] = (double *)calloc(n_total, sizeof(double));
  }

  /* indices - ground state is index 0 */
  for (unsigned i = 1; i <= s1_total; i++) {
    subscripts = ind2sub(i - 1, extents);  
    /* this will return 0, a_1, a_2, so add one to the electronic index */
    subscripts[0]++;
    size_t ind = sub2ind(subscripts, extents);

    for (unsigned j = 1; j <= s1_total; j++) {
      if (i == j) {
        t[i - 1][j - 1] -= decays[i - 1];
        t[0][j - 1] += decays[i - 1];
      } else {
        j_subs = ind2sub(j - 1, extents);

        j_subs[0]++;
        size_t j_ind = sub2ind(j_subs, extents);
        t[j - 1][i - 1] = rates[sub2ind({ind, j_ind}, {n_total, n_total})];
      }

    }

  }
  return t;
}

double**
hybrid_transfer(unsigned n_chl, unsigned n_car, std::vector<VERA> x,
    double *gamma, double **Jij, std::vector<double> k_i_delta,
    double **redfield_rates, double **car_decays)
{
  /* n_vib_tot is the total number of vibrational levels
   * per electronic state; we only need this here because
   * we're not concerned with hot ground states/S2 */
  std::vector<size_t> n_vib_tot(n_car, 0);
  size_t size = n_chl + 1;
  for (unsigned i = 0; i < n_car; i++) {
    n_vib_tot[i] = pow(x[i].n_vib + 1, x[i].n_normal);
    size += n_vib_tot[i] + 1;
  }

  double **k_tot = (double **)calloc(size, sizeof(double*));
  for (unsigned i = 0; i < size; i++) {
    k_tot[i] = (double *)calloc(size, sizeof(double));
  }

  std::vector<std::vector<size_t>> car_extents;
  for (unsigned i = 0; i < n_car; i++) {
    car_extents.push_back(x[i].get_pop_extents());
  }

  /* NB: this assumes n_vib_tot is the same for each carotenoid!
   * this could be fixed pretty easily - make k_i_xa_hybrid return
   * a vector of vectors, one for each carotenoid, then make the 
   * below a vector of vectors in the same way
   */
  std::vector<size_t> chl_car_extents = {n_chl, n_car, n_vib_tot[0], 2};
  std::vector<std::vector<double>> car_rates;
  for (unsigned i = 0; i < n_car; i++) {
    car_rates.push_back(x[i].intra_rates());
  }

  std::vector<size_t> car_gs_indices(n_car, 0);
  car_gs_indices[0] = n_chl + 1;
  for (unsigned k = 1; k < n_car; k++) {
    car_gs_indices[k] = car_gs_indices[k - 1] + n_vib_tot[k - 1];
    fprintf(stdout, "%1u: %2lu\n", k, car_gs_indices[k]);
  }

  std::vector<bool> i_in_car(n_car, false);
  std::vector<bool> j_in_car(n_car, false);
  for (unsigned i = 0; i < size; i++) {
    for (unsigned j = 0; j < size; j++) {
      unsigned rgs = 0;
      for (unsigned k = 0; k < n_car; k++) {
        /* need to reset each time */
        i_in_car[k] = false;
        j_in_car[k] = false;
      }

      for (unsigned k = 0; k < n_car; k++) {
        /* break statements ensure only the first
         * case will be set to true */
        /* cast both to int bc otherwise you get wraparound */
        if ((int)i - (int)car_gs_indices[k] > 0) {
          i_in_car[k] = true;
          break;
        }
      }

      for (unsigned k = 0; k < n_car; k++) {
        if ((int)j - (int)car_gs_indices[k] > 0) {
          j_in_car[k] = true;
          break;
        }
      }
      /* unsigned gs_620 = n_chl + 1; */
      /* unsigned gs_621 = n_vib_tot + 1 + gs_620; */
      bool i_rgs  = (i == rgs);
      bool j_rgs  = (j == rgs);

      bool i_chls = (i > rgs && i <= n_chl);
      bool j_chls = (j > rgs && j <= n_chl);
      /* bool i_620_gs  = (i == gs_620); */
      /* bool j_620_gs  = (j == gs_620); */
      /* bool i_620  = (i > gs_620 && i < gs_621); */
      /* bool j_620  = (j > gs_620 && j < gs_621); */
      /* bool i_621_gs  = (i == gs_621); */
      /* bool j_621_gs  = (j == gs_621); */
      /* bool i_621  = (i > gs_621); */
      /* bool j_621  = (j > gs_621); */

      /* fprintf(stdout, "i = %2u, j = %2u, i_chls = %1d, j_chls = %1d ", */
      /*     i, j, i_chls, j_chls); */
      /* for (unsigned k = 0; k < n_car; k++) { */
      /*   /1* fprintf(stdout, "%1u %1d %1d", k, i_in_car[k], j_in_car[k]); *1/ */
      /*   std::cout << "car " << k << ": i_in_car =" << i_in_car[k] */ 
      /*     << ", j_in_car = " << j_in_car[k]; */
      /* } */
      /* fprintf(stdout, "\n"); */

      if (i_rgs) {
        if (j_rgs) {
          continue;
        }
        if (j_chls) {
          k_tot[i][j] += (1. / (1000 * gamma[j - 1]));
        }
        /* if (j_620) { */
        /*   continue; */
        /* } */
        /* if (j_621) { */
        /*   continue; */
        /* } */
      }

      if (j_rgs) {
        /* nothing comes out of Redfield ground state */
        continue;
      }

      if (i_chls) {
        if (j_rgs) {
          continue;
        }
        if (j_chls) {
          /* redfield rate: i/j - 1 because of ground state */
          /* [i][j] not [j][i] because this is already a transfer matrix */
          k_tot[i][j] = redfield_rates[i - 1][j - 1];
          if (i == j) {
            /* need to subtract the outward rates to the carotenoids */
            for (unsigned car = 0; car < n_car; car++) {
              for (unsigned k = 0; k < n_vib_tot[car]; k++) {
                k_tot[i][j] -= k_i_delta[sub2ind({i - 1, car, k, 0},
                               chl_car_extents)];
              }
            }
          }
        }
      } // i_chls

      for (unsigned k = 0; k < n_car; k++) {
        if (i_chls) {
          if (j_in_car[k]) {
            k_tot[j][i] = k_i_delta[sub2ind({i - 1, 0,
                j - (car_gs_indices[k] + 1), 0}, chl_car_extents)];
          }
        }
        if (i == car_gs_indices[k]) {
          if (j_in_car[k]) {
            k_tot[i][j] += car_decays[k][j - car_gs_indices[k] + 1];
          }
        }
        if (i_in_car[k]) {
          if (j_chls) {
            k_tot[j][i] += k_i_delta[sub2ind({j - 1, 0,
                i - (car_gs_indices[k] + 1), 1}, chl_car_extents)];
          }
          if (j_in_car[k]) {
            if (i == j) {
              k_tot[i][j] -= car_decays[k][j - (car_gs_indices[k] + 1)];
              for (unsigned m = 0; m < n_chl; m++) {
                k_tot[i][j] -= k_i_delta[sub2ind({m, 0,
                    j - (car_gs_indices[k] + 1), 1}, chl_car_extents)];
              }
            }
            std::vector<size_t> i_subs = ind2sub(i - (car_gs_indices[k] + 1),
                car_extents[k]);
            std::vector<size_t> j_subs = ind2sub(j - (car_gs_indices[k] + 1),
                car_extents[k]);
            /* S1! */
            i_subs[0]++; j_subs[0]++;
            size_t i_index = sub2ind(i_subs, car_extents[k]);
            size_t j_index = sub2ind(j_subs, car_extents[k]);
            k_tot[i][j] = car_rates[k][sub2ind({i_index, j_index},
                {x[k].n_total, x[k].n_total})];
          }
        }
      }

      /* if (i_620_gs) {*/
      /*   if (j_rgs) {*/
      /*   }*/
      /*   if (j_chls) {*/
      /*   }*/
      /*   if (j_620) {*/
      /*     k_tot[i][j] += car_decays[j - (gs_620 + 1)];*/
      /*   }*/
      /*   if (j_621) {*/
      /*   }*/
      /* } // i_620_gs*/

      /* if (i_620) {*/
      /*   if (j_rgs) {*/
      /*   }*/
      /*   if (j_chls) {*/
      /*     k_tot[j][i] += k_i_delta[sub2ind({j - 1, 0,*/
      /*         i - (gs_620 + 1), 1}, chl_car_extents)];*/
      /*   }*/
      /*   if (j_620) {*/
      /*     if (i == j) {*/
      /*       k_tot[i][j] -= car_decays[j - (gs_620 + 1)];*/
      /*       for (unsigned k = 0; k < n_chl; k++) {*/
      /*         k_tot[i][j] -= k_i_delta[sub2ind({k, 0, j - (gs_620 + 1), 1},*/
      /*                       chl_car_extents)];*/
      /*       }*/
      /*     }*/
           /* this needs confirming, but:
            * the total rate matrix, VERA.intra_rates(), contains
            * both the IVR and IC rates. but on a given electronic
            * state, only IVR is nonzero, and across electronic
            * states, only IC is nonzero. hence, if we ensure we're
            * on the same electronic state, intra_rates()[i][j]
            * should just be the IVR rate.*/
      /*     std::vector<size_t> i_subs = ind2sub(i - (gs_620 + 1), car_extents);*/
      /*     std::vector<size_t> j_subs = ind2sub(j - (gs_620 + 1), car_extents);*/
           /* S1! */
      /*     i_subs[0]++; j_subs[0]++;*/
      /*     size_t i_index = sub2ind(i_subs, car_extents);*/
      /*     size_t j_index = sub2ind(j_subs, car_extents);*/
      /*     k_tot[i][j] = car_rates[sub2ind({i_index, j_index},*/
      /*         {x->n_total, x->n_total})];*/
      /*   }*/
      /*   if (j_621) {*/
      /*   }*/
      /* }*/

      /* this can be folded into the 620 bit i guess; these blocks
       * will be identical, it's only the chl-car blocks which differ */

      /* if (i_621_gs) {*/
      /*   if (j_rgs) {*/
      /*   }*/
      /*   if (j_chls) {*/
      /*   }*/
      /*   if (j_620) {*/
      /*   }*/
      /*   if (j_621) {*/
      /*     k_tot[i][j] += car_decays[j - (gs_621 + 1)];*/
      /*   }*/
      /* } // i_620_gs*/
      /* if (i_621) {*/
      /*   if (j_rgs) {*/
      /*   }*/
      /*   if (j_chls) {*/
      /*     k_tot[j][i] += k_i_delta[sub2ind({j - 1, 1,*/
      /*         i - (gs_621 + 1), 1}, chl_car_extents)];*/
      /*   }*/
      /*   if (j_620) {*/
      /*   }*/
      /*   if (j_621) {*/
      /*     if (i == j) {*/
      /*       k_tot[i][j] -= car_decays[j - (gs_621 + 1)];*/
      /*       for (unsigned k = 0; k < n_chl; k++) {*/
      /*         k_tot[i][j] -= k_i_delta[sub2ind({k, 1,*/
      /*             j - (gs_621 + 1), 1}, chl_car_extents)];*/
      /*       }*/
      /*     }*/
      /*     std::vector<size_t> i_subs = ind2sub(i - (gs_621 + 1),*/
      /*                                  car_extents);*/
      /*     std::vector<size_t> j_subs = ind2sub(j - (gs_621 + 1),*/
      /*                                  car_extents);*/
      /*     i_subs[0]++; j_subs[0]++;*/
      /*     size_t i_index = sub2ind(i_subs, car_extents);*/
      /*     size_t j_index = sub2ind(j_subs, car_extents);*/
      /*     k_tot[i][j] = car_rates[sub2ind({i_index, j_index},*/
      /*         {x->n_total, x->n_total})];*/
      /*   }*/
      /* }*/

    } // j
  }
  return k_tot;

}

double*
hybrid_boltz(unsigned n_chl, unsigned n_car, double beta,
             double *eigvals, std::vector<VERA> x)
{
  /* get the boltzmann factor for every state and normalise, 
   * then return just the chlorophyll boltzmann factors. we do this
   * because the carotenoid doesn't fluoresce but does have some
   * population on it. note that, as everywhere else in the code,
   * this assumes all the carotenoids are identical */
  std::vector<size_t> n_vib_tot(n_car, 0);
  std::vector<size_t> car_starts(n_car, 0);
  size_t size = n_chl;
  car_starts[0] = n_chl;
  for (unsigned i = 0; i < n_car; i++) {
    n_vib_tot[i] = pow(x[i].n_vib + 1, x[i].n_normal);
    size += n_vib_tot[i];
    if (i > 0) {
      car_starts[i] = car_starts[i - 1] + n_vib_tot[i - 1];
    }
  }
  double *work = (double *)calloc(size, sizeof(double));
  double *boltz = (double *)calloc(n_chl, sizeof(double));
  std::vector<std::vector<size_t>> extents;
  for (unsigned i = 0; i < n_car; i++) {
    extents.push_back(x[i].get_pop_extents());
  }
  double work_sum = 0.;

  std::vector<bool> i_in_car(n_car, false);
  for (unsigned i = 0; i < size; i++) {
    if (i < n_chl) {
      work[i] = exp(-beta * eigvals[i]); // dimensions!
      work_sum += work[i];
    } else {
      for (unsigned k = 0; k < n_car; k++) {
        i_in_car[k] = false;
      }
      for (unsigned k = 0; k < n_car; k++) {
        if (i - car_starts[k] > 0) {
          i_in_car[k] = true;
          break;
        }
      }
      for (unsigned k = 0; k < n_car; k++) {
        if (i_in_car[k]) {
        std::vector<size_t> subs = ind2sub(i - car_starts[k], extents[k]);
        subs[0]++;
        double e_xa = x[k].get_w_elec(subs[0]);
        for (unsigned alpha = 0; alpha < x[k].n_normal; alpha++) {
          e_xa += x[k].get_w_normal(alpha) * subs[alpha + 1];
        }

        work[i] = exp(-beta * e_xa);
        work_sum += work[i];
        }
      }

    }
  }

  for (unsigned i = 0; i < size; i++) {
    work[i] /= work_sum;
  }
  for (unsigned i = 0; i < n_chl; i++) {
    boltz[i] = work[i];
  }
  free(work);
  return boltz;

}
