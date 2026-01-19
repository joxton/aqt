
#include <R.h>
#include <R_ext/Utils.h>
#include <Rmath.h>
#include <assert.h>
#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

void revsorts(double *a, int n);
void quantile_sorted(double *x, int *n, double *tau, int *ntau, double *q);

int compare_doubles(const void *a, const void *b) {
  // Cast the void pointers to double pointers and dereference them
  double val_a = *(const double *)a;
  double val_b = *(const double *)b;

  // Return an integer based on the comparison
  if (val_a < val_b) {
    return -1; // Indicates val_a comes before val_b
  } else if (val_a > val_b) {
    return 1; // Indicates val_a comes after val_b
  } else {
    return 0; // Indicates they are equal
  }
}

void quantile(double *x, int *n, double *tau, int *ntau, double *res) {

  qsort(x, *n, sizeof(double), compare_doubles);
  quantile_sorted(x, n, tau, ntau, res);

  // revsorts(x, *n);
}

int sample_indx(int n) {
  // sample an index from 0:(n-1)
  int i;
  i = (int)((unif_rand()) * ((double)n));
  if (i == n)
    i = n - 1;
  return (i);
}

void bootstrap_quantile(double *x, int *n,
                        int *B) { // , double *tau, int *ntau, double *res){

  int *indx_counter = (int *)malloc(sizeof(int) * (*n));
  double *xb = (double *)malloc(sizeof(double) * (*n));

  qsort(x, *n, sizeof(double), compare_doubles);

  for (int b = 0; b < *B; b++) {
    for (int i = 0; i < *n; i++)
      indx_counter[i] = 0;
    for (int i = 0; i < *n; i++) {
      indx_counter[sample_indx(*n)]++;
    }
    int counter = 0;
    for (int i = 0; i < *n; i++) {
      if (indx_counter[i] > 0) {
        for (int j = 0; j < indx_counter[i]; j++) {
          xb[counter] = x[i];
          counter++;
        }
      }
    }
  }

  free(indx_counter);
  free(xb);
}

void quantile_sorted(double *x, int *n, double *tau, int *ntau, double *q) {
  //  x_1, ... x_n
  //  (k - 1)/(n - 1) <= tau < k/(n - 1) -> (k-1) <= tau*(n-1); -> k-1 =
  //  floor(tau(n-1)) - 1
  double p_k0, p_k1;
  int indx;

  for (int i = 0; i < *ntau; i++) {
    indx = ((int)floor(((double)(*n - 1)) * tau[i]));
    p_k0 = ((double)(indx)) / ((double)(*n - 1));
    p_k1 = ((double)(indx + 1)) / ((double)(*n - 1));
    // q[i] = x[indx] + (x[indx + 1] - x[indx])*(tau[i] - p_k0)/(p_k1 - p_k0);
    q[i] = x[indx] + (x[indx + ((indx < (*n - 1)) ? 1 : 0)] - x[indx]) *
                         (tau[i] - p_k0) / (p_k1 - p_k0);
  }
}

double quantile_srtd(double *x, int n, double tau) {
  // Return quantile for single tau, assuming sorted x (ascending)
  //  x_1, ... x_n
  //  (k - 1)/(n - 1) <= tau < k/(n - 1) -> (k-1) <= tau*(n-1); -> k-1 =
  //  floor(tau(n-1)) - 1
  double p_k0, p_k1;
  int indx;
  double q;

  indx = ((int)floor(((double)(n - 1)) * tau));
  p_k0 = ((double)(indx)) / ((double)(n - 1));
  p_k1 = ((double)(indx + 1)) / ((double)(n - 1));
  // q[i] = x[indx] + (x[indx + 1] - x[indx])*(tau[i] - p_k0)/(p_k1 - p_k0);
  q = x[indx] + (x[indx + ((indx < (n - 1)) ? 1 : 0)] - x[indx]) *
                    (tau - p_k0) / (p_k1 - p_k0);
  return (q);
}

void bootstrap_sample_quantile(double *x, int *n, int *B, double *tau,
                               int *ntau, double *q, double *qvar) {
  // Bootstrap variance for sample quantiles

  int *indx_counter = (int *)malloc(sizeof(int) * (*n));
  double *xb = (double *)malloc(sizeof(double) * (*n));
  double *q_sum = (double *)malloc(sizeof(double) * (*ntau));
  double *q_sumsq = (double *)malloc(sizeof(double) * (*ntau));

  for (int i = 0; i < *ntau; i++) {
    q_sum[i] = 0;
    q_sumsq[i] = 0;
  }

  qsort(x, *n, sizeof(double), compare_doubles);

  for (int b = 0; b < *B; b++) {
    for (int i = 0; i < *n; i++)
      indx_counter[i] = 0;
    for (int i = 0; i < *n; i++) {
      indx_counter[sample_indx(*n)]++;
    }
    int counter = 0;
    for (int i = 0; i < *n; i++) {
      if (indx_counter[i] > 0) {
        for (int j = 0; j < indx_counter[i]; j++) {
          xb[counter] = x[i];
          counter++;
        }
      }
    }
    quantile_sorted(xb, n, tau, ntau, q);
    for (int i = 0; i < *ntau; i++) {
      q_sum[i] += q[i];
      q_sumsq[i] += (q[i] * q[i]);
    }
  }
  for (int i = 0; i < *ntau; i++) {
    qvar[i] = (q_sumsq[i] / ((double)*B)) -
              (q_sum[i] / ((double)*B)) * (q_sum[i] / ((double)*B));
  }
  quantile_sorted(x, n, tau, ntau, q);

  free(indx_counter);
  free(xb);
  free(q_sum);
  free(q_sumsq);
}

void bootstrap_quantile_difference(double *x, int *n, double *y, int *m, int *B,
                                   double *tau, int *ntau, double *q,
                                   double *qvar) {
  // Bootstrap variance for quantiles differences
  // The R function for quantile of binominal: (could use in approximation)
  // double qbinom(double p, double n, double pr, int lower_tail, int log_p)

  // An approximation that uses se(q) = sqrt(p(1-p)/n) / f(q). in estimation
  // could be used (f(q)  is known there)
  double d;
  int counter;
  int mn;
  mn = (*n > *m) ? *n : *m;

  int *indx_counter = (int *)malloc(sizeof(int) * (mn));
  double *xb = (double *)malloc(sizeof(double) * (mn));
  double *q_sum = (double *)malloc(sizeof(double) * (*ntau));
  double *q_sumsq = (double *)malloc(sizeof(double) * (*ntau));
  double *qy = (double *)malloc(sizeof(double) * (*ntau));

  for (int i = 0; i < *ntau; i++) {
    q_sum[i] = 0;
    q_sumsq[i] = 0;
  }

  qsort(x, *n, sizeof(double), compare_doubles);
  qsort(y, *m, sizeof(double), compare_doubles);

  for (int b = 0; b < *B; b++) {
    // Bootstrap x
    for (int i = 0; i < *n; i++)
      indx_counter[i] = 0;
    for (int i = 0; i < *n; i++) {
      indx_counter[sample_indx(*n)]++;
    }
    counter = 0;
    for (int i = 0; i < *n; i++) {
      if (indx_counter[i] > 0) {
        for (int j = 0; j < indx_counter[i]; j++) {
          xb[counter] = x[i];
          counter++;
        }
      }
    }
    quantile_sorted(xb, n, tau, ntau, q);
    // Bootstrap y
    for (int i = 0; i < *m; i++)
      indx_counter[i] = 0;
    for (int i = 0; i < *m; i++) {
      indx_counter[sample_indx(*m)]++;
    }
    counter = 0;
    for (int i = 0; i < *m; i++) {
      if (indx_counter[i] > 0) {
        for (int j = 0; j < indx_counter[i]; j++) {
          xb[counter] = y[i];
          counter++;
        }
      }
    }
    quantile_sorted(xb, m, tau, ntau, qy);

    for (int i = 0; i < *ntau; i++) {
      d = q[i] - qy[i];
      q_sum[i] += d;
      q_sumsq[i] += (d * d);
    }
  }

  quantile_sorted(x, n, tau, ntau, q);
  quantile_sorted(y, m, tau, ntau, qy);

  for (int i = 0; i < *ntau; i++) {
    qvar[i] = (q_sumsq[i] / ((double)*B)) -
              (q_sum[i] / ((double)*B)) * (q_sum[i] / ((double)*B));
    q[i] = (q[i] - qy[i]);
  }

  free(indx_counter);
  free(xb);
  free(q_sum);
  free(q_sumsq);
  free(qy);
}

void permute(double *z, double *x, int n, double *y, int m) {
  // Randomly permute sorted values of z's into x (size n) and y (size m)
  // preserving order
  int total = n + m;
  int my = 0;
  int nx = 0;
  double py, u;

  for (int i = 0; i < (n + m); i++) {

    if (my < m && nx < n) {
      py = ((double)(m - my)) / ((double)(n - nx) + (m - my));
      u = unif_rand();
      // Rprintf("py =%lf u =%lf \n", py, u);
      if (u <= py) {
        y[my] = z[i];
        my++;
      } else {
        x[nx] = z[i];
        nx++;
      }
    } else {
      if (my == m) {
        x[nx] = z[i];
        nx++;
      } else {
        y[my] = z[i];
        my++;
      }
    }
  }
}

void permutation_quantile_difference(double *x, int *n, double *y, int *m,
                                     int *B, double *tau, int *ntau, int *R,
                                     int *all, double *permdist) {
  // Bootstrap variance for quantiles differences

  double a;
  double d;
  int counter;
  int mn;
  mn = (*n > *m) ? *n : *m;

  int *indx_counter = (int *)malloc(sizeof(int) * (mn));
  double *xb = (double *)malloc(sizeof(double) * (mn));
  double *q_sum = (double *)malloc(sizeof(double) * (*ntau));
  double *q_sumsq = (double *)malloc(sizeof(double) * (*ntau));
  double *qy = (double *)malloc(sizeof(double) * (*ntau));
  double *z = (double *)malloc(sizeof(double) * (*n + *m));
  double *xstar = (double *)malloc(sizeof(double) * (*n));
  double *ystar = (double *)malloc(sizeof(double) * (*m));
  double *q = (double *)malloc(sizeof(double) * (*ntau));
  double *qvar = (double *)malloc(sizeof(double) * (*ntau));

  for (int i = 0; i < *R; i++)
    permdist[i] = 0;

  for (int i = 0; i < (*n + *m); i++)
    z[i] = (i < *n) ? x[i] : y[i - *n];

  qsort(z, (*n + *m), sizeof(double), compare_doubles);

  // Permutation loop
  for (int ip = 0; ip < *R; ip++) {
    permute(z, xstar, *n, ystar, *m);
    //		for(int i = 0; i < 10; i++)
    //			Rprintf("xstar[%d] = %lf ystar[%d] = %lf \n", i,
    // xstar[i], i, ystar[i],i);

    for (int i = 0; i < *ntau; i++) {
      q_sum[i] = 0;
      q_sumsq[i] = 0;
    }

    //  qsort(x, *n, sizeof(double), compare_doubles);
    //  qsort(y, *m, sizeof(double), compare_doubles);

    for (int b = 0; b < *B; b++) {
      // Bootstrap x
      for (int i = 0; i < *n; i++)
        indx_counter[i] = 0;
      for (int i = 0; i < *n; i++) {
        indx_counter[sample_indx(*n)]++;
      }
      counter = 0;
      for (int i = 0; i < *n; i++) {
        if (indx_counter[i] > 0) {
          for (int j = 0; j < indx_counter[i]; j++) {
            xb[counter] = xstar[i];
            counter++;
          }
        }
      }
      quantile_sorted(xb, n, tau, ntau, q);
      // Bootstrap y
      for (int i = 0; i < *m; i++)
        indx_counter[i] = 0;
      for (int i = 0; i < *m; i++) {
        indx_counter[sample_indx(*m)]++;
      }
      counter = 0;
      for (int i = 0; i < *m; i++) {
        if (indx_counter[i] > 0) {
          for (int j = 0; j < indx_counter[i]; j++) {
            xb[counter] = ystar[i];
            counter++;
          }
        }
      }
      quantile_sorted(xb, m, tau, ntau, qy);

      for (int i = 0; i < *ntau; i++) {
        d = q[i] - qy[i];
        q_sum[i] += d;
        q_sumsq[i] += (d * d);
      }
    }

    quantile_sorted(xstar, n, tau, ntau, q);
    quantile_sorted(ystar, m, tau, ntau, qy);

    for (int i = 0; i < *ntau; i++) {
      qvar[i] = (q_sumsq[i] / ((double)*B)) -
                (q_sum[i] / ((double)*B)) * (q_sum[i] / ((double)*B));
      q[i] = (q[i] - qy[i]);

      if (*all) {
        permdist[ip * (*ntau) + i] = q[i] / sqrt(qvar[i]);
      } else {
        a = q[i] * q[i] / qvar[i];
        if (a > permdist[ip])
          permdist[ip] = a;
      }
    }
  }
  free(indx_counter);
  free(xb);
  free(q_sum);
  free(q_sumsq);
  free(qy);
  free(z);
  free(xstar);
  free(ystar);
  free(q);
  free(qvar);
}
