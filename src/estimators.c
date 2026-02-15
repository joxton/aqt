
#include <R.h>
#include "aqt.h"
#include "pcdensity.h"
#include <R_ext/Utils.h>
#include <Rmath.h>
#include <assert.h>
#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>


// .......................................................................
// ... from gemini.google.com : generating thread safe uniforms 
// Minimal PCG32 implementation for thread-local RNG
typedef struct { uint64_t state;  uint64_t inc; } pcg32_state;

uint32_t pcg32_random(pcg32_state* rng) {
    uint64_t oldstate = rng->state;
    rng->state = oldstate * 6364136223846793005ULL + (rng->inc | 1);
    uint32_t xorshifted = ((oldstate >> 18u) ^ oldstate) >> 27u;
    uint32_t rot = oldstate >> 59u;
    return (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
}

// Helper to get a double between [0, 1)
double pcg32_double(pcg32_state* rng) {
    return (double)pcg32_random(rng) / (double)UINT32_MAX;
}



void parallel_random_fill( int *nr, int *seedr, double *out) {
	int seed; 
	int n; 
	n = *nr; 
	seed = *seedr; 
    #pragma omp parallel
    {
        // 1. Initialize a private state per thread
        int tid = omp_get_thread_num();
        pcg32_state rng;
        rng.state = (uint64_t)seed + tid; // Unique start per thread
        rng.inc = (uint64_t)tid << 1 | 1;

        #pragma omp for
        for (int i = 0; i < n; i++) {
            // 2. Generate numbers without touching R's global state
            out[i] = pcg32_double(&rng);
        }
    }
}

// .......................................................................


int sample_indx_parallel(int n, pcg32_state* rng) {
  // sample an index from 0:(n-1)
  int i;
  i = (int)floor((pcg32_double(rng)) * ((double)n));
  if (i == n)
    i = n - 1;
  return (i);
}




struct two_sample {

  // x and y will be sorted in ascending order
  double *x;
  int nx;
  double *y;
  int ny;

  // Piecewise constant densities for x and y (if applicable)
  struct pcd *px;
  struct pcd *py;

  // Quantile info (if applicable)
  double *q;    // quantile difference
  double *qvar; // variance estimate of quantile difference
  double *tau;  // probabilities for quantiles
  int ntau;     // length of tau (and q, qvar)

  // Quantile-wise permutation p-values
  double **pvalue;

  int method;
  int combiner;

  // Test specific data
  int B;
  int R;
  int delete_d;
};

void revsortx(double *a, int *ib, int n);
void sort_descending_with_order(double *a, int *ib, int n);
void sort_descending(double *a, int n);
void qdiff_quick_and_dirty(struct two_sample *t);
void qdiff_permutation(struct two_sample *t);
void qdiff_bootstrap(struct two_sample *t);
void qdiff_jackknife_delete_d(struct two_sample *t);
void qdiff_permute_and_combine(struct two_sample *t);
void qdiff_bootstrap_parallell(struct two_sample *t);


void adaptive_fisher(struct two_sample *t);
void minp(struct two_sample *t);
void order_statistic(struct two_sample *t);

void two_sample_aqt(double *x, int *n, double *y, int *m, double *tau,
                    int *ntau, int *R, int *method, double *q, double *qvar) {

  GetRNGstate();
  struct two_sample ts;

  // Set up ts
  qsort(x, *n, sizeof(double), compare_doubles);
  qsort(y, *m, sizeof(double), compare_doubles);

  ts.x = x;
  ts.y = y;
  ts.nx = *n;
  ts.ny = *m;
  ts.q = q;
  ts.qvar = qvar;
  ts.tau = tau;
  ts.ntau = *ntau;
  ts.R = *R;
  ts.method = method[0];
  ts.combiner = method[1];

  if (ts.method == 4) {
    ts.delete_d = R[1];
  }
  // ts.B = *B;
  // Rprintf("method = %d",*method);

  switch (ts.method) {

  case 1:
    qdiff_quick_and_dirty(&ts);
    break;

  case 2:
    // Permutation testing
    qdiff_permutation(&ts);
    break;

  case 3:
    qdiff_bootstrap(&ts);
    break;

  case 4:
    qdiff_jackknife_delete_d(&ts);
    break;

  case 5:
    qdiff_permute_and_combine(&ts);
    break;

  case 6:
    qdiff_bootstrap_parallell(&ts);
    break;

  default:
    qdiff_quick_and_dirty(&ts);
  }
  /**/

  PutRNGstate();

  //
}

void qdiff_permutation(struct two_sample *t) {
  double qx, qy;

  double *z = (double *)malloc(sizeof(double) * (t->nx + t->ny));
  double *xstar = (double *)malloc(sizeof(double) * t->nx);
  double *ystar = (double *)malloc(sizeof(double) * t->ny);

  t->pvalue = (double **)malloc(sizeof(double *) * (t->R));
  t->pvalue[0] = (double *)malloc(sizeof(double) * (t->R * t->ntau));

  for (int i = 0; i < t->R; i++) {
    t->pvalue[i] = &(t->pvalue[0][i * t->ntau]);
  }

  for (int i = 0; i < (t->nx + t->ny); i++)
    z[i] = (i < t->nx) ? t->x[i] : t->y[i - t->nx];

  qsort(z, (t->nx + t->ny), sizeof(double), compare_doubles);

  // Permutation loop
  for (int ip = 0; ip < t->R; ip++) {
    // ip = 0 is for the observed data ip > 0, the permutations
    if (ip > 0) {
      permute(z, xstar, t->nx, ystar, t->ny);
      t->x = xstar;
      t->y = ystar;
    }

    for (int i = 0; i < t->ntau; i++) {

      qx = quantile_srtd(t->x, t->nx, t->tau[i]);
      qy = quantile_srtd(t->y, t->ny, t->tau[i]);
      t->q[ip + i * t->R] = fabs(qx - qy);
      t->qvar[ip + i * t->R] = t->q[ip + i * t->R];
    }
  }

  int *ib = (int *)malloc(sizeof(int) * t->R);

  // Sort t->[k*R + 1:R]  to sort the k+1-th quantile difference
  for (int i = 0; i < t->ntau; i++) {

    for (int j = 0; j < t->R; j++) {
      ib[j] = j;
    }

    // Sort in descending order,
    revsortx(&(t->q[i * t->R]), ib, t->R);
    for (int j = 0; j < t->R; j++) {
      t->pvalue[ib[j]][i] = ((double)(j)) / ((double)t->R);
    }
  }

  for (int i = 0; i < t->R; i++)
    for (int j = 0; j < t->ntau; j++)
      t->q[i * t->ntau + j] = t->pvalue[i][j];

  free(t->pvalue[0]);
  free(t->pvalue);

  free(ib);
  free(xstar);
  free(ystar);
  free(z);
}

void qdiff_permute_and_combine(struct two_sample *t) {
  double qx, qy;

  double *z = (double *)malloc(sizeof(double) * (t->nx + t->ny));
  double *xstar = (double *)malloc(sizeof(double) * t->nx);
  double *ystar = (double *)malloc(sizeof(double) * t->ny);

  t->pvalue = (double **)malloc(sizeof(double *) * (t->R));
  t->pvalue[0] = (double *)malloc(sizeof(double) * (t->R * t->ntau));

  for (int i = 0; i < t->R; i++) {
    t->pvalue[i] = &(t->pvalue[0][i * t->ntau]);
  }

  for (int i = 0; i < (t->nx + t->ny); i++)
    z[i] = (i < t->nx) ? t->x[i] : t->y[i - t->nx];

  qsort(z, (t->nx + t->ny), sizeof(double), compare_doubles);

  // Permutation loop
  for (int ip = 0; ip < t->R; ip++) {
    // ip = 0 is for the observed data ip > 0, the permutations
    if (ip > 0) {
      permute(z, xstar, t->nx, ystar, t->ny);
      t->x = xstar;
      t->y = ystar;
    }

    for (int i = 0; i < t->ntau; i++) {

      qx = quantile_srtd(t->x, t->nx, t->tau[i]);
      qy = quantile_srtd(t->y, t->ny, t->tau[i]);
      t->q[ip + i * t->R] = fabs(qx - qy);
      t->qvar[ip + i * t->R] = t->q[ip + i * t->R];
    }
  }

  int *ib = (int *)malloc(sizeof(int) * t->R);

  // Sort t->[k*R + 1:R]  to sort the k+1-th quantile difference
  for (int i = 0; i < t->ntau; i++) {

    for (int j = 0; j < t->R; j++) {
      ib[j] = j;
    }

    // Sort in descending order,
    revsortx(&(t->q[i * t->R]), ib, t->R);
    for (int j = 0; j < t->R; j++) {
      t->pvalue[ib[j]][i] = ((double)(j)) / ((double)t->R);
    }
  }

  for (int i = 0; i < t->R; i++)
    for (int j = 0; j < t->ntau; j++)
      t->q[i * t->ntau + j] = t->pvalue[i][j];

  // Combine across quantile comparisons
  if (t->combiner == 1)
    minp(t);

  if (t->combiner == 2)
    adaptive_fisher(t);

  if (t->combiner == 3)
    order_statistic(t);

  free(t->pvalue[0]);
  free(t->pvalue);

  free(ib);
  free(xstar);
  free(ystar);
  free(z);
}

void minp(struct two_sample *t) {

  for (int r = 0; r < t->R; r++) {
    t->qvar[r] = 1;
    for (int i = 0; i < t->ntau; i++) {
      if (t->pvalue[r][i] < t->qvar[r])
        t->qvar[r] = t->pvalue[r][i];
    }
  }

  double pval = 0;
  for (int r = 0; r < t->R; r++) {
    if (t->qvar[r] <= t->qvar[0]) {
      pval += 1.0;
    }
  }

  pval = pval / ((double)t->R);
  t->qvar[0] = pval;
}

void adaptive_fisher(struct two_sample *t) {

  for (int r = 0; r < t->R; r++) {
    for (int i = 0; i < t->ntau; i++) {
      t->pvalue[r][i] = -log(t->pvalue[r][i]);
    }
  }

  for (int r = 0; r < t->R; r++) {
    // Descending order sort
    sort_descending(t->pvalue[r], t->ntau);

    // Form successive combinations
    for (int i = 1; i < t->ntau; i++)
      t->pvalue[r][i] += t->pvalue[r][i - 1];
  }

  //
  int *ib = (int *)malloc(sizeof(int) * t->R);
  double *a = (double *)malloc(sizeof(double) * t->R);

  for (int i = 0; i < t->ntau; i++) {

    for (int j = 0; j < t->R; j++) {
      ib[j] = j;
      a[j] = t->pvalue[j][i];
    }

    sort_descending_with_order(a, ib, t->R);
    for (int j = 0; j < t->R; j++) {
      t->pvalue[ib[j]][i] = ((double)j) / ((double)t->R);
    }
  }

  minp(t);
  free(a);
  free(ib);
}

void order_statistic(struct two_sample *t) {
  // Combine quantile difference p-values using most extreme p-value order
  // statistic

  // Order p-values
  for (int r = 0; r < t->R; r++) {
    // Descending order sort
    sort_descending(t->pvalue[r], t->ntau);
  }

  // Find p-value of order statistics (ordered p-values)
  int *ib = (int *)malloc(sizeof(int) * t->R);
  double *a = (double *)malloc(sizeof(double) * t->R);

  for (int i = 0; i < t->ntau; i++) {
    for (int j = 0; j < t->R; j++) {
      ib[j] = j;
      a[j] = t->pvalue[j][i];
    }

    // P-values of i-th ordered p-value
    sort_descending_with_order(a, ib, t->R);
    for (int j = 0; j < t->R; j++) {
      t->pvalue[ib[j]][i] = ((double)(t->R - j)) / ((double)t->R);
    }
  }

  minp(t);
  free(a);
  free(ib);
}

void qdiff_bootstrap(struct two_sample *t) {
  // Bootstrap variance for quantiles differences
  // The R function for quantile of binominal: (could use in approximation)
  // double qbinom(double p, double n, double pr, int lower_tail, int log_p)

  // An approximation that uses se(q) = sqrt(p(1-p)/n) / f(q). in estimation
  // could be used (f(q)  is known there)
  /* double *x, int *n, double *y, int *m, int *B,
                                   double *tau, int *ntau, double *q,
                                   double *qvar) {
  */
  double *x, *y, *tau, *q, *qvar;
  int *m, *n, *B, *ntau;

  x = t->x;
  y = t->y;
  tau = t->tau;
  q = t->q;
  qvar = t->qvar;
  m = &(t->ny);
  n = &(t->nx);
  B = &(t->R);
  ntau = &(t->ntau);

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

  /*
    qsort(x, *n, sizeof(double), compare_doubles);
    qsort(y, *m, sizeof(double), compare_doubles);
  */

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




void qdiff_bootstrap_parallell(struct two_sample *t) {
  // Bootstrap variance for quantiles differences
  // The R function for quantile of binominal: (could use in approximation)
  // double qbinom(double p, double n, double pr, int lower_tail, int log_p)

  // An approximation that uses se(q) = sqrt(p(1-p)/n) / f(q). in estimation
  // could be used (f(q)  is known there)
  /* double *x, int *n, double *y, int *m, int *B,
                                   double *tau, int *ntau, double *q,
                                   double *qvar) {
  */
  double *x, *y, *tau, *q, *qvar;
  int *m, *n, *B, *ntau;

  int nthreads = omp_get_max_threads();
  omp_set_num_threads(nthreads);


  x = t->x;
  y = t->y;
  tau = t->tau;
  q = t->q;
  qvar = t->qvar;
  m = &(t->ny);
  n = &(t->nx);
  B = &(t->R);
  ntau = &(t->ntau);

 // double d;
  //int counter;
  int mn;
  mn = (*n > *m) ? *n : *m;

 // int *indx_counter = (int *)malloc(sizeof(int) * (mn));
 // double *xb = (double *)malloc(sizeof(double) * (mn));
  double *q_sum = (double *)malloc(sizeof(double) * (*ntau));
  double *q_sumsq = (double *)malloc(sizeof(double) * (*ntau));
//  double *qy = (double *)malloc(sizeof(double) * (*ntau));

  for (int i = 0; i < *ntau; i++) {
    q_sum[i] = 0;
    q_sumsq[i] = 0;
  }

  /*
    qsort(x, *n, sizeof(double), compare_doubles);
    qsort(y, *m, sizeof(double), compare_doubles);
  */
int N;
N = *ntau;

int Bt = (int) *B/nthreads; 

int seed = sample_indx(999999)+1000; 

#pragma omp parallel for reduction(+:q_sum[0:N], q_sumsq[0:N])
for(int nt = 0; nt < nthreads; nt++){ 

  int *indx_counter = (int *)malloc(sizeof(int) * (mn));
  double *xb = (double *)malloc(sizeof(double) * (mn));
  double *qx = (double *)malloc(sizeof(double) * (*ntau));
  double *qy = (double *)malloc(sizeof(double) * (*ntau));

        int tid = omp_get_thread_num();
        pcg32_state rng;
        rng.state = (uint64_t)seed + tid; // Unique start per thread
        rng.inc = (uint64_t)tid << 1 | 1;


  for (int b = 0; b < Bt; b++) {
  
// sample_indx_parallel(int n, pcg32_state* rng)
  
    // Bootstrap x
    for (int i = 0; i < *n; i++)
      indx_counter[i] = 0;
    for (int i = 0; i < *n; i++) {
    	
      indx_counter[sample_indx_parallel(*n,&rng)]++;
    }
  int  counter = 0;
    for (int i = 0; i < *n; i++) {
      if (indx_counter[i] > 0) {
        for (int j = 0; j < indx_counter[i]; j++) {
          xb[counter] = x[i];
          counter++;
        }
      }
    }
    quantile_sorted(xb, n, tau, ntau, qx);
    // Bootstrap y
    for (int i = 0; i < *m; i++)
      indx_counter[i] = 0;
    for (int i = 0; i < *m; i++) {
      indx_counter[sample_indx_parallel(*m,&rng)]++;
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


	double d; 
    for (int i = 0; i < N; i++) {
      d = qx[i] - qy[i];
      q_sum[i] += d;
      q_sumsq[i] += (d * d);
    }

  }
  free(indx_counter);
  free(xb);
  free(qy);  
  free(qx);
  }

 // double *qx = (double *)malloc(sizeof(double) * (*ntau));
  double *qy = (double *)malloc(sizeof(double) * (*ntau));


  quantile_sorted(x, n, tau, ntau, q);
  quantile_sorted(y, m, tau, ntau, qy);


  for (int i = 0; i < *ntau; i++) {
    qvar[i] = (q_sumsq[i] / ((double)*B)) -
              (q_sum[i] / ((double)*B)) * (q_sum[i] / ((double)*B));
    q[i] = (q[i] - qy[i]);
  }

 // free(indx_counter);
 //  free(xb);
  free(q_sum);
  free(q_sumsq);
  free(qy);  
  //free(qx);
}











void qdiff_jackknife_delete_d(struct two_sample *t) {
  // Bootstrap variance for quantiles differences
  // The R function for quantile of binominal: (could use in approximation)
  // double qbinom(double p, double n, double pr, int lower_tail, int log_p)

  // An approximation that uses se(q) = sqrt(p(1-p)/n) / f(q). in estimation
  // could be used (f(q)  is known there)
  /* double *x, int *n, double *y, int *m, int *B,
                                   double *tau, int *ntau, double *q,
                                   double *qvar) {
  */
  double *x, *y, *tau, *q, *qvar;
  int *m, *n, *B, *ntau;
  int delete_d;

  delete_d = t->delete_d;
  x = t->x;
  y = t->y;
  tau = t->tau;
  q = t->q;
  qvar = t->qvar;
  m = &(t->ny);
  n = &(t->nx);
  B = &(t->R);
  ntau = &(t->ntau);

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

  /*
    qsort(x, *n, sizeof(double), compare_doubles);
    qsort(y, *m, sizeof(double), compare_doubles);
  */
  int delete_count;
  double delete_prob;

  for (int b = 0; b < *B; b++) {
    // Bootstrap x
    delete_count = 0;
    for (int i = 0; i < *n; i++)
      indx_counter[i] = 0;
    for (int i = 0; i < *n; i++) {
      delete_prob = ((double)(delete_d - delete_count)) / ((double)(*n - i));
      if (unif_rand() <= delete_prob) {
        indx_counter[i] = 1;
        delete_count++;
      }
      if (delete_count == delete_d)
        break;
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
    quantile_sorted(xb, &counter, tau, ntau, q);
    // Bootstrap y
    delete_count = 0;
    for (int i = 0; i < *m; i++)
      indx_counter[i] = 0;
    for (int i = 0; i < *m; i++) {
      delete_prob = ((double)(delete_d - delete_count)) / ((double)(*m - i));
      if (unif_rand() <= delete_prob) {
        indx_counter[i] = 1;
        delete_count++;
      }
      if (delete_count == delete_d)
        break;
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
    quantile_sorted(xb, &counter, tau, ntau, qy);

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

    qvar[i] =
        qvar[i] * ((double)(*n + *m - 2 * delete_d)) / ((double)(2 * delete_d));
    q[i] = (q[i] - qy[i]);
  }

  free(indx_counter);
  free(xb);
  free(q_sum);
  free(q_sumsq);
  free(qy);
}

void qdiff_quick_and_dirty(struct two_sample *t) {

  // Approximate quantile difference statistic using conversion that:
  // length exact 95% CI for quantile is (assuming normal distribution)
  // approximately 4 times the SE
  //
  // Quantile function for binomial (from R source, in Rmath.h)
  // double qbinom(double p, double n, double pr, int lower_tail, int log_p)
  int nl, nu;
  double sx, sy;

  for (int i = 0; i < t->ntau; i++) {

    nl = (int)qbinom(0.025, t->nx, t->tau[i], 1, 0);
    nu = (int)qbinom(0.975, t->nx, t->tau[i], 1, 0);
    if (nu >= t->nx)
      nu = t->nx - 1;

    sx = ((t->x[nu] - t->x[nl]) / 4);

    nl = (int)qbinom(0.025, t->ny, t->tau[i], 1, 0);
    nu = (int)qbinom(0.975, t->ny, t->tau[i], 1, 0);
    if (nu >= t->ny)
      nu = t->ny - 1;
    sy = ((t->y[nu] - t->y[nl]) / 4);

    t->qvar[i] = (sx * sx + sy * sy);

    sx = quantile_srtd(t->x, t->nx, t->tau[i]);
    sy = quantile_srtd(t->y, t->ny, t->tau[i]);
    t->q[i] = sx - sy;
  }
}

void qdiff_approx(struct two_sample *t) {
  double qx, qy, fx, fy;

  for (int i = 0; i < t->ntau; i++) {
    qx = quantile_srtd(t->x, t->nx, t->tau[i]);
    qy = quantile_srtd(t->y, t->ny, t->tau[i]);
    t->q[i] = qx - qy;
    fx = density_pcd(qx, t->px);
    fy = density_pcd(qy, t->py);
    t->qvar[i] = (t->tau[i] * (1 - t->tau[i]) *
                  (1 / (t->nx * fx * fx) + 1 / (t->ny * fy * fy)));
  }
}


