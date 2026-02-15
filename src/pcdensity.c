
#include "pcdensity.h"
#include "aqt.h"
#include <R.h>
#include <R_ext/Utils.h>
#include <Rmath.h>
#include <assert.h>
#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

double density_pcd(double x, struct pcd *p);
double sample_pcd_table(struct pcd *p);
double rejection_afisher(double *v, double *w, int n, struct comp_pcd *p,
                         struct method *m);
double rejection_minp(double *v, double *w, int n, struct comp_pcd *p,
                      struct method *m);
double rejection_ttest(double *v, double *w, int n, struct comp_pcd *p,
                       struct method *m);
double rejection(double *v, double *w, int n, struct comp_pcd *p,
                 struct method *m);
double power_difference_estimate(struct comp_pcd *p, int n, int B,
                                 struct method *m1, struct method *m2, double *pwr);

double llike(double *x, int n, struct pcd *p);
double density_pcd(double x, struct pcd *p);
void free_pcd(struct pcd *p);

void test_pcd(double *x, int *nx, double *y, int *ny, double *tau, int *ntau,
              double *knotsx, int *nknotsx, double *betax, double *knotsy,
              int *nknotsy, double *betay, int *table_length, double *q,
              double *qvar) {
  GetRNGstate();
  double qx, qy, dx, dy;
  struct comp_pcd p;
  setup_pcd(knotsx, *nknotsx, betax, *table_length, &(p.p1));
  setup_pcd(knotsy, *nknotsy, betay, *table_length, &(p.p2));

  for (int b = 0; b < 200; b++) {

    for (int i = 0; i < *nx; i++)
      x[i] = sample_pcd_table(&(p.p1));
    for (int i = 0; i < *ny; i++)
      y[i] = sample_pcd_table(&(p.p2));

    qsort(x, *nx, sizeof(double), compare_doubles);
    qsort(y, *ny, sizeof(double), compare_doubles);

    for (int i = 0; i < *ntau; i++) {

      qx = quantile_srtd(x, *nx, tau[i]);
      qy = quantile_srtd(y, *ny, tau[i]);
      q[i] = qx - qy; // se(q) = sqrt(p(1-p)/n) / f(q).

      dx = density_pcd(qx, &(p.p1));
      // Rprintf("dx = %lf x = %lf \n", dx, qx );
      dy = density_pcd(qy, &(p.p2));
      qvar[i] = tau[i] * (1 - tau[i]) *
                (1 / ((*nx) * dx * dx) + 1 / ((*ny) * dy * dy));
    }
  }
  PutRNGstate();
}

void normalizing_pcd(struct pcd *p) {

  // Normalizing constant
  p->nconst = 0;
  for (int i = 1; i < p->nknots; i++) {
    p->nconst += (p->knots[i] - p->knots[i - 1]) * p->betax[i - 1];
  }
}

void optim_pcd(double *x, int *n, int *niter, double *knots, int *nknots,
               double *beta, double *like) {

  struct pcd p;
  double delta = .1;
  double epsilon = 0.1;
  double l0, l1, d;
  double b0, b1;

  int table_length = 100;
  setup_pcd(knots, *nknots, beta, table_length, &p);

  for (int i = 0; i < *niter; i++) {

    for (int j = 0; j < (*nknots - 1); j++) {
      l0 = llike(x, *n, &p);

      p.betax[j] = p.betax[j] + delta;
      normalizing_pcd(&p);
      // p.nconst += (p.knots[j+1] - p.knots[j])*delta;
      l1 = llike(x, *n, &p);
      p.betax[j] = p.betax[j] - delta;
      normalizing_pcd(&p);
      // p.nconst += - (p.knots[j+1] - p.knots[j])*delta;
      d = (l1 - l0) / delta;

      b0 = p.betax[j];
      for (int k = 0; k < 10; k++) {
        // Half-stepping
        p.betax[j] = p.betax[j] + d * epsilon;
        if (p.betax[j] < 0)
          p.betax[j] = .00001;

        normalizing_pcd(&p);
        l1 = llike(x, *n, &p);
        if (l1 > l0) {
          break;
        } else {
          d = d / 2;
          p.betax[j] = b0;
          normalizing_pcd(&p);
        }
      }

      // p.nconst += (p.knots[j+1] - p.knots[j])*(b1 - b0);
      if (p.nconst < 0)
        Rprintf("Negative normalizing constant \n");
    }
    like[i] = l0;
    // Rprintf("ll[%d] = %lf \n", i, llike(x, *n, &p));
  }

  free_pcd(&p);
}

double llike(double *x, int n, struct pcd *p) {

  double ll, d;
  ll = 0;

  for (int i = 0; i < n; i++) {
    // Rprintf(" x[%d] = %lf d = %lf", i, x[i], density_pcd(x[i], p));
    d = density_pcd(x[i], p);
    if (d < 0)
      Rprintf("Negative density\n");
    ll += log(d);
  }
  ll = ll / ((double)n);
  return (ll);
}

double density_pcd_p(double x, struct pcd *p) {
  int upper_indx, lower_indx, mid_indx;
  double *knots;
  double mid_point;

  knots = p->knots;
  upper_indx = p->nknots - 1;
  lower_indx = 0;

  while ((upper_indx - lower_indx) > 1) {
    mid_indx = lower_indx + ((int)(upper_indx - lower_indx) / 2);
    mid_point = knots[mid_indx];
    if (x >= mid_point) {
      lower_indx = mid_indx;
    } else {
      upper_indx = mid_indx;
    }
  }
  return (p->betax[lower_indx] / p->nconst);
}


double density_pcd(double x, struct pcd *p) {
  int upper_indx, lower_indx, mid_indx;
  double *knots;
  double mid_point;

  knots = p->knots;
  upper_indx = p->nknots - 1;
  lower_indx = 0;
	for(int i = 0; i < (p->nknots - 1); i++)
	{
		if(x >= knots[i] && x < knots[i+1]){
			lower_indx = i; 	
			break;
		}
	}

  return (p->betax[lower_indx] / p->nconst);
}


void setup_pcd(double *knots, int nknots, double *betax, int table_length,
               struct pcd *p) {

  p->knots = knots;
  p->nknots = nknots;
  p->betax = betax;
  p->cprob_at_knot = (double *)malloc(sizeof(double) * nknots);

  // Normalizing constant
  p->nconst = 0;
  for (int i = 1; i < nknots; i++) {
    p->nconst += (p->knots[i] - p->knots[i - 1]) * p->betax[i - 1];
  }

  p->cprob_at_knot[0] = 0;
  for (int i = 1; i < nknots; i++) {
    p->cprob_at_knot[i] =
        p->cprob_at_knot[i - 1] +
        (p->knots[i] - p->knots[i - 1]) * p->betax[i - 1] / p->nconst;
  }

  // Look-up table
  p->table_length = table_length;
  p->table = (int *)malloc(sizeof(double) * table_length);

  int cknot = 0;
  double cprob = 0;

  p->table[0] = 0;
  for (int i = 1; i < p->table_length; i++) {
    cprob = (double)((double)i) / ((double)table_length);
    if (cprob < p->cprob_at_knot[cknot + 1]) {
      p->table[i] = cknot;
    } else {
      while (cprob >= p->cprob_at_knot[cknot + 1]) {
        cknot++;
      }
      p->table[i] = cknot;
    }
  }
}

void free_pcd(struct pcd *p) {
  free(p->table);
  free(p->cprob_at_knot);
}

double sample_pcd_table(struct pcd *p) {

  int indx;
  int cknot;
  double u, dp, du, dx, s;
  int counter;
  // counter = 0;

  u = unif_rand();
  indx = (int)floor(u * ((double)p->table_length));
  // Rprintf("u = %lf indx = %d tablen = %d \n", u, indx, p->table_length);

  cknot = p->table[indx];
  // Rprintf("p->tableindx = %d \n", p->table[indx]);

  while (p->cprob_at_knot[cknot + 1] < u) {
    cknot++;
    // counter++;
  }

  // Rprintf("Search counter = %d\n", counter);
  dp = (p->cprob_at_knot[cknot + 1] - p->cprob_at_knot[cknot]);
  du = u - p->cprob_at_knot[cknot];
  dx = p->knots[cknot + 1] - p->knots[cknot];
  s = du / dp * dx + p->knots[cknot];

  return (s);
}

double sample_pcd(struct pcd *p) {
  // Placeholder versjon
  // Could have look-up table to start search: table[floor(u*table_length)]
  // ... gives knot such that cprob_at_knot[knot]<=
  // floor(u*table_length)/table_length

  double u, dp, du, dx, s;

  u = unif_rand();

  for (int i = 0; i < (p->nknots - 1); i++) {
    if (u >= p->cprob_at_knot[i] && u < p->cprob_at_knot[i + 1]) {
      dp = (p->cprob_at_knot[i + 1] - p->cprob_at_knot[i]);
      du = u - p->cprob_at_knot[i];
      dx = p->knots[i + 1] - p->knots[i];
      s = du / dp * dx + p->knots[i];
      break;
    }
  }
  return (s);
}

void pcdR(int *n, double *knots, int *nknots, double *betax, double *s) {
  GetRNGstate();
  struct pcd p;
  int table_length = 1000;

  setup_pcd(knots, *nknots, betax, table_length, &p);

  for (int i = 0; i < *n; i++) {
    s[i] = sample_pcd_table(&p); // sample_pcd(&p); //
  }

  free_pcd(&p);
  PutRNGstate();
}

void setup_method(int type, double *tau, int ntau, double *critical_value,
                  struct method *m) {

  m->type = type;
  m->tau = tau;
  m->ntau = ntau;
  m->critical_value = critical_value;
}

void power_pcd(int *t1, double *c1, int *t2, double *c2, int *n, int *B,
               double *tau, int *ntau, double *knotsx, int *nknotsx,
               double *betax, double *knotsy, int *nknotsy, double *betay,
               int *table_length, double *pwr) {
  GetRNGstate();
  double pd;
  struct comp_pcd p;
  struct method m1;
  struct method m2;

  setup_pcd(knotsx, *nknotsx, betax, *table_length, &(p.p1));
  setup_pcd(knotsy, *nknotsy, betay, *table_length, &(p.p2));

  setup_method(*t1, tau, *ntau, c1, &m1);
  setup_method(*t2, tau, *ntau, c2, &m2);

  pwr[0] = power_difference_estimate(&p, *n, *B, &m1, &m2, &(pwr[1]) );
  // Rprintf("power = %lf \n", pd);

  free_pcd(&(p.p1));
  free_pcd(&(p.p2));

  PutRNGstate();
}

double power_difference_estimate(struct comp_pcd *p, int n, int B,
                                 struct method *m1, struct method *m2, double *pwr) {

  double *v, *w;
  double r1, r2;

  v = (double *)malloc(sizeof(double) * n);
  w = (double *)malloc(sizeof(double) * n);

  int nreject = 0;

  r1 = 0;
  r2 = 0;
  for (int b = 0; b < B; b++) {
    for (int i = 0; i < n; i++) {
      v[i] = sample_pcd_table(&(p->p1));
      w[i] = sample_pcd_table(&(p->p2));
    }

    qsort(v, n, sizeof(double), compare_doubles);
    qsort(w, n, sizeof(double), compare_doubles);

    r1 += rejection(v, w, n, p, m1);
    r2 += rejection(v, w, n, p, m2);

    // nreject += test_wrapper(v, w, method, critical_value);
  }

  pwr[0] = r1/((double)B);
  pwr[1] = r2/((double)B); 

  return ((r1 - r2) / ((double)B));
}

double rejection(double *v, double *w, int n, struct comp_pcd *p,
                 struct method *m) {

  double r;

  switch (m->type) {

  case 1:
    r = rejection_minp(v, w, n, p, m);
    break;

  case 2:
    r = rejection_afisher(v, w, n, p, m);
    break;

  case 3:
    r = rejection_ttest(v, w, n, p, m);
    break;

  case -1:
    r = 0;
    break;
  }

  return (r);
}

double rejection_ttest(double *v, double *w, int n, struct comp_pcd *p,
                       struct method *m) {

  double crit, tstat, dgfr, s;
  double a1, s1, a2, s2;
  a1 = a2 = s1 = s2 = 0;
  for (int i = 0; i < n; i++) {
    a1 += v[i];
    s1 += v[i] * v[i];
    a2 += w[i];
    s2 += w[i] * w[i];
  }
  a1 = a1 / ((double)n);
  a2 = a2 / ((double)n);
  s1 = s1 / ((double)n) - a1 * a1;
  s2 = s2 / ((double)n) - a2 * a2;

  s1 = ((double)n) / ((double)n - 1) * s1;
  s2 = ((double)n) / ((double)n - 1) * s2;

  s1 = sqrt(s1 / ((double)n));
  s2 = sqrt(s2 / ((double)n));

  s = sqrt(s1 * s1 + s2 * s2);
  dgfr = pow(s, 4.0) /
         (pow(s1, 4.0) / ((double)n - 1) + pow(s2, 4.0) / ((double)n - 1));
  tstat = fabs((a1 - a2) / s);
  crit = qt(0.975, dgfr, 1, 0);

  return ((tstat > crit) ? 1 : 0);
}

double rejection_minp(double *v, double *w, int n, struct comp_pcd *p,
                      struct method *m) {

  double qx, qy, dx, dy, q, qvar;
  double reject;

  reject = 0.0;
  for (int i = 0; i < m->ntau; i++) {

    qx = quantile_srtd(v, n, m->tau[i]);
    qy = quantile_srtd(w, n, m->tau[i]);
    q = qx - qy; // se(q) = sqrt(p(1-p)/n) / f(q).

    dx = density_pcd(qx, &(p->p1));
    // Rprintf("dx = %lf x = %lf \n", dx, qx );
    dy = density_pcd(qy, &(p->p2));
    qvar =
        m->tau[i] * (1 - m->tau[i]) * (1 / ((n)*dx * dx) + 1 / ((n)*dy * dy));
    q = q * q / qvar;
    if (m->critical_value[0] < (q)) {
      reject = 1;
      break;
    }
  }
  return (reject);
}

double rejection_afisher(double *v, double *w, int n, struct comp_pcd *p,
                         struct method *m) {

  double qx, qy, dx, dy, q, qvar;
  double reject;
  
  double *qsum; 
  qsum = (double*)malloc(sizeof(double)*m->ntau); 

 
  for (int i = 0; i < m->ntau; i++) {

    qx = quantile_srtd(v, n, m->tau[i]);
    qy = quantile_srtd(w, n, m->tau[i]);
    q = qx - qy; // se(q) = sqrt(p(1-p)/n) / f(q).

    dx = density_pcd(qx, &(p->p1));
    // Rprintf("dx = %lf x = %lf \n", dx, qx );
    dy = density_pcd(qy, &(p->p2));
    qvar =
        m->tau[i] * (1 - m->tau[i]) * (1 / ((n)*dx * dx) + 1 / ((n)*dy * dy));
    q = q * q / qvar;
    qsum[i] = q; 
  } 
   
   // Sort in descending order
    qsort(qsum, m->ntau, sizeof(double), compare_doubles_descending);
  reject = 0.0; 
  q = 0;  
    for (int i = 0; i < m->ntau; i++) {
    q += qsum[i]; 
    if (m->critical_value[i] < (q)) {
      reject = 1;
      break;
    }
  }
  free(qsum); 
  return (reject);
}



