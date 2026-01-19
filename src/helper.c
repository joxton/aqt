#include <R.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

void revsorts(double *a, int n) {
  /* Sort a[] into descending order by "heapsort";
   */

  int l, j, ir, i;
  double ra;
  int ii;

  if (n <= 1)
    return;

  a--;

  l = (n >> 1) + 1;
  ir = n;

  for (;;) {
    if (l > 1) {
      l = l - 1;
      ra = a[l];
    } else {
      ra = a[ir];
      a[ir] = a[1];

      if (--ir == 1) {
        a[1] = ra;
        return;
      }
    }
    i = l;
    j = l << 1;
    while (j <= ir) {
      if (j < ir && a[j] > a[j + 1])
        ++j;
      if (ra > a[j]) {
        a[i] = a[j];
        j += (i = j);
      } else
        j = ir + 1;
    }
    a[i] = ra;
  }
}

void sort_descending(double *a, int n) { revsorts(a, n); }

void revsortx(double *a, int *ib, int n) {
  /* Sort a[] into descending order by "heapsort";
   * sort ib[] alongside;
   * if initially, ib[] = 1...n, it will contain the permutation finally
   */

  int l, j, ir, i;
  double ra;
  int ii;

  if (n <= 1)
    return;

  a--;
  ib--;

  l = (n >> 1) + 1;
  ir = n;

  for (;;) {
    if (l > 1) {
      l = l - 1;
      ra = a[l];
      ii = ib[l];
    } else {
      ra = a[ir];
      ii = ib[ir];
      a[ir] = a[1];
      ib[ir] = ib[1];
      if (--ir == 1) {
        a[1] = ra;
        ib[1] = ii;
        return;
      }
    }
    i = l;
    j = l << 1;
    while (j <= ir) {
      if (j < ir && a[j] > a[j + 1])
        ++j;
      if (ra > a[j]) {
        a[i] = a[j];
        ib[i] = ib[j];
        j += (i = j);
      } else
        j = ir + 1;
    }
    a[i] = ra;
    ib[i] = ii;
  }
}

void sort_descending_with_order(double *a, int *ib, int n) {
  revsortx(a, ib, n);
}

void sort_test(int *niter, double *a, int *n, double *res, int *b) {

  // void revsort(double *a, int *ib, int n)

  double *ax;
  int *ib;
  int k, i;
  double x;

  ib = (int *)malloc(sizeof(int) * n[0]);
  ax = (double *)malloc(sizeof(double) * n[0]);

  // for(k = 0; k < n[0]; k++)
  //	ib[k] = 0;

  //	Rprintf(" n = %d \n", n[0]);
  for (k = 0; k < niter[0]; k++) {

    for (i = 0; i < n[0]; i++) {
      // ax[i]=a[i];
      // ib[i]=i;
      x = (double)k * i;
      if (a[i] < x) {
        ib[i] = 1;
      } else {
        ib[i] = 0;
      }
    }
    //  revsort(ax,ib,n[0]);
  }

  for (i = 0; i < n[0]; i++) {
    res[i] = ax[i];
    b[i] = ib[i];
  }
}
