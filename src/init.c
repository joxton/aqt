#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
// htree_metal 
extern void quantile(double *x, int *n, double *tau, int *ntau, double *res); 
extern void bootstrap_quantile(double *x, int *n, int *B);
extern void bootstrap_sample_quantile(double *x, int *n, int *B, double *tau, int *ntau, double *q, double *qvar);
extern void bootstrap_quantile_difference(double *x, int *n, double *y, int *m, int *B, double *tau, int *ntau, double *q, double *qvar); 

extern void permutation_quantile_difference(double *x, int *n, double *y, int *m, int *B, double *tau, int *ntau, int *R, int *all, double *permdist);
extern void two_sample_aqt(double *x, int *n, double *y, int *m, double *tau, int *ntau, int *R,int *method, double *q, double *qvar); 
extern void optim_pcd(double *x, int *n, int *niter, double *knots, int *nknots, double *beta);
        

static const R_CMethodDef CEntries[] = {
  {"quantile",         (DL_FUNC) &quantile, 5},
  {"bootstrap_quantile",         (DL_FUNC) &bootstrap_quantile, 3},
  {"bootstrap_sample_quantile",         (DL_FUNC) &bootstrap_sample_quantile, 7},
  {"bootstrap_quantile_difference",         (DL_FUNC) &bootstrap_quantile_difference, 9},
  {"permutation_quantile_difference",         (DL_FUNC) &permutation_quantile_difference, 10},
  {"two_sample_aqt",         (DL_FUNC) &two_sample_aqt, 10},
  {"optim_pcd",         (DL_FUNC) &optim_pcd, 6},
    {NULL, NULL, 0}
};



void R_init_htree(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
