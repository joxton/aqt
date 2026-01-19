

double quantile_srtd(double *x, int n, double tau);
int compare_doubles(const void *a, const void *b);

void permute(double *z, double *x, int n, double *y, int m);
void quantile_sorted(double *x, int *n, double *tau, int *ntau, double *q);
int sample_indx(int n);
