

struct pcd {
  // ---- piecewise_constant_density --------
  // Support of density is min(knots) to max(knots)
  // Between knots[k] and knots[k+1] the density is propotional to beta[k]

  double *knots;
  int nknots;
  double *betax;
  // Normalizing constant
  double nconst;
  double *cprob_at_knot;
  int *table;
  int table_length;
};

struct comp_pcd {

  struct pcd p1;
  struct pcd p2;
};

struct method {

  double *critical_value;
  double *tau;
  int ntau;
  int type;
  int combiner;
};

double density_pcd(double x, struct pcd *p);
void setup_pcd(double *knots, int nknots, double *beta, int table_length,
               struct pcd *p);
