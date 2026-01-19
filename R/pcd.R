sample_pcd <- function(n, knots, beta) {
  # void pcdR(int *n, double *knots, int *nknots, double *beta, double *s)
  h <- .C("pcdR", as.integer(n), as.double(knots), as.integer(length(knots)), as.double(beta), s = double(n))
  return(h$s)
}

qdiff_pcd <- function(x, y, tau, knotsx, betax, knotsy, betay) {
  # void pcdR(int *n, double *knots, int *nknots, double *beta, double *s)

  # void test_pcd(double *x, int *nx, double *y, int *ny, double *tau, int *ntau,
  #               double *knotsx, int *nknotsx, double *betax,
  #                double *knotsy, int *nknotsy, double *betay, int *table_length, double *q, double *qvar){

  table_length <- 100
  h <- .C("test_pcd", as.double(x), as.integer(length(x)), as.double(y), as.integer(length(y)),
    as.double(tau), as.integer(length(tau)), as.double(knotsx), as.integer(length(knotsx)), as.double(betax),
    as.double(knotsy), as.integer(length(knotsy)), as.double(betay), as.integer(table_length),
    q = double(length(tau)), qvar = double(length(tau))
  )

  list(qdiff = h$q, se = sqrt(h$qvar))
}

optim_pcd <- function(x, niter = 100, knots, beta = NULL) {
  if (is.null(beta)) {
    beta <- rep(1, length(knots) - 1)
  }

  h <- .C("optim_pcd", as.double(x), as.integer(length(x)), as.integer(niter),
    as.double(knots), as.integer(length(knots)), as.double(beta),
    llike = double(niter)
  )
  return(list(beta = h[[6]], knots = h[[4]], llike = h$llike))
}


power_pcd <- function(n, B, control = list()) {
  # method_1, method_2, crit_1, crit_2, tau, knotsx, betax, knotsy, betay){){

  m1 <- control$method_1
  m2 <- control$method_2
  c1 <- control$c1
  c2 <- control$c2
  tau <- control$tau
  ntau <- length(tau)
  knotsx <- control$knots
  knotsy <- control$knots
  betax <- control$betax
  betay <- control$betay
  table_length <- 100

  # void power_pcd(int *t1, double *c1, int *t2, double *c2, int *n, int *B, double *tau, int *ntau,
  #                double *knotsx, int *nknotsx, double *betax,
  #                double *knotsy, int *nknotsy, double *betay, int *table_length){

  h <- .C("power_pcd", as.integer(m1), as.double(c1), as.integer(m2), as.double(c2), as.integer(n), as.integer(B),
    as.double(tau), as.integer(length(tau)), as.double(knotsx), as.integer(length(knotsx)), as.double(betax),
    as.double(knotsy), as.integer(length(knotsy)), as.double(betay), as.integer(table_length),
    pwr = double(1)
  )

  return(h$pwr)
}

# void test_pcd(double *x, int *nx, double *y, int *ny, double *tau, int *ntau,
#               double *knotsx, int *nknotsx, double *betax,
#                double *knotsy, int *nknotsy, double *betay, int *table_length, double *q, double *qvar){
