# ===========================================================================================================
# Functions for sampling/optimizing piecewise-constant-densities (pcd) 
# ===========================================================================================================

sample_pcd <- function(n, knots, beta) {
  h <- .C("pcdR", as.integer(n), as.double(knots), as.integer(length(knots)), as.double(beta), s = double(n))
  return(h$s)
}

qdiff_pcd <- function(x, y, tau, knotsx, betax, knotsy, betay) {
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
  knotsx <- control$knotsx
  knotsy <- control$knotsy
  betax <- control$betax
  betay <- control$betay
  table_length <- 100

  # void power_pcd(int *t1, double *c1, int *t2, double *c2, int *n, int *B, double *tau, int *ntau,
  #                double *knotsx, int *nknotsx, double *betax,
  #                double *knotsy, int *nknotsy, double *betay, int *table_length){

  h <- .C("power_pcd", as.integer(m1), as.double(c1), as.integer(m2), as.double(c2), as.integer(n), as.integer(B),
    as.double(tau), as.integer(length(tau)), as.double(knotsx), as.integer(length(knotsx)), as.double(betax),
    as.double(knotsy), as.integer(length(knotsy)), as.double(betay), as.integer(table_length),
    pwr = double(3)
  )

  return(h$pwr)
}


# ====================================================================================================================
# Critical value for asymptotic statistic
# ====================================================================================================================
asymptotic_test_normal <- function(x, y, tau) {
  qx <- as.numeric(quantile(x, prob = tau))
  qy <- as.numeric(quantile(y, prob = tau))
  dx <- dnorm(qx)
  dy <- dnorm(qy)

  vari <- tau * (1 - tau) * (1 / (length(x) * dx * dx) + 1 / (length(y) * dy * dy))
  z <- sort((qx - qy)^2 / vari, decreasing = TRUE)
  stat_afisher <- cumsum(z)
  stat_minp <- z[1]
  stat_berkjones <- z
  list(afisher = stat_afisher, minp = stat_minp, berkjones = stat_berkjones)
}

null_distribution_normal_approx <- function(B = 10000, n = 200, tau = c(1:99) / 100) {
  # n = 200; B = 10000; tau = c(1:99)/100

  h <- lapply(1:B, function(b) {
    return(asymptotic_test_normal(x = rnorm(n), y = rnorm(n), tau = tau)$afisher)
  })
  afish <- do.call(rbind, h)

  h <- lapply(1:B, function(b) {
    return(asymptotic_test_normal(x = rnorm(n), y = rnorm(n), tau = tau)$minp)
  })
  minp <- do.call(rbind, h)

  h <- lapply(1:B, function(b) {
    return(asymptotic_test_normal(x = rnorm(n), y = rnorm(n), tau = tau)$berkjones)
  })
  berkj <- do.call(rbind, h)


  a <- (B - apply(berkj, 2, rank)) / B
  ac <- quantile(apply(a, 1, min), prob = 0.05)
  berkj_critical <- apply(berkj, 2, quantile, prob = 1 - ac)

  a <- (B - apply(afish, 2, rank)) / B
  ac <- quantile(apply(a, 1, min), prob = 0.05)
  afish_critical <- apply(afish, 2, quantile, prob = 1 - ac)

  minp_critical <- as.numeric(quantile(minp, prob = .95))

  list(minp_critical = minp_critical, berkj_critical = berkj_critical, afish_critical = afish_critical)
}
