qntl <- function(x, tau) {
  # quantile(double *x, int *n, double *tau, int *ntau, double *res){


  h <- .C("quantile", as.double(x), as.integer(length(x)), as.double(tau), as.integer(length(tau)), q = double(length(tau)))
  return(h$q)
}

sed_approx <- function(p, v, w, z = 1) {
  ## Approximate SE-diff
  v <- sort(v)
  w <- sort(w)
  nv <- length(v)
  nw <- length(w)
  tail_prob <- 1 - pnorm(z)
  x <- p
  nl <- qbinom(tail_prob, size = nv, prob = x)
  nu <- qbinom(1 - tail_prob, size = nv, prob = x)
  vl <- v[max(c(nl, 1))]
  vu <- v[min(c(nv, nu))]
  vse <- (vu - vl) / (2 * z)

  nl <- qbinom(tail_prob, size = nw, prob = x)
  nu <- qbinom(1 - tail_prob, size = nw, prob = x)
  wl <- w[max(c(nl, 1))]
  wu <- w[min(c(nw, nu))]
  wse <- (wu - wl) / (2 * z)

  se_diff <- sqrt(wse^2 + vse^2)
  return(se_diff)
}


qntl_qd <- function(x, y, tau, method = 1, R = 1000, delete_d = NULL, combiner = 1) {
  if (method == 4) {
    if (is.null(delete_d)) {
      delete_d <- round((length(x) + length(y)) / 2 * 0.5)
    }
  } else {
    delete_d <- 1
  }

  if (method == 1) {
    R <- 1
  }
  # void two_sample_aqt(double *x, int *n, double *y, int *m, double *tau, int *ntau, double *q, double *qvar){
  ntau <- length(tau)
  h <- .C("two_sample_aqt", as.double(x), as.integer(length(x)), as.double(y), as.integer(length(y)),
    as.double(tau), as.integer(ntau), as.integer(c(R, delete_d)), as.integer(c(method, combiner)),
    q = double(ntau * ifelse(method %in% c(1, 3, 4), 2, R)), qvar = double(ntau * ifelse(method %in% c(1, 3, 4), 1, R))
  )
  if (method %in% c(2)) {
    return(list(qdiff = matrix(h$q, nrow = R, byrow = TRUE), q2 = matrix(h$qvar, nrow = R, byrow = FALSE)))
  }
  if (method %in% c(1, 3, 4)) {
    return(list(qdiff = h$q[1:ntau], se = sqrt(h$qvar)))
  }
  if (method %in% c(5)) {
    return(list(qdiff_pvalues = h$q[1:ntau], pvalue = h$qvar[1]))
  }
}

qntl_c <- function(x, y, control = list()) {
  if (is.null(control$tau)) {
    n <- round((length(x) + length(y)) / 2)
    tau <- c(1:(n - 1)) / n
    control$tau <- tau
  }
  if (is.null(control$method)) {
    control$method <- "permutation"
  }
  if (is.null(control$combiner)) {
    control$combiner <- "afisher"
  }
  if (is.null(control$R)) {
    control$R <- 200
    if (control$method == "permutation") {
      control$R <- 1000
    }
  }

  ## Need change later ...
  control$statistic <- "pvalue"
  if (is.null(control$statistic)) {
    control$statistic <- "pvalue"
  }
  if (is.null(control$alpha)) {
    control$alpha <- 0.05
  }

  # Number of draws from asymptotic distribution
  if (is.null(control$B)) {
    control$B <- 10000
  }


  # Generate pvalue function, if not present
  if (control$method != "permutation" & control$combiner != "none") {
    if (is.null(control$cregion[[control$statistic]][[control$combiner]])) {
      if (is.null(control$cregion)) {
        control$cregion <- list()
      }
      if (is.null(control$cregion[[control$statistic]])) {
        control$cregion[[control$statistic]] <- list()
      }
      pvfun <- asymptotic_null(
        combiner = control$combiner, statistic = control$statistic,
        tau = control$tau, B = control$B, alpha = control$alpha
      )
      control$cregion[[control$statistic]][[control$combiner]] <- pvfun
    }
  }


  method_map <- list()
  method_map[["bootstrap"]] <- 3
  method_map[["jackknife"]] <- 4
  method_map[["qd"]] <- 1
  method_map[["permutation"]] <- 5

  combiner_map <- list()
  combiner_map[["minp"]] <- 1
  combiner_map[["afisher"]] <- 2
  combiner_map[["berkjones"]] <- 3

  q <- qntl_qd(
    x = x, y = y, tau = control$tau, method = method_map[[control$method]],
    R = control$R, delete_d = NULL, combiner = combiner_map[[control$combiner]]
  )
  return_list <- q
  if (method_map[[control$method]] < 5) {
    z <- q$qdiff / q$se
    function.pvalue <- control$cregion[[control$statistic]][[control$combiner]]
    return_list <- list(quantiles = data.frame(tau = control$tau, qdiff = q$qdiff, se = q$se), pvalue = function.pvalue(z), function.pvalue = function.pvalue)
  }
  return_list$control <- control
  return(return_list)
}


asymptotic_null <- function(combiner, statistic, tau, B = 10000, alpha = 0.05) {
  # ====================================================================================
  # Return a function for finding p-value of combination statistic
  # defined by "combiner" and "statistic"
  # p-value is based on monte-carlo approximation using asymptotic distribution
  # of sample quantiles (Brownian bridge)
  # ====================================================================================

  require(mvtnorm)
  ## Simulate Brownian bridge process evaluated at tau
  s <- sapply(tau, function(x) {
    g <- x * tau
    b <- x
    p <- sapply(tau, function(x) {
      return(min(b, x))
    })
    return(p - g)
  })
  d <- diag(s)
  m <- diag(1 / sqrt(d))
  sigma <- m %*% s %*% m
  draw_bridge <- rmvnorm(B, sigma = sigma)

  if (statistic == "zsquared") {
    draw_bridge <- (draw_bridge)^2
  }
  if (statistic == "pvalue") {
    draw_bridge <- -log(2 * (1 - pnorm(abs(draw_bridge))))
  }

  if (combiner == "minp") {
    cfun <- function(x) {
      return(max(x))
    }
  }

  if (combiner == "afisher") {
    cfun <- function(x) {
      return(cumsum(sort(x, decreasing = TRUE)))
    }
  }

  if (combiner == "berkjones") {
    cfun <- function(x) {
      return((sort(x, decreasing = TRUE)))
    }
  }

  cstat <- do.call(rbind, lapply(split(draw_bridge, 1:nrow(draw_bridge)), cfun))

  rank.statistics <- apply(cstat, 2, rank)
  p.value.statistics <- 1 - rank.statistics / B
  minp.combination <- apply(p.value.statistics, 1, min)
  critical.minp <- as.numeric(quantile(minp.combination, prob = alpha))
  rejection.region.statistic <- apply(cstat, 2, quantile, prob = 1 - critical.minp)

  pvalue_function <- function(x) {
    if (statistic == "pvalue") {
      x <- -log(2 * (1 - pnorm(abs(x))))
    }
    if (statistic == "zsquared") {
      x <- x^2
    }

    cf <- cfun(x)
    obs_statistic <- min(sapply(1:length(cf), function(k) {
      return(mean(cstat[, k] >= cf[k]))
    }))
    pvalue <- mean(obs_statistic >= minp.combination)

    return(pvalue)
  }
  return(pvalue_function)
}


qntl_boot <- function(x, B = 200, tau) {
  # h =.C("bootstrap_quantile", as.double(x), as.integer(length(x)), as.integer(B))
  # bootstrap_sample_quantile(double *x, int *n, int *B, double *tau, int *ntau, double *q, double *qvar);

  ntau <- length(tau)
  h <- .C("bootstrap_sample_quantile", as.double(x), as.integer(length(x)), as.integer(B),
    as.double(tau), as.integer(ntau),
    q = double(ntau), qvar = double(ntau)
  )

  return(list(quantile = h$q, se = sqrt(h$qvar)))
}


qntld_boot <- function(x, y, tau, B = 200) {
  # h =.C("bootstrap_quantile", as.double(x), as.integer(length(x)), as.integer(B))
  # bootstrap_sample_quantile(double *x, int *n, int *B, double *tau, int *ntau, double *q, double *qvar);
  #         bootstrap_quantile_difference(double *x, int *n, double *y, int *m, int *B, double *tau, int *ntau, double *q, double *qvar);

  ntau <- length(tau)
  h <- .C("bootstrap_quantile_difference", as.double(x), as.integer(length(x)), as.double(y), as.integer(length(y)), as.integer(B),
    as.double(tau), as.integer(ntau),
    q = double(ntau), qvar = double(ntau)
  )

  return(list(qdiff = h$q, se = sqrt(h$qvar)))
}


qntld_perm <- function(x, y, tau, B = 200, R = 1000, all = FALSE) {
  #  permutation_quantile_difference(double *x, int *n, double *y, int *m, int *B, double *tau, int *ntau, int *R, double *permdist);

  ntau <- length(tau)
  a <- ifelse(all, 1, 0)

  h <- .C("permutation_quantile_difference", as.double(x), as.integer(length(x)), as.double(y), as.integer(length(y)), as.integer(B),
    as.double(tau), as.integer(ntau), as.integer(R), as.integer(a),
    permdist = double(R * ifelse(all, length(tau), 1))
  )

  permdist <- h$permdist
  if (all) {
    permdist <- matrix(h$permdist, ncol = R, byrow = FALSE)
  }

  return(permdist)
}
