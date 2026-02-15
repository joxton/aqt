get_range <- function(j, beta) {
  # ==========================================================================
  # Return sampling range for beta[j] such that unimodality of distribution is preserved 
  # ==========================================================================

  j <- j
  min_j <- NA
  max_j <- NA
  if (j == 1 | j == length(beta)) {
    if (j == 1) {
      min_j <- -beta[1]
      max_j <- (beta[2] - beta[1])
    }

    if (j == length(beta)) {
      min_j <- -beta[j]
      max_j <- (beta[j - 1] - beta[j])
    }
  } else {
    if (beta[j - 1] >= beta[j] & beta[j + 1] <= beta[j]) {
      min_j <- (beta[j + 1] - beta[j])
      max_j <- (beta[j - 1] - beta[j])
    }


    if (beta[j - 1] <= beta[j] & beta[j + 1] >= beta[j]) {
      min_j <- (beta[j - 1] - beta[j])
      max_j <- (beta[j + 1] - beta[j])
    }

    if (beta[j - 1] < beta[j] & beta[j + 1] < beta[j]) {
      if (beta[j - 1] >= beta[j + 1]) {
        min_j <- (beta[j + 1] - beta[j])
        max_j <- 1
      }
      if (beta[j - 1] < beta[j + 1]) {
        min_j <- (beta[j - 1] - beta[j])
        max_j <- 1
      }
    }
  }
  return(list(min_j = min_j, max_j = max_j))
}


get_range_stdom <- function(j, betax, betay) {
  # ==========================================================================
  # Return sampling range for beta[j] such that stochastic dominance of distribution 
  # (proportional to) betay over distribution (proportional to betax is preserved
  # (j > length(betax) => updates betay
  # ==========================================================================

  
  nbeta <- length(betax)

  cbetax <- cumsum(betax)
  cbetay <- cumsum(betay)

  ## cumulative distributions
  dbetax <- cbetax / cbetax[nbeta]
  dbetay <- cbetay / cbetay[nbeta]


  if (j <= nbeta) {
    ## Sampling for the dominated distribution

    # Get range constraining unimodality
    rang <- get_range(j = j, beta = betax)

    lower_bound <- NULL
    if (j > 1) {
      lower_bound <- max(cbetax[1:(j - 1)] / dbetay[1:(j - 1)] - cbetax[nbeta])
    }

    upper_bound <- NULL
    if (j < nbeta) {
      upper_bound <- min((dbetay[j:(nbeta - 1)] * cbetax[nbeta] - cbetax[j:(nbeta - 1)]) / (1 - dbetay[j:(nbeta - 1)]))
    }

    ## max_sd_constraint = (cbetax[nbeta]/cbetay[nbeta])*cbetay[j] - ifelse(j == 1, 0, cbetax[j-1])
    rang$max_j <- min(c(rang$max_j, upper_bound))
    rang$min_j <- max(c(rang$min_j, lower_bound))
  } else {
    ## Sampling for the dominating distribution

    # Get range constraining unimodality
    jstar <- j - nbeta
    rang <- get_range(j = jstar, beta = betay)

    upper_bound <- NULL
    if (jstar > 1) {
      upper_bound <- min(cbetay[1:(jstar - 1)] / dbetax[1:(jstar - 1)] - cbetay[nbeta])
    }

    lower_bound <- NULL
    if (jstar < nbeta) {
      lower_bound <- max((dbetax[jstar:(nbeta - 1)] * cbetay[nbeta] - cbetay[jstar:(nbeta - 1)]) / (1 - dbetax[jstar:(nbeta - 1)]))
    }


    rang$max_j <- min(c(rang$max_j, upper_bound))
    rang$min_j <- max(c(rang$min_j, lower_bound))
  }

  return(rang)
}


power_pcdw <- function(n, B, control) {
  # ==========================================================================
  # Return estimated power difference between control$test1 and control$test2 
  # on data from piecewise-constant-densities (pcd) control$betax and control$betay  
  # ==========================================================================

  v <- NULL
  w <- NULL
  if (!is.null(control$seed)) {
    set.seed(control$seed)
  }
  for (k in 1:B) {
    x <- sample_pcd(n = n, knots = control$knotsx, beta = control$betax)
    y <- sample_pcd(n = n, knots = control$knotsy, beta = control$betay)

    h <- control$test1(x, y) 
    v <- c(v, h) 
    h <- control$test2(x, y) 
    w <- c(w, h)
  }

  mv <- mean(v <= .05)
  mw <- mean(w <= .05)

  return(c(mv - mw, mv, mw))
}



pwr_optim <- function(control) {
  # ==========================================================================
  # Maximize power of two-sample test control$test1 over test control$test2 
  # over piecewise-constant-densitites. 
  # ==========================================================================
  if(class(control) != "list") 
	stop("control is not a list.") 
  if(class(control$test1) != "function") 
	stop("control$test1 is not a function.") 
  if(class(control$test2) != "function") 
	stop("control$test2 is not a function.") 

  
  knots <-  seq(-4, 4, length = 10)
  beta = c(0.00402, 0.07901, 0.57987, 1.42218, 1.88708, 1.45477, 0.57430, 0.06555, 0.00483) 

  control.default <- list(
		max_doubling_stagnant = 2,	## If two successive knot-doublings don't increase power difference, terminate
		max_doubling = 5,		## Maximum number of knot doublings	     	 
		max_stagnant = 5,		## Maximum number of iterations without improvement (for given knots setting)  
		max_iter = 100, 		## Maximum number of iterations (for given knots setting) 
		n = 200, 			## Sample sizes for the tests
		B = 200, 			## Number of runs for each power estimation
		controlB = 500, 		## Number of runs for power estimation in objective function 
		controlseed = 898989, 		## seed for objective function 
  		verbose = TRUE, 		## 
		location_shift = TRUE, 		## Location shift model? 
		constraint = "unimodal", 	## density constraint, can be "none", "unimodal", "dominance" 
		knotsx = knots, 		## Knots for density of X 
		knotsy = knots, 		## Knots for density of Y 
		betax = beta, 			## Density heights between knots for X
		betay = beta  			## Density heights between knots for Y
		)

  for(k in names(control.default)){
	if(is.null(control[[k]] ))
		control[[k]] = control.default[[k]] 
  }


  pwr_diff <- 0
  stagnant_counter <- 0
  doubling_counter <- 0
  k <- 1
  while (TRUE) {
    if (k > 1) {
      doubling_counter <- doubling_counter + 1
      cat("Knot doubling\n")
      control <- knot_doubler(control = control)
    }

    if (control$location_shift) {
      control <- cm_powerdiff_lshift(control = control)
    } else {
      control <- cm_powerdiffw(control = control)
    }

    if (control$pwr_diff <= pwr_diff) {
      stagnant_counter <- stagnant_counter + 1
    } else {
      stagnant_counter <- 0
    }
    pwr_diff <- control$pwr_diff

    if (doubling_counter >= control$max_doubling | stagnant_counter >= control$max_doubling_stagnant) {
      break
    }
    k <- k + 1
  }
  return(control)
}

knot_doubler <- function(control) {
  nknots <- length(control$knotsx)
  control$knotsx <- sort(c(control$knotsx, (control$knotsx[2:nknots] + control$knotsx[1:(nknots - 1)]) / 2))
  control$knotsy <- sort(c(control$knotsy, (control$knotsy[2:nknots] + control$knotsy[1:(nknots - 1)]) / 2))
  control$betax <- control$betax[sort(rep(1:(nknots - 1), 2))]
  control$betay <- control$betay[sort(rep(1:(nknots - 1), 2))]
  return(control)
}

cm_powerdiffw <- function(control) {
  pwr_est <- ifelse(is.null(control$pwr_diff), 0, control$pwr_diff)

  counter_stagnant <- 0
  iter <- 0
  n <- control$n
  B <- control$B
  while (TRUE) {
    iter <- iter + 1
    set.seed(as.integer(Sys.time()))
    cprev <- control
    p0 <- power_pcdw(n, B, control = control)[1]
    for (j in 1:length(control$betax)) {
      b0 <- control$betax[j]
      if (control$constraint == "none") {
        r <- list(min_j = -b0, max_j = 0.5)
      }
      if (control$constraint == "unimodal") {
        r <- get_range(j, control$betax)
      }
      if (control$constraint == "dominance") {
        r <- get_range_stdom(j, control$betax, control$betay)
      }
      if (!(is.na(r$min_j) | is.na(r$max_j))) {
        if (r$min_j < r$max_j) {
          control$betax[j] <- b0 + runif(1, min = r$min_j, max = r$max_j)
          p1 <- power_pcdw(n, B, control = control)[1]
          if (p1 < p0) {
            control$betax[j] <- b0
          } else {
            p0 <- p1
          }
        }
      }
    }
    for (j in 1:length(control$betay)) {
      b0 <- control$betay[j]
      if (control$constraint == "none") {
        r <- list(min_j = -b0, max_j = 0.5)
      }
      if (control$constraint == "unimodal") {
        r <- get_range(j, control$betay)
      }
      if (control$constraint == "dominance") {
        r <- get_range_stdom(j + length(control$betax), control$betax, control$betay)
      }
      if (!(is.na(r$min_j) | is.na(r$max_j))) {
        if (r$min_j < r$max_j) {
          control$betay[j] <- b0 + runif(1, min = r$min_j, max = r$max_j)
          p1 <- power_pcdw(n, B, control = control)[1]
          if (p1 < p0) {
            control$betay[j] <- b0
          } else {
            p0 <- p1
          }
        }
      }
    }

    cc <- control
    cc$seed <- control$controlseed
    cc$R <- control$controlR
    pwr_new <- power_pcdw(n, B = control$controlB, control = cc)[1]
    if (pwr_new <= pwr_est) {
      control <- cprev
      counter_stagnant <- counter_stagnant + 1
    } else {
      pwr_est <- pwr_new
      counter_stagnant <- 0
    }

    if (counter_stagnant > control$max_stagnant | iter > control$max_iter) {
      break
    }
    if (control$verbose) {
      cat(" Iteration ", iter, " power : ", pwr_est, "\n")
    }
  }
  control$pwr_diff <- pwr_est
  return(control)
}


cm_powerdiff_lshift <- function(control) {
  # Power difference optimization under location shift model

  pwr_est <- ifelse(is.null(control$pwr_diff), 0, control$pwr_diff)
  counter_stagnant <- 0
  iter <- 0
  n <- control$n
  B <- control$B
  while (TRUE) {
    iter <- iter + 1
    set.seed(as.integer(Sys.time()))
    cprev <- control
    p0 <- power_pcdw(n, B, control = control)[1]

    ## Sample location shift
    for (j in 1:min(c(10, length(control$betax)))) {
      knotsy <- control$knotsy
      control$knotsy <- control$knotsy + runif(1, min = -0.5, max = 0.5)
      p1 <- power_pcdw(n, B, control = control)[1]

      if (p1 < p0) {
        control$knotsy <- knotsy
      } else {
        p0 <- p1
      }
    }


    for (j in 1:length(control$betax)) {
      b0 <- control$betax[j]
      if (control$constraint == "none") {
        r <- list(min_j = -b0, max_j = 0.5)
      }
      if (control$constraint == "unimodal") {
        r <- get_range(j, control$betax)
      }
      if (!(is.na(r$min_j) | is.na(r$max_j))) {
        if (r$min_j < r$max_j) {
          control$betax[j] <- b0 + runif(1, min = r$min_j, max = r$max_j)
          control$betay[j] <- control$betax[j]
          p1 <- power_pcdw(n, B, control = control)[1]

          if (p1 < p0) {
            control$betax[j] <- b0
            control$betay[j] <- control$betax[j]
          } else {
            p0 <- p1
          }
        }
      }
    }

    cc <- control
    cc$seed <- control$controlseed
    cc$R <- control$controlR
    pwr_new <- power_pcdw(n, B = control$controlB, control = cc)[1]
    if (pwr_new <= pwr_est) {
      control <- cprev
      counter_stagnant <- counter_stagnant + 1
    } else {
      pwr_est <- pwr_new
      counter_stagnant <- 0
    }

    if (counter_stagnant > control$max_stagnant | iter > control$max_iter) {
      break
    }
    if (control$verbose) {
      cat(" Iteration ", iter, " power : ", pwr_est, "\n")
    }
  }

  control$pwr_diff <- pwr_est
  return(control)
}
