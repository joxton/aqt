## style_file("~/Desktop/mood/aqt/R/pcd_sampling.R")


get_range <- function(j, beta) {
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
  # Add range constraint on update to preserve stochastic dominance
  # of distribution defined by par[(npar/2 +1):npar] over density defined by
  # par[1:(npar/2)]
  # ==========================================================================

  #
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
