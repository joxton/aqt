.onAttach <- function(...) {
  if (!(exists(".atqEnv") && is.environment(atqEnv()))) {
    cat("Loading null-distribution data \n")

    ## ... adaptive quantile difference tests
    data(aqtNullDistribution, envir = .aqtEnv)

    ## ... adaptive quantile mood tests
    data(amtNullDistribution, envir = .aqtEnv)
  }
}


.aqtEnv <- new.env(parent = emptyenv())

aqtEnv <- function() .aqtEnv
