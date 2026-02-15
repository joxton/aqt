\name{qd.test}
\alias{qd.test}
\title{Adaptive quantile-based two-sample tests}
\description{Performs adaptive quantile-based two-sample tests.}
\usage{

qd.test(x, y, control = list()) 

}
\arguments{
\item{x}{A numeric vector of data values}
\item{y}{A numeric vector of data values}
\item{control}{A list of arguments for controlling the test (see below).}
}

\details{

 \code{qd.test} performs adaptive quantile-based two-sample tests. The sample quantiles of  \code{x} and \code{y} are compared across the quantiles \code{control$tau}, and a combined test is formed pooling the individual 
comparisons. The type of combination is specified by \code{control$combiner} which can be \code{"minp"}, \code{"berkjones"} and \code{"afisher"}.  
The  \code{"minp"} variant combines the individual quantile comparisons using the one whose squared standardized difference is the greatest (like the 'minp' method for p-value combination); the \code{"berkjones"} variant 
forms order statistics of the individual comparisons and uses the most extreme order statistic as the summary test. The \code{"afisher"} combiner uses an approach similar to the adaptive Fisher p-value combination (Yu et.al. (2009)), basing the 
test on the most extreme element of the cumulative sum of ordered (from largest to smallest) squared standardized quantile differences.  

Quantile difference and rank based quantile comparisons can be selected. Quantile difference tests are obtained by setting \code{control$method} equal to \code{"bootstrap"} (or to \code{"bootp"} for parallelized estimation), while \code{control$method = "mood"} gives rank based comparisons .

The quantile difference based tests are as follows. 
Let  \eqn{ \hat{q}_X(\tau)} be the \eqn{\tau}-th sample quantile of \eqn{X}, 
  \eqn{ \hat{q}_Y(\tau)} for \eqn{Y}, and \eqn{\Delta(\tau) =  ( \hat{q}_X(\tau) - \hat{q}_Y(\tau)) }. The  tests are based on the 
  standardized differences \eqn{Z(\tau) = \Delta(\tau) /\hat{se}(\Delta(\tau)) }. The \code{"minp"} combination forms a combined test-statistic 
  using \eqn{\max_{\tau} = Z^2(\tau) }; the \code{"berkjones"} is based on the ordered standardized differences, i.e \eqn{ Z_{(1)}^2 \geq Z_{(2)}^2 \geq \ldots  
  \geq Z_{(m)}^2}, using the ordered difference with the smallest p-value as the combined test statistic; the \code{"afisher"} statistic is based on 
  the sums \eqn{ \sum_{j = 1}^k  Z_{(j)}^2} for \code{j = 1,...,m},  and use the smallest of their respective p-values as the combined test statistic.  
  
The rank based quantile tests use the following. Given sample values of \eqn{X} and \eqn{Y}, testing for a difference in their  \eqn{\tau}-th quantile can be done by finding the \eqn{\tau}-th quantile of the their combined sample and testing for a difference between the proportions of  \eqn{X}s and \eqn{Y}s falling above/below the combined sample quantile. When \eqn{\tau = 0.5} this is Mood's median test. Applying this test across the range of \eqn{\tau} values, combined tests as above can be formed. 

Finally, setting \code{control$method = "permutation"}, permutation is used to compare the individual quantiles and form the null distributions of the 
combined tests. The number of bootstrap/permutation samples are set by \code{control$R}. 




} 


\value{ Returns  \code{quantiles}, a \code{data.frame} with the quantile differences and their estimated standard errors, and \code{pvalue}, the p-value of the combined test. }  

\references{
Yu, K., Li, Q., Bergen, A.W., Pfeiffer, R.M., Rosenberg, P.S., Caporaso, N., Kraft, P., Chatterjee, N. (2009) \dQuote{Pathway analysis by adaptive combination of p-values.} \emph{ Genetic Epidemiology } 33, 700–709.

}

\author{Joe Sexton \email{joesexton0@gmail.com}}

\seealso{ }

\examples{

\dontrun{


  x <- rnorm(200)
  y <- rnorm(200, mean = .3)


# -- Rank based quantile tests 
  qd.test(x, y, control = list(R = 200, method = "mood", combiner = "minp"))$pvalue
  qd.test(x, y, control = list(R = 200, method = "mood", combiner = "afisher"))$pvalue
  # --- not implemented ---#  qd.test(x, y, control = list(R = 200, method = "mood", combiner = "berkjones"))$pvalue

# -- Quantile difference tests with bootstrap standard errors
  qd.test(x, y, control = list(R = 200, method = "bootstrap", combiner = "minp"))$pvalue
  qd.test(x, y, control = list(R = 200, method = "bootstrap", combiner = "afisher"))$pvalue
  qd.test(x, y, control = list(R = 200, method = "bootstrap", combiner = "berkjones"))$pvalue

# -- Quantile difference tests with bootstrap standard errors, bootstrapping in parallel
  qd.test(x, y, control = list(R = 200, method = "bootp", combiner = "minp"))$pvalue
  qd.test(x, y, control = list(R = 200, method = "bootp", combiner = "afisher"))$pvalue
  qd.test(x, y, control = list(R = 200, method = "bootp", combiner = "berkjones"))$pvalue

# -- Permutation based tests 
  qd.test(x, y, control = list(R = 2000, method = "permutation", combiner = "minp"))$pvalue
  qd.test(x, y, control = list(R = 2000, method = "permutation", combiner = "afisher"))$pvalue
  qd.test(x, y, control = list(R = 2000, method = "permutation", combiner = "berkjones"))$pvalue



}
}


