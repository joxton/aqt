\name{pwr_optim }
\alias{pwr_optim }
\title{Power difference optimizer}
\description{Estimate the maximal power difference between two-sample tests .}
\usage{

pwr_optim (control = list())  

}
\arguments{
\item{control}{A list of arguments for controlling the optimization (see below). } 
}

\details{
\code{pwr_optim} maximizes the difference in power between the two tests in \code{control},  \code{test1} and \code{test2}. Both \code{test1} and \code{test2} must be functions of two arguments (the vectors of sample values from the two groups), 
and return the p-value of the test. For example \code{test1 = function(x, y) return(t.test(x, y)$p.value)}.  The optimization maximizes the power of \code{test1} minus the power of \code{test2} with respect to piecewise constant densities. 
The density in one group is propotional to \code{control$betax} over the knots \code{control$knotsx},  the other to \code{control$betay} over the knots \code{control$knotsy}, which are optimized over. The optimization uses knot doubling, 
when it plateaus for a give set of knots, and it terminated once either the number of doublings exceeds \code{control$max_doubling} or doesn't improve after \code{control$max_doubling_stagnant} doublings. 


Different constraints can be placed on the densities being optimized over. If \code{control$location_shift} is \code{TRUE} and \code{control$constraint} is \code{"unimodal"} the optimization is over location-shift variants of unimodal densities. 
 If \code{control$location_shift} if \code{FALSE} and \code{control$constraint = "dominance"} the optimization imposes stochastic dominance of one distribution over the other. 

The optimization uses sampling based approach and multiple runs should be used. The optimization takes time.    

The full list of control elements and their default settings
\preformatted{
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
}


}

\value{ Returns a \code{list} with the control arguments used along with the optimizing densities and the estimated power difference. }  

\references{

}


\seealso{}

\examples{

\dontrun{

# Maximize the power difference of the t.test over median test 
median.test = function(x, y){ z = c(x, y);  zm = as.numeric(z > median(z)); 
		gr = as.numeric(1:length(z)<=length(x));  return(chisq.test( table(zm, gr))$p.value)} 

test1 = function(x, y) return( t.test(x, y)$p.value)
test2 = median.test 

r = pwr_optim(control = list(test1 = test1, test2 = test2) ) 

xlim <- c(min(c(r$knotsx, r$knotsy)), max(c(r$knotsx, r$knotsy)))
plot(r$knotsx[-1], r$betax, type = "l", xlim = xlim, xlab = "", ylab = "density")
points(r$knotsy[-1], r$betay, type = "l", col = "blue")




# Maximize the power difference of a rank based adaptive quantile test and the t-test,  

test1 = function(x, y) return( qd.test(x, y, control = list(method = "mood", combiner = "afisher"))$pvalue)
test2 = function(x, y) return( t.test(x, y)$p.value)

r = pwr_optim(control = list(test1 = test1, test2 = test2, constraint = "dominance") ) 

xlim <- c(min(c(r$knotsx, r$knotsy)), max(c(r$knotsx, r$knotsy)))
plot(r$knotsx[-1], r$betax, type = "l", xlim = xlim, xlab = "", ylab = "density")
points(r$knotsy[-1], r$betay, type = "l", col = "blue")


}
}


