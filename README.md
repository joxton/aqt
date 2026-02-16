# aqt
Adaptive quantile based two-sample tests in R 

## The tests 
The function 
```
qd.test(x, y, control = list())
```
performs a two-sample test, comparing sample values `x` with that of `y`.

The tests implemented by `qd.test` are each based on comparing the `x` and `y` across the quantiles 
specified by `control$tau`, and then pooling the quantile comparisons to produce the two sample test. The individual quantile comparisons can be rank-based (setting `control$method = 'mood'`),  or based on quantile differences, if `control$method = 'bootstrap'`. The standard error of the quantile differences are estimated by the non-parametric bootstrap (with `control$R` bootstrap samples).  

Setting `control$combiner` determines how the individual quantile comparisons are pooled to form the two sample test, an can be one of `minp`, `berkjones` and `afisher`. 


Running 
```
qd.test(x, y)
```
i.e. without specifying `control`, 

## Power difference optimization 

The function 
```
pwr_optim(control = list())
```
maximizes the power difference between two-sample tests `control$test1` and `control$test2`.  

For example, 
```
test1 = function(x, y) return(t.test(x,y)$p.value)
test2 = function(x, y) return(wilcox.test(x,y)$p.value)

res = pwr_optim(control = list(test1 = test1, test2 = test2))
```
attempts to find a density function for `x` and one for `y` that maximizes the power of `test1` over `test2`. The optimization is over desities that are piecewice constant between `control$knotsx` and `control$knotsy`, and by default, both are initialized to `seq(-4,4,length = 10)`. The density values are the parameters `control$betax`and `control$betay`, which are optimized over. The optimization uses a stochastic cyclic ascent approach, sampling candidate values of the components of the densities and accepting those leading to an increase in estimated power. When the opimization plateaus, the knot-sequences are doubled, and the optimization continued. Once knot-doubling no longer increases the power difference, the optimization terminates. Multiple runs should be used for reliable results. 

Constraints on the optimization are set by `control$constraint` and `control$location_shift`. By default, `control$location_shift = TRUE` and  `control$constraint = 'unimodal'`, and the power difference is optimized is under a location-shift model with the densities constrained to be unimodal. Setting 
```
res = pwr_optim(control = list(test1 = test1, test2 = test2,
                location_shift = FALSE, constraint = "dominance"))
```
maximizes the power difference under the constraint that one of the distributions stochastically dominates the other. 
 
