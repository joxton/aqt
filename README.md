# aqt
Adaptive quantile based two-sample tests in R 

## The tests 
The function 
```
qd.test(x, y, control = list())
```
performs a two-sample test, comparing sample values `x` with that of `y`.

Different tests are implemented by the function, each based on comparing the samples across a range of quantiles 
specified by `control$tau`. The individual quantile comparisons are pooled to produce the two sample test. The individual quantile comparisons can be rank based, if `control$method = 'mood'` is chosen, or based on quantile differences, if `control$method = 'bootstrap'`. The standard error of the quantile differences are estimated by the non-parametric bootstrap (with `control$R` bootstrap samples).  


## Power difference optimization 

