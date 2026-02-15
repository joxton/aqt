# aqt
Adaptive quantile based two-sample tests in R 

The function 
```
qd.test(x, y, control = list())
```
performs adaptive twp-sample quantile-based tests, comparing sample values `x` with that of `y`. All tests compare sample quantiles across the quantiles specified by `control$tau`, and produce a test combining across the individual quantile comparisons. 
