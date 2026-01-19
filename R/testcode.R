if(0){
 
 
library(aqt)

source("~/Desktop/mood/code_clean/sa_power_clean.R")
source("~/Desktop/mood/code_clean/sa_lspline.R")




knots = seq(-4,4,length = 10) 
x = rnorm(1000,sd = .5)
par = c(1:5, 4:1)
h= list(par = par, knots = knots) 
#h = sdescent_cspline (niter = 100, par, x, knots, verbose = FALSE)


#z = sample_cspline(n = 10000, beta = h$par, knots = knots, nint = 1000)
#hist(z) 

library(aqt)

knots = seq(-4,4,length = 10) 
par = c(1:5, 4:1)
h= list(par = par, knots = knots) 

control = list() 
control$method_1 = 1
control$method_2 = 3
control$c1 = 3^2
control$c2 = 3^2
control$tau = c(5:95)/100
control$knots = knots
control$betax = par 
control$betay = par 
control$table_length = 100 

 library(aqt)
x = rnorm(500)
x = rexp(1000)
h = optim_pcd (x, niter = 500, knots = seq(0,max(x)+.5,length = 20)  , beta = NULL)


x = sample_pcd (n = 10000, knots = h$knots, beta = h$beta) 
hist(x, nclass = 100)

system.time({
for(k in 1:1000)
power_pcd(n = 200, B = 200 , control = control)  
})


x = sample_pcd (n = 100, knots = knots, beta = par)
y = sample_pcd (n = 100, knots = knots, beta = par)


system.time({
for(k in 1:1000)
h = qdiff_pcd(x, y, tau = c(1:99)/100, knotsx = knots, betax = par, knotsy = knots, betay = par)
})


w = qntl_qd (x, y, tau = c(1:99)/100,, method = 3, R = 200, combiner = 1)

h = qdiff_pcd(x, y, tau = c(1:99)/100, knotsx = knots, betax = par, knotsy = knots, betay = par)
d = data.frame(qd = w$qdiff, seb = w$se, se_pcd = h$se)


knots = seq(-4,4,length = 10) 
par = c(1:5, 4:1)
h= list(par = par, knots = knots) 

v = NULL 
for(k in 1:1000){
x = sample_pcd (n = 100, knots = knots, beta = h$par)
y = sample_pcd (n = 100, knots = knots, beta = h$par)
hh = qdiff_pcd(x, y, tau = c(1:99)/100, knotsx = knots, betax = h$par, knotsy = knots, betay = h$par)
v = cbind(v, hh$q/hh$se)
}
apply(abs(v), 1, quantile , prob = .95)


# 0.677
system.time({
for(k in 1:100)
m = sample_pcd (n = 200000, knots = knots, beta = h$par)
}) 

system.time({
for(k in 1:100)
z = sample_cspline(n = 10000, beta = h$par, knots = knots, nint = 1000)
})


library(aqt)
x = rexp(100)
y = rexp(100)
h = qntl_qd(x, y , tau = c(1:100)/100, method = 4, R = 1000)

m = qntl_qd(x, y , tau = c(1:100)/100, method = 3, R = 1000)

d = data.frame(se_j = h$se, se_b = m$se)


library(aqt)
n = 100 
x = rnorm(n)
y = rnorm(n)+.3

qntl_c(x = x, y = y)$pvalue
qntl_c(x = x, y = y, control = list(method = "bootstrap", tau = c(5:95)/100))$pvalue 


control = list(method = "bootstrap", combiner = "minp", statistic = "pvalue", R = 200, cregion = NULL)
h  = qntl_c(x = x, y = y, tau = c(1:(n-1))/n, control = control) 

control = list(method = "permutation", combiner = "afisher", statistic = "pvalue", R = 200, cregion = NULL)
h = qntl_c(x = x, y = y, tau = c(1:(n-1))/n, control = control) 

n = 100
system.time({
    v = NULL ; w = NULL
for(k in 1:1000){
x = rnorm(n)
y = rnorm(n)
#x = rexp(n)
#y = rexp(n)
h = qntl_qd(x, y , tau = c(1:(n-1))/n, method = 3, R = 1000)
#z = h$pvalue
#z = h$qdiff/h$se
z = h$se
# h = qntl_qd(x, y , tau = c(1:100)/100, method = 1, R = 200)
# z = h$qdiff/sqrt(h$q2)
#w = c(w, t.test(x,y)$p.value)
v = cbind(v, array(z,dim=c(length(z),1)))
#v = c(v,max(abs(z)),1)
}
}) 
apply(abs(v),1, quantile, prob = .95, na.rm =T)


apply(abs(v),1, mean) #quantile, prob = .95, na.rm =T)


}


