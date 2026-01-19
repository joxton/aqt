if(0){

library(aqt)

knots = seq(-4,4,length = 20) 
par = c(1:10, 9:1)
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
n = 200 
B = 500



for(k in 1:10){
	for(j in 1:length(control$betax)){ 
		p0 = power_pcd(n, B, control = control)
		b0 = control$betax[j] 
		#r = get_range(j, control$betax)
		
		r = get_range_stdom(j, control$betax, control$betay) 
		control$betax[j] = b0 + runif(1,min = r$min_j, max = r$max_j)
		p1 = power_pcd(n, B, control = control)
		if(p1 < p0){
		 	control$betax[j] = b0
		}
	}
	for(j in 1:length(control$betay)){ 
		p0 = power_pcd(n, B, control = control)
		b0 = control$betay[j] 
		# r = get_range(j, control$betay)
		r = get_range_stdom(j+length(control$betax), control$betax, control$betay) 
		control$betay[j] = b0 +  runif(1,min = r$min_j, max = r$max_j)
		p1 = power_pcd(n, B, control = control)
		if(p1 < p0){
		 	control$betay[j] = b0
		}
	}
	
	print(power_pcd(n, B, control = control))
}



x = sample_pcd (n = 10000, knots = control$knots, beta = control$betax)
y = sample_pcd (n = 10000, knots = control$knots, beta = control$betay)
hist(x)
hist(y)
plot(ecdf(x)) 
fy = ecdf(y)(seq(-4,4,length= 100)) 
points(seq(-4,4,length= 100), fy ,type = "l", col = "blue")

}

