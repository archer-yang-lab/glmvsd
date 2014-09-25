expnd=function(x){
	n = nrow(x)
	q = ncol(x)
	x <- scale(x)
z = matrix(NA,n,5*q) 
for(k in 1:q){
	i1 = 5*k - 4
	i2 = 5*k - 3
	i3 = 5*k - 2
	i4 = 5*k - 1
	i5 = 5*k
	v = bs(x[,k], df = 5)
	z[,i1] = v[,1]
	z[,i2] = v[,2]
	z[,i3] = v[,3]
	z[,i4] = v[,4]
	z[,i5] = v[,5]
}
return(z)
}


rand_tweedie<- function(mu) {
	Y=rtweedie(1, mu = mu, xi = 1.5, phi = 1)
	Y
}