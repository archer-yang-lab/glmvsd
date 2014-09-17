ARMPweight<-function(x,y,n_rep=100,psi=1,candidate_model,n_train=ceiling(n/2)){
	y <- drop(y)
	x <- as.matrix(x)
	p<-NCOL(x)
	n<-length(y)

	if (n != NROW(x)) 
	    stop("x and y have different number of observations")
	if (n_train >= n) 
	    stop("Training size must be less than the number of observations")

	   if(missing(candidate_model)) stop("Users must supply a candidate model.")
	   if(is.matrix(candidate_model)!=TRUE) stop("Supplied model must be a matrix.")
	   if(NCOL(candidate_model)!=NCOL(x)) stop("Number of variables in candidate model and x does not match.")


cand.nonzero<-rowSums(candidate_model)
o<-order(cand.nonzero)
model.ordered<-candidate_model[o,]
model<-model.ordered[rowSums(model.ordered)<n_train, ]
nonzero<-rowSums(model)

m <- NROW(model)
d1<-matrix(NA, n_rep, m)
s1<-matrix(NA, n_rep, m)

for (i in 1:n_rep){
	tindex<-sample(n,n_train,replace=F)
	if (any(model[1,]==1))
    {
			for (j in 1:m)
			{
	           LSL<-lm(y[tindex]~x[tindex, model[j,]==1])
	           d1[i,j]<-sum((y[-tindex]-cbind(1,x[-tindex, model[j,]==1])%*%LSL$coef)^2)
	           s1[i,j]<-summary(LSL)$sigma
			}
		
	}
	else{
			d1[i,1]<-sum((y[-tindex]-mean(y[tindex]))^2)
	        s1[i,1]<-sd(y[tindex])
		    for (j in 2:m)
			{
	           LSL<-lm(y[tindex]~x[tindex, model[j,]==1])
	           d1[i,j]<-sum((y[-tindex]-cbind(1,x[-tindex, model[j,]==1])%*%LSL$coef)^2)
	           s1[i,j]<-summary(LSL)$sigma
			}
        }
}


D2<-NULL
s2<-NULL
candidate_model2<-NULL

for (j in 1:m)  #removing the model which create NA results
    {
     if (all(is.na(d1[,j]))==FALSE&all(is.na(s1[,j]))==FALSE){
     D2<-cbind(D2,d1[,j])
     s2<-cbind(s2,s1[,j])
     candidate_model2<-rbind(candidate_model2,model[j,])
    }
}

pstar<-apply(candidate_model2,1,sum)   #non-zero variables in candidate model after removing model with NA results

prior.arm<-rep(0,dim(D2)[2])

for (i in 1:dim(D2)[2]){
 if(pstar[i]==0) {prior.arm[i]<-2*log(pstar[i]+2)/choose(p,pstar[i])}
 if(pstar[i]!=0){
 prior.arm[i]<-pstar[i]*log(exp(1)*p/pstar[i])+2*log(pstar[i]+2)}
 }

E.prior<-matrix(rep(0,n_rep*dim(D2)[2]),nrow=n_rep)

for(i in 1:n_rep)
  {
   for (j in 1:dim(D2)[2])
    {
     E.prior[i,j]<-(-n/2)*log(s2[i,j])-(s2[i,j]^(-2))*D2[i,j]/2-psi*prior.arm[j]
    }
  }

for(i in 1:n_rep)
  {
   E.prior.max<-max(E.prior[i,])
   for (j in 1:dim(E.prior)[2])
    {
    	E.prior[i,j]<-E.prior[i,j]-E.prior.max
    }
  }

numerator.prior<-matrix(rep(0,n_rep*dim(D2)[2]),nrow=n_rep)

 for (i in 1:n_rep){
 for(j in 1:dim(E.prior)[2]){
 if (abs(E.prior[i,j])>700) (numerator.prior[i,j]=0)
 else numerator.prior[i,j]=exp(E.prior[i,j])
 }
 }

denominator.prior<-apply(numerator.prior,1,sum)

w.prior<-numerator.prior/denominator.prior

weight.arm.prior<-round(apply(w.prior,2,mean),5)

outlist<-list(weight.ARM.Prior=weight.arm.prior,ending_candidate_model=candidate_model2)

}





