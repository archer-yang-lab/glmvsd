BICPweight<-function(x,y,cand.mod,psi){
	
	n<-length(y)
	p<-ncol(x)
	
	if(is.matrix(x) == "FASLE") stop("x must be matrix with n rows")
	if(is.vector(y)=="FALSE") stop("y must be a vector")
	
	if(missing(cand.mod)) stop("missing candidate model")
	if(!missing(psi)) psi<-psi
	  else psi<-1
	

	cand.nonzero<-apply(cand.mod,1,sum)
	o<-order(cand.nonzero)
	model<-cand.mod[o,]
    nonzero<-apply(model,1,sum)
    m<-dim(model)[1]
	
	prior<-rep(0,m)

	for (i in 1:m){
	 if(nonzero[i]==0) {prior[i]<-2*log(nonzero[i]+2)/choose(p,nonzero[i])}
	 if(nonzero[i]!=0){
	 prior[i]<-nonzero[i]*log(exp(1)*p/nonzero[i])+2*log(nonzero[i]+2)}
	 }	
	

BIC.prior<-rep(0,m)

for (i in 1:m){
if(nonzero[i]==0){
LSL<-lm(y~1) 
rss<-sum(summary(LSL)$res^2)
BIC.prior[i]<-n*log(rss/n)+nonzero[i]*log(n)+psi*prior[i]
}
if(nonzero[i]!=0){
x.new<-matrix(x[,which(model[i,]==1)],ncol=nonzero[i])
LSL<-lm(y~x.new)
rss<-sum(summary(LSL)$res^2)
BIC.prior[i]<-n*log(rss/n)+nonzero[i]*log(n)+psi*prior[i]

}
}

 m2<-NULL
 BIC.prior2<-NULL

 for(i in 1:m){
 if (is.na(BIC[i])==FALSE&abs(BIC[i])!=Inf){
 m2<-rbind(m2,model[i,])
 BIC.prior2<-c(BIC.prior2,BIC.prior[i])
 }
 }

BIC.new.prior<-BIC.prior2-min(BIC.prior2)

weight.BIC.prior<-round(exp(-BIC.new.prior/2)/sum(exp(-BIC.new.prior/2)),5)

outlist<-list(weight.BIC.Prior=weight.BIC.prior,ending_candidate_model=m2)
}