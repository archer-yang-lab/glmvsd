BICweight<-function(x,y,cand.mod){
	
	n<-length(y)
	p<-ncol(x)
	
	
	if(is.matrix(x) == "FASLE") stop("x must be matrix with n rows")
	if(is.vector(y)=="FALSE") stop("y must be a vector")
	
	if(missing(cand.mod)) stop("missing candidate model")


	cand.nonzero<-apply(cand.mod,1,sum)
	o<-order(cand.nonzero)
	model<-cand.mod[o,]
    nonzero<-apply(model,1,sum)
    m<-dim(model)[1]


BIC<-rep(0,m)

for (i in 1:m){
if(nonzero[i]==0){
LSL<-lm(y~1) 
rss<-sum(summary(LSL)$res^2)
BIC[i]<-n*log(rss/n)+nonzero[i]*log(n)
}
if(nonzero[i]!=0){
x.new<-matrix(x[,which(model[i,]==1)],ncol=nonzero[i])
LSL<-lm(y~x.new)
rss<-sum(summary(LSL)$res^2)
BIC[i]<-n*log(rss/n)+nonzero[i]*log(n)
}
}


 BIC2<-NULL
 m2<-NULL

 for(i in 1:m){
 if (is.na(BIC[i])==FALSE&abs(BIC[i])!=Inf){
 BIC2<-c(BIC2,BIC[i])
 m2<-rbind(m2,model[i,])
 }
 }

BIC.new<-BIC2-min(BIC2)

weight.BIC<-round(exp(-BIC.new/2)/sum(exp(-BIC.new/2)),5)

outlist<-list(weight.BIC=weight.BIC,ending_candidate_model=m2)
}