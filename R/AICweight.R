AICweight<-function(x,y,candidate_model){
	
	n<-length(y)
	p<-NCOL(x)
	
	
	if(is.matrix(x) == "FASLE") stop("x must be matrix with n rows")
	if(is.vector(y)=="FALSE") stop("y must be a vector")
	
	if(missing(candidate_model)) stop("missing candidate model")


	cand.nonzero<-apply(candidate_model,1,sum)
	o<-order(cand.nonzero)
	model<-candidate_model[o,]
    nonzero<-apply(model,1,sum)
    m<-dim(model)[1]

	

AIC<-rep(0,m)

for (i in 1:m){
if(nonzero[i]==0){
LSL<-lm(y~1) 
rss<-sum(summary(LSL)$res^2)
AIC[i]<-n*log(rss/n)+nonzero[i]*2
}
if(nonzero[i]!=0){
x.new<-matrix(x[,which(model[i,]==1)],ncol=nonzero[i])
LSL<-lm(y~x.new)
rss<-sum(summary(LSL)$res^2)
AIC[i]<-n*log(rss/n)+nonzero[i]*2
}
}


 AIC2<-NULL
 m2<-NULL

 for(i in 1:m){
 if (is.na(AIC[i])==FALSE&abs(AIC[i])!=Inf){
 AIC2<-c(AIC2,AIC[i])
 m2<-rbind(m2,model[i,])
 }
 }

AIC.new<-AIC2-min(AIC2)

weight.AIC<-round(exp(-AIC.new/2)/sum(exp(-AIC.new/2)),5)

outlist<-list(weight.AIC=weight.AIC,ending_candidate_model=m2)
}