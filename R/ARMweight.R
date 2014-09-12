ARMweight<-function(x,y,n.sim,cand.mod,n.train){
	
	#cand.mod: m*p matrix, list of candiate model selected. 
	#n.sim: number of replication in data split
	#n.train:number of observation in the training set

	n<-length(y)
	p<-ncol(x)
	
	if(is.matrix(x) == "FASLE") stop("x must be matrix with n rows")
	if(is.vector(y)=="FALSE") stop("y must be a vector")
	
	if(missing(cand.mod)) stop("missing candidate model")
	if(missing(n.train)) stop("missing n.train")


cand.nonzero<-apply(cand.mod,1,sum)
o<-order(cand.nonzero)
model.ordered<-cand.mod[o,]

model<-model.ordered[apply(model.ordered,1,sum)<n.train,]

nonzero<-apply(model,1,sum)

m<-dim(model)[1]

one<-matrix(rep(1,n-n.train),ncol=1)

D1<-matrix(rep(0,n.sim*m),ncol=m)

s1<-matrix(rep(0,n.sim*m),ncol=m)

for (i in 1:n.sim){
	train<-sample(n,n.train,replace=F)
	x.test<-x[-train,]
	y.test<-y[-train]
	x.train<-x[train,]
	y.train<-y[train]

      for (j in 1:m)
	  {
		if (sum(model[j,])==0) 
           {
           LSL<-lm(y.train~1)
           coef.train<-LSL$coef
           pred<-one%*%coef.train
           D1[i,j]<-sum((y.test-pred)^2)
           s1[i,j]<-summary(LSL)$sigma
           }
           if (sum(model[j,])!=0)
           {
           mi.train<-matrix(x.train[,which(model[j,]==1)],ncol=nonzero[j])
           mi.test<-matrix(x.test[,which(model[j,]==1)],ncol=nonzero[j])
           LSL<-lm(y.train~mi.train)
           coef.train<-LSL$coef
           xnew<-cbind(one,mi.test)
           pred<-xnew%*%coef.train
           D1[i,j]<-sum((y.test-pred)^2)
           s1[i,j]<-summary(LSL)$sigma
		}
	  }
}


D2<-NULL
s2<-NULL
cand.mod2<-NULL

for (j in 1:m)    #removing the model which create NA results
    {
     if (all(is.na(D1[,j]))==FALSE&all(is.na(s1[,j]))==FALSE){
     D2<-cbind(D2,D1[,j])
     s2<-cbind(s2,s1[,j])
     cand.mod2<-rbind(cand.mod2,model[j,])
    }
}

E<-matrix(rep(0,n.sim*dim(cand.mod2)[1]),nrow=n.sim)

for(i in 1:n.sim)
  {
   for (j in 1:dim(cand.mod2)[1])
    {
     E[i,j]<-(-n/2)*log(s2[i,j])-(s2[i,j]^(-2))*D2[i,j]/2
    }
  }

for(i in 1:n.sim)
  {
   E.max<-max(E[i,])
   for (j in 1:dim(E)[2])
    {
    	E[i,j]<-E[i,j]-E.max
    }
  }


numerator<-matrix(rep(0,n.sim*dim(cand.mod2)[1]),nrow=n.sim)

 for (i in 1:n.sim){
 for(j in 1:dim(E)[2]){
 if (abs(E[i,j])>150) (numerator[i,j]==0)
 else numerator[i,j]=exp(E[i,j])
 }
 }

denominator<-apply(numerator,1,sum)

w.ARM<-numerator/denominator

weight.ARM<-round(apply(w.ARM,2,mean),5)

outlist<-list(weight.ARM=weight.ARM,cand.model=cand.mod2)
}
