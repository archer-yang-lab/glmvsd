vsd<-function(x,y,n.train,n.sim=100,k=10,best,psi=1,candidate=c("union","single.path","supplied"),model.supplied,penalty=c("Lasso","SCAD","MCP"),weight.function=c("ARM","ARM.Prior","AIC","AIC.Prior","BIC","BIC.Prior"))

{
# x = n by p covariate matrix
# y = size n vector
# n.train = the size of training set in MARM
# n.s = number of regularization in penalization procedure
# n.sim = number of data spliting in MARM
# k = number of fold in cross validation
# candidate = way to get candidate model; 
#           "union" is union of solution path of Lasso, SCAD, and MCP; 
#           "single.path" is the solution path of single penalization procedure. Penalty                
#            function must be specified. 
# best = base model used in calculating indice.
# psi = the size of penalty in calculting the non-uniformed AIC or BIC vsd

penalty=match.arg(penalty)

n<-length(y)
np<-dim(x)
p<-np[2]

if(missing(n.train)){n.train<-ceiling(n/2)}
if(is.matrix(x) == "FASLE") stop("x must be matrix with n rows")
if(is.vector(y)=="FALSE") stop("y must be a vector")
if(missing(weight.function)) stop("missing weight function!")
   else weight.function=match.arg(weight.function)

if(missing(candidate)) stop("missing the method of getting candidate model")
   else candidate=match.arg(candidate)

if(!missing(psi)) psi<-psi
   else psi<-1

if(!missing(n.sim)) n.sim<-n.sim
   else n.sim<-100



lassofit<-glmnet(x=x,y=y,alpha=1,maxit=100000000)
scadfit<-ncvreg(X=x,y=y,family="gaussian",penalty="SCAD",max.iter=10000000)
mcpfit<-ncvreg(X=x,y=y,family="gaussian",penalty="MCP",max.iter=10000000)


################Calculate differece d_n.sim#############



if(candidate=="union"){
      lasso.path<-t(as.matrix(lassofit$beta))
      scad.path<-t(as.matrix(scadfit$beta[-1,]))
      mcp.path<-t(as.matrix(mcpfit$beta[-1,]))

      beta.path<-rbind(lasso.path,scad.path,mcp.path)
	  ind.path<-matrix((1-as.numeric(beta.path==0)),ncol=p)
	
	##############################################################
	####### unique is to find the intersection of models index
	###############################################################
	  cand.mod=unique(ind.path)
     }

if(candidate=="single.path") {
      if(penalty=="Lasso"){
		ind.path<-matrix(as.numeric(t(as.matrix(lassofit$beta))!=0),ncol=p)
		cand.mod=unique(ind.path)}
      if(penalty=="SCAD"){
		ind.path<-matrix(as.numeric(t(as.matrix(scadfit$beta))!=0),ncol=p)
		cand.mod=unique(ind.path)}
      if(penalty=="MCP"){
		ind.path<-matrix(as.numeric(t(as.matrix(mcpfit$beta))!=0),ncol=p)
		cand.mod=unique(ind.path)}
      }

if(candidate=="supplied"){
   if(missing(model.supplied)) stop("Supplied model is missing. ")
   if(is.matrix(model.supplied)!=TRUE) stop("Supplied model must be a matrix. ")
   if(ncol(model.supplied)!=ncol(x)) stop("Number of variables in candidate model and x does not match. ")		
   cand.mod<-model.supplied
}


########################   weight   ###############


if(weight.function=="ARM") {
	ARM.weight<-ARMweight(x=x,y=y,n.sim=n.sim,cand.mod=cand.mod,n.train=n.train)
	weight<-ARM.weight$weight.ARM
    model<-ARM.weight$cand.model
    }

if(weight.function=="ARM.Prior") {
	ARMP.weight<-ARMPweight(x=x,y=y,n.sim=n.sim,cand.mod=cand.mod,psi=psi,n.train=n.train)
    weight<-ARMP.weight$weight.ARM.Prior
    model<-ARMP.weight$cand.model
}

if(weight.function=="AIC") {
	AIC.weight<-AICweight(x=x,y=y,cand.mod=cand.mod)
	weight<-AIC.weight$weight.AIC
    model<-AIC.weight$cand.model
    }

if(weight.function=="AIC.Prior") {
	AICP.weight<-AICPweight(x=x,y=y,cand.mod=cand.mod,psi=psi)
	weight<-AICP.weight$weight.AIC.Prior
    model<-AICP.weight$cand.model
    }

if(weight.function=="BIC") {
	BIC.weight<-BICweight(x=x,y=y,cand.mod=cand.mod)
	weight<-BIC.weight$weight.BIC
    model<-BIC.weight$cand.model
	    }

if(weight.function=="BIC.Prior") {
	BICP.weight<-BICPweight(x=x,y=y,cand.mod=cand.mod,psi=psi)
	weight<-BICP.weight$weight.BIC.Prior
    model<-BICP.weight$cand.model
    }

####################################################

if(!missing(best) & missing(penalty)) stop("either base model or penalty need to defined. ")

if (k < 3) 
    stop("k must be bigger than 3; k=10 recommended")

if(missing(best) & penalty=="Lasso"){
   	cv.lasso<-cv.glmnet(x=x,y=y,nfolds=k,alpha=1,type.measure="mse",maxit=100000000)
	coef.lasso<-coef(cv.lasso,s="lambda.min")[-1]
	best<-(1-as.numeric(coef.lasso==0))
}


if(missing(best) & penalty=="SCAD"){
	cv.scad<-cv.ncvreg(X=x,y=y,nfolds=k,penalty="SCAD",family="gaussian",max.iter=10000000)
	coef.scad<-scadfit$beta[,cv.scad$min]
	best<-matrix((1-as.numeric(coef.scad==0)),ncol=p)
}


if(missing(best) & penalty=="MCP"){
	cv.mcp<-cv.ncvreg(X=x,y=y,nfolds=k,penalty="MCP",family="gaussian",max.iter=10000000)
	coef.mcp<-mcpfit$beta[,cv.mcp$min]
	best<-matrix((1-as.numeric(coef.mcp==0)),ncol=p)
}
	
diff<-symdiff(model,best)
missing1<-NULL
missing2<-NULL

for(i in 1:dim(model)[1])
    {
	missing1[i]<-sum(as.numeric(best-model[i,]==-1))
	missing2[i]<-sum(as.numeric(best-model[i,]==1))
	}
	
size<-sum(best)

index<-weight%*%diff                 # vsd value

index.missing1<-weight%*%missing1    #false positive
index.missing2<-weight%*%missing2    #false negative
	
	

object<-list(vsd.index=index,selection.size=size,weight=weight,best=best,difference=diff,false.positive=index.missing1,false.negative=index.missing2)
class(object)<-"vsd"
object

}

