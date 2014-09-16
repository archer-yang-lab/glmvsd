vsd<-function(x,y,
	n_train = ceiling(n/2),
	n_rep=100,
	k=10,
	base_model,
	psi=1,
	candidate=c("union","supplied"),
	candidate_model, 
	weight_fun=c("ARM","ARM.Prior","AIC","AIC.Prior","BIC","BIC.Prior"))
{
	
y <- drop(y)
x <- as.matrix(x)
p<-NCOL(x)
n<-length(y)

candidate=match.arg(candidate)
weight_fun=match.arg(weight_fun)

if (n != NROW(x)) 
    stop("x and y have different number of observations")
if (n_train >= n) 
    stop("Training size must be less than the number of observations")
if(missing(psi)) psi<-1
if(missing(n_rep)) n_rep<-100
if(missing(base_model)) stop("User must provide a base model.")


################Calculate differece d_n_rep#############


if(candidate=="union"){
	
	lassofit<-glmnet(x=x,y=y,alpha=1,maxit=100000000)
	scadfit<-ncvreg(X=x,y=y,family="gaussian",penalty="SCAD",max.iter=10000000)
	mcpfit<-ncvreg(X=x,y=y,family="gaussian",penalty="MCP",max.iter=10000000)


      lasso.path<- as.matrix(lassofit$beta)
      scad.path<- as.matrix(scadfit$beta[-1,])
      mcp.path<- as.matrix(mcpfit$beta[-1,])

      beta.path<-t(cbind(lasso.path,scad.path,mcp.path))
	  ind.path<-(1-(beta.path==0))
	
	##############################################################
	####### unique is to find the intersection of models index
	###############################################################
	  candidate_model=unique(ind.path)
	  rownames(candidate_model) <- NULL
}

if(candidate=="supplied"){
   if(missing(candidate_model)) stop("Users must supply a candidate model.")
   if(is.matrix(candidate_model)!=TRUE) stop("Supplied model must be a matrix.")
   if(ncol(candidate_model)!=ncol(x)) stop("Number of variables in candidate model and x does not match.")		
   
########################   
########################   
##  add 0-1 check
########################   
########################   
}


########################   weight   ###############


if(weight_fun=="ARM") {
	ARM.weight<-ARMweight(x=x,y=y,n_rep=n_rep,candidate_model=candidate_model,n_train=n_train)
	weight<-ARM.weight$weight.ARM
    ending_candidate_model<-ARM.weight$ending_candidate_model
    }

if(weight_fun=="ARM.Prior") {
	ARMP.weight<-ARMPweight(x=x,y=y,n_rep=n_rep,candidate_model=candidate_model,psi=psi,n_train=n_train)
    weight<-ARMP.weight$weight.ARM.Prior
    ending_candidate_model<-ARMP.weight$ending_candidate_model
}

if(weight_fun=="AIC") {
	AIC.weight<-AICweight(x=x,y=y,candidate_model=candidate_model)
	weight<-AIC.weight$weight.AIC
    ending_candidate_model<-AIC.weight$ending_candidate_model
    }

if(weight_fun=="AIC.Prior") {
	AICP.weight<-AICPweight(x=x,y=y,candidate_model=candidate_model,psi=psi)
	weight<-AICP.weight$weight.AIC.Prior
    ending_candidate_model<-AICP.weight$ending_candidate_model
    }

if(weight_fun=="BIC") {
	BIC.weight<-BICweight(x=x,y=y,candidate_model=candidate_model)
	weight<-BIC.weight$weight.BIC
    ending_candidate_model<-BIC.weight$ending_candidate_model
	    }

if(weight_fun=="BIC.Prior") {
	BICP.weight<-BICPweight(x=x,y=y,candidate_model=candidate_model,psi=psi)
	weight<-BICP.weight$weight.BIC.Prior
    ending_candidate_model<-BICP.weight$ending_candidate_model
    }

	
DIFF<-symdiff(ending_candidate_model,base_model)

TMP_matrix = sweep(candidate_model,MARGIN=2,base_model,"-")
DIFF_minus<-rowSums(TMP_matrix==-1)
DIFF_plus<-rowSums(TMP_matrix==1)
	

VSD<-weight%*%DIFF                 # vsd value
VSD_minus<-weight%*%DIFF_minus    #false positive
VSD_plus<-weight%*%DIFF_plus    #false negative
	

object<-list(VSD=VSD,VSD_minus=VSD_minus,VSD_plus=VSD_plus,weight=weight,difference=DIFF)
class(object)<-"vsd"
object
}

