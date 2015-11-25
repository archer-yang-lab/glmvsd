modelfit <- function(x, y, nfolds, 
	penalty = c("lasso", "scad", "mcp"),
	family=c("gaussian","binomial")) {
  penalty <- match.arg(penalty)
  family <- match.arg(family)
  y <- drop(y)
  y <- as.numeric(y)
  x <- as.matrix(x)
  if(penalty == "lasso"){					
    cvfit <- cv.glmnet(x=x,y=y,nfolds=nfolds,alpha=1,type.measure="mse",maxit=1e6,family=family)
    coefit<-coef(cvfit,s="lambda.min")}
  if(penalty == "mcp"){
    cvfit<-cv.ncvreg(X=x,y=y,nfolds=nfolds,penalty="MCP",family=family,max.iter=1e4)
    coefit<-cvfit$fit$beta[,cvfit$min]}
  if(penalty == "scad"){
    cvfit<-cv.ncvreg(X=x,y=y,nfolds=nfolds,penalty="SCAD",family=family,max.iter=1e4)
    coefit<-cvfit$fit$beta[,cvfit$min]}
  modelfit<- 1-(coefit[-1]==0)
  list(coefit = coefit, modelfit = modelfit)
}