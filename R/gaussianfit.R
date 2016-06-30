gaussianfit <- function(x, y) {
  m1 <- glmnet(x = x, y = y, family = "gaussian", alpha = 1, maxit = 1e+06)
  m2 <- ncvreg(X = x, y = y, family = "gaussian", penalty = "SCAD",
                   warn = FALSE, max.iter = 1e+04)
  m3 <- ncvreg(X = x, y = y, family = "gaussian", penalty = "MCP",
                   warn = FALSE, max.iter = 1e+04)
  #############################			   
  ## adaptive lasso routine  ##			   
  ## The adaptive lasso needs a first stage that is consistent. 
  ## Zou (2006) recommends OLS or ridge
  thelasso.cv<-cv.glmnet(x,y,family = "gaussian",alpha=1) ## first stage ridge
  ## Second stage weights from the coefficients of the first stage
  bhat<-as.matrix(coef(thelasso.cv,s="lambda.1se"))[-1,1] ## coef() is a sparseMatrix
  if(all(bhat==0)){
    ## if bhat is all zero then assign very close to zero weight to all.
    ## Amounts to penalizing all of the second stage to zero.
    bhat<-rep(.Machine$double.eps*2,length(bhat))
  }
  adpen<-(1/pmax(abs(bhat),.Machine$double.eps)) ## the adaptive lasso weight
  ## Second stage lasso (the adaptive lasso)
  m4 <- glmnet(x,y,family = "gaussian",alpha=1,exclude=which(bhat==0),penalty.factor=adpen)
  #############################					 				   
  m1.path <- as.matrix(m1$beta)
  m2.path <- as.matrix(m2$beta[-1, ])
  m3.path <- as.matrix(m3$beta[-1, ])
  m4.path <- as.matrix(m4$beta)
  
  beta.path <- t(cbind(m1.path, m2.path, m3.path, m4.path))
  candidate_models <- (1 - (beta.path == 0))
  candidate_models
}