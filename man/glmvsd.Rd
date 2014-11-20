\name{glmvsd}
\alias{glmvsd}
\title{Variable Selection Deviation (VSD)}
\description{
The package calculate the variable selection deviation (VSD) to measure the uncertainty of the selection in terms of inclusion of predictors in the model.
}
\usage{
glmvsd(x, y, n_train = ceiling(n/2), 
no_rep = 100, model_check, 
psi = 1, family = c("gaussian", "binomial"), 
method = c("union", "customize"), 
candidate_models, 
weight_function = c("ARM", "BIC"), prior = TRUE)
}
\arguments{
  \item{x}{Matrix of predictors.}

  \item{y}{Response variable.}

  \item{n_train}{Size of training set when the weight function is ARM or ARM with prior. The default value is \code{n_train}=ceiling(n/2).}

  \item{no_rep}{Number of replications when the weight function is ARM and ARM with prior. The default value is \code{no_rep=100}.}

  \item{model_check}{The index of the model to be assessed by calculating the VSD measures.}

  \item{psi}{A positive number to control the improvement of the prior weight. The default value is 1.}
  \item{family}{Choose the family for GLM models. So far only \code{gaussian}, \code{binomial} and \code{tweedie} are implemented. The default is \code{gaussian}.}
  \item{method}{User chooses one of the \code{union} and \code{customize}. If \code{method=="union"}, then the program automatically provides the candidate models as a union of solution paths of Lasso, SCAD, and MCP; If \code{method="customize"}, the user must provide their own set of candidate models in the input argument \code{candidate_models} as a matrix, each row of which is a 0/1 index vector representing whether each variable is included/excluded in the model.}

  \item{candidate_models}{Only available when \code{method="customize"}. It is a matrix of candidate models, each row of which is a 0/1 index vector representing whether each variable is included/excluded in the model.}

  \item{weight_function}{Options for computing weights for VSD measure. User chooses one of the \code{ARM}, \code{BIC}. The default is \code{ARM}.}
  \item{prior}{Whether use prior in the weight function. The default is \code{TRUE}.}
}

\details{
See Reference section. 
}

\value{
A "glmvsd" object is retured. The components are:
 \item{VSD}{Variable selection deviation (VSD) value.}
 \item{VSD_minus}{The lower VSD value of \code{model_check}, representing the number of predictors in the model (\code{model_check}) not quite justified at the present sample size.}
 \item{VSD_plus}{The upper VSD value of \code{model_check} model, representing the number of predictors missed by the model (\code{model_check}).}
 \item{weight}{The weight for each candidate model.}
 \item{DIFF}{Counting the variable differences between candidate models and \code{model_check}.}
\item{candidate_models_cleaned}{Cleaned candidate models: the duplicated candidate models are cleaned; When computing VSD weights using BIC, the models with more than n-2 variables are removed (n is the number of observaitons); When computing VSD weights using ARM, the models with more than n_train-2 variables  are removed (n_train is the number of training observations).}
}

\references{
Nan, Y. and Yang, Y. (2013), "Variable Selection Diagnostics Measures for High-dimensional Regression," \emph{Journal of Computational and Graphical Statistics}, 23:3, 636-656.\cr
\url{http://dx.doi.org/10.1080/10618600.2013.829780}\cr
BugReport: \url{https://github.com/emeryyi/glmvsd}\cr
}

\examples{
# REGRESSION CASE

# generate simulation data
n <- 50
p <- 8
beta <- c(3,1.5,0,0,2,0,0,0)
sigma <- matrix(0,p,p)
for(i in 1:p){
   for(j in 1:p) sigma[i,j] <- 0.5^abs(i-j)
}
x <- mvrnorm(n, rep(0,p), sigma)
e <- rnorm(n)
y <- x \%*\% beta + e

# user provide a model to be checked
model_check <- c(0,1,1,1,0,0,0,1)

# compute VSD for model_check using ARM with prior
v_ARM <- glmvsd(x, y, n_train = ceiling(n/2), 
no_rep=50, model_check = model_check, psi=1, 
family = "gaussian", method = "union", 
weight_function = "ARM", prior = TRUE)

# compute VSD for model_check using BIC
v_BIC <- glmvsd(x, y, 
model_check = model_check,
family = "gaussian", method = "union", 
weight_function = "BIC", prior = TRUE)

# user supplied candidate models
candidate_models = rbind(c(0,0,0,0,0,0,0,1), 
c(0,1,0,0,0,0,0,1), c(0,1,1,1,0,0,0,1), 
c(0,1,1,0,0,0,0,1), c(1,1,0,1,1,0,0,0), 
c(1,1,0,0,1,0,0,0))

v1_BIC <- glmvsd(x, y, 
model_check = model_check, psi=1, 
family = "gaussian",
method = "customize", 
candidate_models = candidate_models, 
weight_function = "BIC", prior = TRUE)

# CLASSIFICATION CASE

# generate simulation data
n = 300 
p = 8
b <- c(1,1,1,-3*sqrt(2)/2) 
x=matrix(rnorm(n*p, mean=0, sd=1), n, p) 
feta=x[, 1:4]\%*\%b 
fprob=exp(feta)/(1+exp(feta))
y=rbinom(n, 1, fprob)

# user provide a model to be checked
model_check <- c(0,1,1,1,0,0,0,1)

# compute VSD for model_check using BIC with prior
b_BIC <- glmvsd(x, y, n_train = ceiling(n/2),
family = "binomial",
no_rep=50, model_check = model_check, psi=1, 
method = "union", weight_function = "BIC", 
prior = TRUE)

candidate_models = 
rbind(c(0,0,0,0,0,0,0,1), 
c(0,1,0,0,0,0,0,1), 
c(1,1,1,1,0,0,0,0), 
c(0,1,1,0,0,0,0,1), 
c(1,1,0,1,1,0,0,0), 
c(1,1,0,0,1,0,0,0), 
c(0,0,0,0,0,0,0,0),
c(1,1,1,1,1,0,0,0))

# compute VSD for model_check using BIC
# user supplied candidate models
b_BIC1 <- glmvsd(x, y, 
family = "binomial",
model_check = model_check, psi=1, 
method = "customize", 
candidate_models = candidate_models, 
weight_function = "BIC")
}

