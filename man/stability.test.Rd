\name{stability.test}
\alias{stability.test}
\title{
Instability tests
}
\description{
This function calculate the sequential, parametric bootstrap and perturbation instability measures for linear regression with Lasso, SCAD and MCP penalty. 
}
\usage{
stability.test(x, y, 
method = c("seq", "bs", "perturb"), 
penalty = c("LASSO", "SCAD", "MCP"), 
nrep = 50, remove = 0.2, tau = 0.5, nfolds = 5,
family=c("gaussian","binomial"))
}

\arguments{
  \item{x}{Matrix of predictors. }
  \item{y}{Response variable. }
  \item{method}{Type of instability measures. \code{seq} = sequential instability, \code{bs} = parametric bootstrap instability, and \code{perturb} = perturbation instability. }
  \item{penalty}{Penalty function. }
  \item{nrep}{Number of repetition for calculating instability, default is 50. }
  \item{remove}{The portion of observation to be removed when the sequential instability is calculated, default is 0.2.  }
  \item{tau}{The size of perturbation when perturbation instability is calculated. The range of \code{tau} is (0,1), default is 0.5 }
  \item{nfolds}{number of folds - default is 5. }
  \item{family}{Choose the family for the instability test. So far only \code{gaussian}, \code{binomial} and \code{tweedie} are implemented. The default is \code{gaussian}.}
}
\details{
See Reference section. 
}
\value{
Return the instability index according to the type of instability measures. 
}
\references{
Nan, Y. and Yang, Y. (2013), "Variable Selection Diagnostics Measures for High-dimensional Regression," \emph{Journal of Computational and Graphical Statistics}, 23:3, 636-656.
BugReport: \url{https://github.com/emeryyi/glmvsd}\cr
}


\examples{

# generate simulation data
n <- 50
p <- 8
beta<-c(2.5,1.5,0.5,rep(0,5))
sigma<-matrix(0,p,p)
for(i in 1:p){
   for(j in 1:p) sigma[i,j] <- 0.5^abs(i-j)
}
x <- mvrnorm(n, rep(0,p), sigma)
e <- rnorm(n)
y <- x \%*\% beta + e

ins_seq <- stability.test(x, y, method = "seq", 
penalty = "SCAD", nrep = 20, 
remove = 0.1, tau = 0.2, nfolds = 5)
}

