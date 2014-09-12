\name{vsd}
\alias{vsd}
\title{Variable Selection Deviation}
\description{
The package calculate the variable selection deviation (VSD) to measure the uncertainty of the selection in terms of inclusion of predictors in the model.
}
\usage{
vsd(x, y, n.train, n.sim=100, k = 10, best, psi=1, candidate = c("union", "single.path","supplied"),model.supplied, penalty = c("Lasso", "SCAD", "MCP"), weight.function = c("ARM", "ARM.Prior", "AIC", "AIC.Prior", "BIC", "BIC.Prior"))
}
\arguments{
  \item{x}{Matrix of predictors}
  \item{y}{Response}
  \item{n.train}{Size of training set when the weight function is ARM and ARM with prior. If \code{n.train} is not specified, then \code{n.train}=ceiling(n/2)}
  \item{n.sim}{Number of replications when the weight function is ARM and ARM with prior. The default value is \code{n.sim=100}}
  \item{k}{Number of folds. The default value is \code{k=10}}
  \item{best}{The index of base model for calculating the difference. If the best is missing, then the penalty type should be specified. The procedure will fit the base model based on the penalty type. }
  \item{psi}{The size of un-normalized prior weight. The default value \code{psi=1}}
  \item{candidate}{Matrix of candidate model to calculate the weight; candidate model can be supplied by user. If \code{candidate=="union"}, then the candidate model is union of solution path of Lasso, SCAD, and MCP; If \code{candidate=="single.path"}, then the candidate model is the solution path of single penalization procedure. If the \code{candidate="single.path"}, penalty function must be specified. User can set \code{candidate="supplied"} to use own set of candidate mode, and input candidate model as matrix.}
  \item{model.supplied}{Candidate model used by users. When \code{candidate="supplied"}, input the set of candidate model as matrix. }
  \item{penalty}{One of  "Lasso", "SCAD" or "MCP". Default is "Lasso"}
  \item{weight.function}{One of "ARM","ARM with prior","AIC","AIC with prior", "BIC" and "BIC with prior".}
}

\details{
The detail is described in Nan and Yang (2012). 
}
\value{
A "vsd" object is retured. The components are:
 \item{vsd.index}{VSD value}
 \item{selection.size}{If \code{best} is not specified, the selected model based on the \code{penalty} will be provided and the \code{selection.size} is the number of selected variables. Else, \code{selection.size} is the number of selected variables in \code{best}}
 \item{weight}{The weight for each candidate model}
 \item{best}{If \code{best} is not specified, the selected model based on the \code{penalty} will be provided.}
 \item{difference}{The number different variables in the candidate model from the \code{best} model.}
 \item{false.positive}{The upper VSD value of \code{best} model}
 \item{false.negative}{The lower VSD value of \code{best} model}
}
\references{Nan and Yang (2012) "VARIABLE SELECTION DIAGNOSTICS MEASURES FOR HIGH-DIMENSIONAL REGRESSION"}
\author{
Ying Nan and Yuhong Yang
}
\examples{

x<-matrix(rnorm(50*8),ncol=8)
y<-rnorm(50,0,1.5)
fit<-vsd(x,y,n.train=40,n.sim=50,candidate="union",penalty="Lasso",weight.function="ARM")
}

