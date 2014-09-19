\name{glmvsd}
\alias{glmvsd}
\title{Variable Selection Deviation}
\description{
The package calculate the variable selection deviation (VSD) to measure the uncertainty of the selection in terms of inclusion of predictors in the model.
}
\usage{
glmvsd(x, y, n_train, n_rep=100, k = 10, model_check, psi=1, candidate = c("union","supplied"),candidate_models, weight_fun = c("ARM", "ARM.Prior","BIC", "BIC.Prior"))
}
\arguments{
  \item{x}{Matrix of predictors}
  \item{y}{Response variable}
  \item{n_train}{Size of training set when the weight function is ARM or ARM with prior. The default value is \code{n_train}=ceiling(n/2).}
  \item{n_rep}{Number of replications when the weight function is ARM and ARM with prior. The default value is \code{n_rep=100}.}
  \item{model_check}{The index of base model for calculating the difference. If the \code{model_check} is missing, then the penalty type should be specified. The procedure will fit the base model based on the penalty type.}
  \item{psi}{The size of un-normalized prior weight. The default value is 1}
  \item{candidate}{Matrix of candidate model to calculate the weight. If \code{candidate=="union"}, then the candidate models are a union of solution path of Lasso, SCAD, and MCP; If \code{candidate="supplied"}, users can use their own set of candidate models, and input candidate model as a matrix, each row of which is a 0/1 index vector representing whether the variable included/excluded in the model.}
  \item{candidate_models}{Candidate model used by users. When \code{candidate="supplied"}, input the set of candidate model as matrix. }
  \item{weight_fun}{One of the "ARM", "ARM with prior", "BIC" and "BIC with prior".}
}

\details{
The detail is described in Nan and Yang (2012). 
}
\value{
A "glmvsd" object is retured. The components are:
 \item{glmvsd.index}{VSD value}
 \item{weight}{The weight for each candidate model}
 \item{model_check}{If \code{model_check} is not specified, the selected model based on the \code{penalty} will be provided.}
 \item{difference}{The number different variables in the candidate model from the \code{model_check} model.}
 \item{false.positive}{The upper VSD value of \code{model_check} model}
 \item{false.negative}{The lower VSD value of \code{model_check} model}
}
\references{

}
\examples{

x<-matrix(rnorm(50*8),ncol=8)
y<-rnorm(50,0,1.5)
fit<-glmvsd(x,y,n_train=40,n_rep=50,candidate="union",penalty="Lasso",weight_fun="ARM")
}

