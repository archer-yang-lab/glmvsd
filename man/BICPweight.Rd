\name{BICPweight}
\alias{BICPweight}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
BIC weight with un-normalized prior
}
\description{
Calculate un-normalized prior weight based on BIC for each candidate model.
}
\usage{
BICPweight(x, y, cand.mod, psi)
}

\arguments{
  \item{x}{
Matrix of predictors
}
  \item{y}{
Response vector
}
  \item{cand.mod}{
The list of candidate models for calculating the weight. The candidate model is the form of matrix with indicators of 0 and 1. 1 indicates the nonzero variables and 0 indicates the null variables. 
}
\item{psi}{
Size of un-normalized prior. If \code{psi} is missing, the default value of \code{psi}=1
}
}

\details{
Given x, y and candidate models with size less than \code{n.train}, \code{BICPweight} calculates the BIC weight with un-normalized prior. \code{psi} controls the size of un-normalized prior. 
}
\value{
The weight and the candidate model are returned. 
 \item{weight.BIC.Prior}{the BIC weight with un-normalized prior}
 \item{cand.model}{The list of the candidate model}
}

