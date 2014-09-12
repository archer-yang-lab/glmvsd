\name{ARMPweight}
\alias{ARMPweight}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
ARM weight with un-normalized prior
}
\description{
Calculate un-normalized prior weight based on ARM for each candidate model.
}
\usage{
ARMPweight(x, y, n.sim,cand.mod,n.train,psi)
}

\arguments{
  \item{x}{
Matrix of predictors
}
  \item{y}{
Response
}
\item{n.sim}{n.sim: number of replication in data split}
\item{cand.mod}{
The list of candidate models for calculating the weight. The candidate model is the form of matrix with indicators of 0 and 1. 1 indicates the nonzero variables and 0 indicates the null variables. 
}
\item{n.train}{Size of training set when the weight function is ARM and ARM with prior. If \code{n.train} is not specified, then \code{n.train}=ceiling(n)}
  \item{psi}{
Size of un-normalized prior
}
}
\details{
Given x, y, psi and candidate models, \code{ARMPweight} calculates the ARM weight with un-normalized prior. \code{psi} controls the size of un-normalized prior. 
}
\value{
The weight and the candidate model are returned. 
 \item{weight.ARM.Prior}{the ARM weight with un-normalized prior}
 \item{cand.model}{The list of the candidate model}
}

