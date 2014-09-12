\name{ARMweight}
\alias{ARMweight}
\title{
ARM weight
}
\description{
Calculate weight based on ARM for each candidate model
}
\usage{
ARMweight(x, y, n.sim, cand.mod, n.train)
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
}
\details{
Given x, y and candidate models, \code{ARMweight} calculates the ARM weight for each candidate model. 
}
\value{
The weight and the candidate model are returned. 
 \item{weight.ARM}{the ARM weight.}
 \item{cand.model}{The list of the candidate model.}
}
