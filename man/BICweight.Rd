\name{BICweight}
\alias{BICweight}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
BIC weight
}
\description{
Calculate weight based on BIC
}
\usage{
BICweight(x, y, cand.mod)
}

\arguments{
  \item{x}{Matrix of predictors}
  \item{y}{Response}
  \item{cand.mod}{The list of candidate models for calculating the weight. The candidate model is the form of matrix with indicators of 0 and 1. 1 indicates the nonzero variables and 0 indicates the null variables. }
}
\details{
Given x, y and candidate models with size less than \code{n.train}, \code{BICweight} calculates the BIC weight for each candidate model. 
}
\value{
The weight and the candidate model are returned. 
 \item{weight.BIC}{the BIC weight.}
 \item{cand.model}{The list of the candidate model.}
}

