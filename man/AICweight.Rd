\name{AICweight}
\alias{AICweight}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
AIC weight
}
\description{
Calculate weight based on AIC
}
\usage{
AICweight(x, y, cand.mod)
}

\arguments{
  \item{x}{Matrix of predictors}
  \item{y}{Response}
  \item{cand.mod}{The list of candidate models for calculating the weight. The candidate model is the form of matrix with indicators of 0 and 1. 1 indicates the nonzero variables and 0 indicates the null variables. }
}
\details{
Given x, y and candidate models with size less than \code{n.train}, \code{AICweight} calculates the AIC weight for each candidate model. 
}
\value{
The weight and the candidate model are returned. 
 \item{weight.AIC}{the AIC weight.}
 \item{cand.model}{The list of the candidate model.}
}

