library(glmvsd)
library(HDtweedie)
library(splines)
library(tweedie)
library(statmod)
source("Tweedie_stepwise.R", chdir = TRUE)
source("tools.R", chdir = TRUE)

n = 300 
p = 8

b <- c(1,1,1,-3*sqrt(2)/2) 
x=matrix(rnorm(n*p, mean=0, sd=1), n, p) 
feta=x[, 1:4]%*%b 
y = sapply(exp(feta), rand_tweedie)


family = "tweedie"
# user provide a model to be checked
model_check <- c(0,1,1,1,0,0,0,1)

candidate_models = 
rbind(c(0,0,0,0,0,0,0,1), 
c(0,1,0,0,0,0,0,1), c(1,1,1,1,0,0,0,0), 
c(0,1,1,0,0,0,0,1), c(1,1,0,1,1,0,0,0), 
c(1,1,0,0,1,0,0,0), c(0,0,0,0,0,0,0,0),
c(1,1,1,1,1,0,0,0))


n_train = ceiling(n/2)
n_rep = 100
candidate_models <- unique(candidate_models)
rownames(candidate_models) <- NULL
candidate_models <- candidate_models[order(rowSums(candidate_models)), ]
candidate_models <- candidate_models[rowSums(candidate_models) < (n_train-2), ]

psi = 1
i = 3
j = 30
prior = TRUE