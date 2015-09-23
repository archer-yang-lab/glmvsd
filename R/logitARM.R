logitARM <- function(x, y, candidate_models, n_train, no_rep, psi, prior = TRUE, reduce_bias = FALSE) {
  p <- NCOL(x)
  n <- length(y)
  n_mo <- NROW(candidate_models)
  sk <- rowSums(candidate_models)
  
  wt_calc <- function(rep_id) {
    wts <- rep(NA, n_mo)
    tindex <- sample(n, n_train, replace = FALSE)
    if (any(candidate_models[1, ] == 1)) {
      for (j in seq(n_mo)) {
        varindex <- (candidate_models[j, ] == 1)
        glmfit <- if(reduce_bias==TRUE) brglm(y[tindex] ~ x[tindex, varindex], family = binomial) else glm(y[tindex] ~ x[tindex, varindex], family = binomial)
        if (any(is.na(glmfit$coef))) {
          wts[j] <- 0
        } else {
          gk <- as.vector(cbind(1, x[-tindex, varindex]) %*% glmfit$coef)
          fk <- ifelse(gk < 0, exp(gk)/(1 + exp(gk)), 1/(1 + exp(-gk)))
          wts[j] <- prod(fk^(y[-tindex]) * (1 - fk)^(1 - y[-tindex]))
        }  
      }
    } else {
      wts[1] <- prod(mean(y[tindex])^(y[-tindex]) * (1 - mean(y[tindex]))^(1 - y[-tindex]))
      for (j in seq(2, n_mo)) {
        varindex <- (candidate_models[j, ] == 1)
        glmfit <- if(reduce_bias==TRUE) brglm(y[tindex] ~ x[tindex, varindex], family = binomial) else glm(y[tindex] ~ x[tindex, varindex], family = binomial) 
        if(any(is.na(glmfit$coef))) {
          wts[j] <- 0
        } else {
          gk <- as.vector(cbind(1, x[-tindex, varindex]) %*% glmfit$coef)
          fk <- ifelse(gk < 0, exp(gk)/(1 + exp(gk)), 1/(1 + exp(-gk)))
          wts[j] <- prod(fk^(y[-tindex]) * (1 - fk)^(1 - y[-tindex]))
        }
      }
    }
    return(wts)
  }
  library(parallel)
  w_num <- matrix(unlist(mclapply(seq(no_rep), wt_calc)), nrow = no_rep, ncol = n_mo, byrow = TRUE)
  
  if (prior == TRUE) {
    ck <- ck_compute(n_mo, sk, p)
    w_num <- sweep(w_num, MARGIN = 2, exp(-psi * ck), "*")
  }
  weight <- colMeans(w_num/rowSums(w_num))
  list(weight = weight)
}
