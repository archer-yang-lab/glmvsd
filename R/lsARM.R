lsARM <- function(x, y, candidate_models, n_train, no_rep, psi, prior = TRUE) {
  p <- NCOL(x)
  n <- length(y)
  n_mo <- NROW(candidate_models)
  sk <- rowSums(candidate_models)
  
  wt_calc <- function(rep_id) {
    dk <- rep(NA, n_mo)
    sigma_k <- rep(NA, n_mo)
    tindex <- sample(n, n_train, replace = FALSE)
    if (any(candidate_models[1, ] == 1)) {
      for (j in seq(n_mo)) {
        varindex <- (candidate_models[j, ] == 1)
        LSL <- lm(y[tindex] ~ x[tindex, varindex])
        sigma_k[j] <- summary(LSL)$sigma
        if(any(is.na(LSL$coef))) {
          dk[j] <- Inf
        } else {
          dk[j] <- sum((y[-tindex] - cbind(1, x[-tindex, varindex]) %*% LSL$coef)^2)
        }
      }
    } else {
      dk[1] <- sum((y[-tindex] - mean(y[tindex]))^2)
      sigma_k[1] <- sd(y[tindex])
      for (j in seq(2, n_mo)) {
        varindex <- (candidate_models[j, ] == 1)
        LSL <- lm(y[tindex] ~ x[tindex, varindex])
        sigma_k[j] <- summary(LSL)$sigma  	
        if (any(is.na(LSL$coef))) {
          dk[j] <- Inf
        } else {
          dk[j] <- sum((y[-tindex] - cbind(1, x[-tindex, varindex]) %*% LSL$coef)^2)
        }
      }
    }
    return(c(dk, sigma_k))
  }
  library(parallel)
  ARM_ret <- matrix(unlist(mclapply(seq(no_rep), wt_calc)), nrow = no_rep, ncol = 2*n_mo, byrow = TRUE)
  dk <- ARM_ret[, seq(n_mo)]
  sigma_k <- ARM_ret[, seq(n_mo + 1, 2*n_mo)]
  
#   dk <- matrix(NA, no_rep, n_mo)
#   sigma_k <- matrix(NA, no_rep, n_mo)
#   for (i in 1:no_rep) {
#     tindex <- sample(n, n_train, replace = FALSE)
#     if (any(candidate_models[1, ] == 1)) {
#       for (j in 1:n_mo) {
#         LSL <- lm(y[tindex] ~ x[tindex, candidate_models[j, ] == 
#                                   1])
#         sigma_k[i, j] <- summary(LSL)$sigma
#         if(any(is.na(LSL$coef))) {
#           dk[i, j] <- Inf
#         } else {
#           dk[i, j] <- sum((y[-tindex] - cbind(1, x[-tindex, candidate_models[j, ] == 1]) %*% LSL$coef)^2)
#         }
#       }
#     } else {
#       dk[i, 1] <- sum((y[-tindex] - mean(y[tindex]))^2)
#       sigma_k[i, 1] <- sd(y[tindex])
#       for (j in 2:n_mo) {
#         LSL <- lm(y[tindex] ~ x[tindex, candidate_models[j, ] == 1])
#         sigma_k[i, j] <- summary(LSL)$sigma		
#         if (any(is.na(LSL$coef))) {
#           dk[i, j] <- Inf
#         } else {
#           dk[i, j] <- sum((y[-tindex] - cbind(1, x[-tindex, candidate_models[j, ] == 1]) %*% LSL$coef)^2)
#         }
#       }
#     }
#   }

  lw_num <- (-n/2) * log(sigma_k) - ((sigma_k)^(-2)) * dk/2
  if (prior == TRUE) {
    ck <- ck_compute(n_mo, sk, p)
    lw_num <- sweep(lw_num, MARGIN = 2, psi * ck, "-")
  }
  lw_num <- sweep(lw_num, MARGIN = 1, apply(lw_num, 1, max), "-")
  w_num <- apply(lw_num, c(1, 2), function(x) ifelse(abs(x) > 700, 0, exp(x)))
  weight <- colMeans(w_num/rowSums(w_num))
  list(weight = weight)
}