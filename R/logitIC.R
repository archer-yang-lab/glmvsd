logitIC <- function(x, y, candidate_models, psi,
                    type = c("BIC", "AIC"), prior = TRUE) {
    p <- NCOL(x)
    n <- length(y)
	type <- match.arg(type)
    n_mo <- NROW(candidate_models)
    sk <- rowSums(candidate_models)
    ik <- rep(NA, n_mo)
    if (any(candidate_models[1, ] == 1)) {
      for (i in seq(n_mo)) {
        glmfit <- logistf(y ~ x[, candidate_models[i, ] == 1])
        if (type == "BIC") ik[i] <- extractAIC(glmfit, k=log(n))[2]
        else ik[i] <- extractAIC(glmfit)[2]
      }   
    } else {
      ik[1] <- -2 * log(prod(mean(y)^y * (1 - mean(y))^(1 - y))) + sk[1] * log(n)
      for (i in seq(2, n_mo)) {
        glmfit <- logistf(y ~ x[, candidate_models[i, ] == 1])
        if(type == "BIC") ik[i] <- extractAIC(glmfit, k=log(n))[2]
        else ik[i] <- extractAIC(glmfit)[2]
      }
    }
    if (prior == TRUE) {
      ck <- ck_compute(n_mo, sk, p)
      ik <- ik + 2 * psi * ck
    }
    weight <- exp(-ik/2)/sum(exp(-ik/2))
    list(weight = weight)
}
