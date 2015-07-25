logitIC <- function(x, y, candidate_models, n_train, no_rep, psi, 
	type = c("BIC", "AIC"), prior = TRUE) {
    p <- NCOL(x)
    n <- length(y)
    n_mo <- NROW(candidate_models)
    sk <- rowSums(candidate_models)
    ik <- rep(NA, n_mo)
    if (any(candidate_models[1, ] == 1)) {
        for (i in 1:n_mo) {
            glmfit <- glm(y ~ x[, candidate_models[i, ] == 1], family = "binomial")
            if(type == "BIC") ik[i] <- glmfit$deviance + sk[i] * log(n)
			else ik[i] <- glmfit$deviance + sk[i] * 2
        }
    } else {
        ik[1] <- -2 * log(prod(mean(y)^y * (1 - mean(y))^(1 - y))) + sk[1] * 
            log(n)
        for (i in 2:n_mo) {
            glmfit <- glm(y ~ x[, candidate_models[i, ] == 1], family = "binomial")
            if(type == "BIC") ik[i] <- glmfit$deviance + sk[i] * log(n)
			else ik[i] <- glmfit$deviance + sk[i] * 2
        }
    }
    if (prior == TRUE) {
        ck <- ck_compute(n_mo, sk, p)
        ik <- ik + 2 * psi * ck
    }
    weight <- exp(-ik/2)/sum(exp(-ik/2))
    list(weight = weight)
}
