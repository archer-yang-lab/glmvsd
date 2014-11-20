logitARM <- function(x, y, candidate_models, n_train, no_rep, psi, prior = TRUE) {
    p <- NCOL(x)
    n <- length(y)
    n_mo <- NROW(candidate_models)
    sk <- rowSums(candidate_models)
    w_num <- matrix(NA, no_rep, n_mo)
    for (i in 1:no_rep) {
        tindex <- sample(n, n_train, replace = F)
        if (any(candidate_models[1, ] == 1)) {
            for (j in 1:n_mo) {
                glmfit <- glm(y[tindex] ~ x[tindex, candidate_models[j, 
                  ] == 1], family = "binomial", control = list(maxit = 1e7))
                gk <- cbind(1, x[-tindex, candidate_models[j, ] == 1]) %*% 
                  glmfit$coef
                fk <- exp(gk)/(exp(gk) + 1)
                fk[gk > 700] <- 1
                fk[gk < -700] <- 0
                w_num[i, j] <- prod(fk^y[-tindex] * (1 - fk)^(1 - y[-tindex]))
            }
        } else {
            w_num[i, 1] <- prod(mean(y[tindex])^y[-tindex] * (1 - mean(y[tindex]))^(1 - 
                y[-tindex]))
            for (j in 2:n_mo) {
                glmfit <- glm(y[tindex] ~ x[tindex, candidate_models[j, 
                  ] == 1], family = "binomial", control = list(maxit = 1e7))
                gk <- cbind(1, x[-tindex, candidate_models[j, ] == 1]) %*% 
                  glmfit$coef
                fk <- exp(gk)/(exp(gk) + 1)
                fk[gk > 700] <- 1
                fk[gk < -700] <- 0
                w_num[i, j] <- prod(fk^y[-tindex] * (1 - fk)^(1 - y[-tindex]))
            }
        }
    }
    if (prior == TRUE) {
        ck <- ck_compute(n_mo, sk, p)
        w_num <- sweep(w_num, MARGIN = 2, exp(-psi * ck), "*")
    }
    weight <- apply(w_num/rowSums(w_num), 2, mean)
    list(weight = weight)
}
