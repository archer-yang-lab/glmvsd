lsARM <- function(x, y, candidate_models, n_train, n_rep, psi, prior = TRUE) {
    p <- NCOL(x)
    n <- length(y)
    n_mo <- NROW(candidate_models)
    sk <- rowSums(candidate_models)
    dk <- matrix(NA, n_rep, n_mo)
    sigma_k <- matrix(NA, n_rep, n_mo)
    for (i in 1:n_rep) {
        tindex <- sample(n, n_train, replace = F)
        if (any(candidate_models[1, ] == 1)) {
            for (j in 1:n_mo) {
                LSL <- lm(y[tindex] ~ x[tindex, candidate_models[j, ] == 
                  1])
                dk[i, j] <- sum((y[-tindex] - cbind(1, x[-tindex, candidate_models[j, 
                  ] == 1]) %*% LSL$coef)^2)
                sigma_k[i, j] <- summary(LSL)$sigma
            }
        } else {
            dk[i, 1] <- sum((y[-tindex] - mean(y[tindex]))^2)
            sigma_k[i, 1] <- sd(y[tindex])
            for (j in 2:n_mo) {
                LSL <- lm(y[tindex] ~ x[tindex, candidate_models[j, ] == 
                  1])
                dk[i, j] <- sum((y[-tindex] - cbind(1, x[-tindex, candidate_models[j, 
                  ] == 1]) %*% LSL$coef)^2)
                sigma_k[i, j] <- summary(LSL)$sigma
            }
        }
    }
    lw_num <- (-n/2) * log(sigma_k) - (1/sqrt(sigma_k)) * dk/2
    if (prior == TRUE) {
        ck <- ck_compute(n_mo, sk, p)
        lw_num <- sweep(lw_num, MARGIN = 2, psi * ck, "-")
    }
    lw_num <- sweep(lw_num, MARGIN = 1, apply(lw_num, 1, max), "-")
    w_num <- apply(lw_num, c(1, 2), function(x) ifelse(abs(x) > 700, 0, 
        exp(x)))
    weight <- apply(w_num/rowSums(w_num), 2, mean)
    list(weight = weight)
}