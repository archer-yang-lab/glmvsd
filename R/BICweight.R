BICPweight <- function(x, y, candidate_model, psi, prior = TRUE) {
    p <- NCOL(x)
    n <- length(y)
    n_mo <- NROW(candidate_model)
    sk <- rowSums(candidate_model)
    ik <- rep(NA, n_mo)
    if (any(candidate_model[1, ] == 1)) {
        for (i in 1:n_mo) {
            LSL <- lm(y ~ x[, candidate_model[i, ] == 1])
            rss <- crossprod(summary(LSL)$res, summary(LSL)$res)
            ik[i] <- n * log(rss/n) + sk[i] * log(n)
        }
    } else {
        rss <- sum((y - mean(y))^2)
        ik[1] <- n * log(rss/n) + sk[1] * log(n)
        for (i in 2:n_mo) {
            LSL <- lm(y ~ x[, candidate_model[i, ] == 1])
            rss <- crossprod(summary(LSL)$res, summary(LSL)$res)
            ik[i] <- n * log(rss/n) + sk[i] * log(n)
        }
    }
    if (prior == TRUE) {
        ck <- rep(NA, n_mo)
        if (sk[1] == 0) {
            ck[1] <- 2 * log(sk[1] + 2)/choose(p, sk[1])
            ck[2:n_mo] <- sk[2:n_mo] * log(exp(1) * p/sk[2:n_mo]) + 2 * 
                log(sk[2:n_mo] + 2)
        } else {
            ck <- sk * log(exp(1) * p/sk) + 2 * log(sk + 2)
        }
        ik <- ik + psi * ck
    }
    ik <- ik - min(ik)
    weight <- exp(-ik/2)/sum(exp(-ik/2))
    outlist <- list(weight = weight)
}