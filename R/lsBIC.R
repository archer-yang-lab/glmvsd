lsBIC <- function(x, y, candidate_models, psi, prior = TRUE) {
    p <- NCOL(x)
    n <- length(y)
    n_mo <- NROW(candidate_models)
    sk <- rowSums(candidate_models)
    ik <- rep(NA, n_mo)
    if (any(candidate_models[1, ] == 1)) {
        for (i in 1:n_mo) {
            LSL <- lm(y ~ x[, candidate_models[i, ] == 1])
            rss <- crossprod(summary(LSL)$res, summary(LSL)$res)
            ik[i] <- n * log(rss/n) + sk[i] * log(n)
        }
    } else {
        rss <- sum((y - mean(y))^2)
        ik[1] <- n * log(rss/n) + sk[1] * log(n)
        for (i in 2:n_mo) {
            LSL <- lm(y ~ x[, candidate_models[i, ] == 1])
            rss <- crossprod(summary(LSL)$res, summary(LSL)$res)
            ik[i] <- n * log(rss/n) + sk[i] * log(n)
        }
    }
    if (prior == TRUE) {
        ck <- ck_compute(n_mo, sk, p)
        ik <- ik + psi * ck
    }
    ik <- ik - min(ik)
    weight <- exp(-ik/2)/sum(exp(-ik/2))
    outlist <- list(weight = weight)
}