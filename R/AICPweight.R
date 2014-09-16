AICPweight <- function(x, y, candidate_model, psi) {
    
    n <- length(y)
    p <- NCOL(x)
    
    if (is.matrix(x) == "FASLE") 
        stop("x must be matrix with n rows")
    if (is.vector(y) == "FALSE") 
        stop("y must be a vector")
    
    if (missing(candidate_model)) 
        stop("missing candidate model")
    if (!missing(psi)) 
        psi <- psi else psi <- 1
    
    cand.nonzero <- apply(candidate_model, 1, sum)
    o <- order(cand.nonzero)
    model <- candidate_model[o, ]
    nonzero <- apply(model, 1, sum)
    m <- dim(model)[1]
    
    prior <- rep(0, m)
    
    for (i in 1:m) {
        if (nonzero[i] == 0) {
            prior[i] <- 2 * log(nonzero[i] + 2)/choose(p, nonzero[i])
        }
        if (nonzero[i] != 0) {
            prior[i] <- nonzero[i] * log(exp(1) * p/nonzero[i]) + 2 * log(nonzero[i] + 
                2)
        }
    }
    
    AIC.prior <- rep(0, m)
    
    for (i in 1:m) {
        if (nonzero[i] == 0) {
            LSL <- lm(y ~ 1)
            rss <- sum(summary(LSL)$res^2)
            AIC.prior[i] <- n * log(rss/n) + nonzero[i] * 2 + psi * prior[i]
        }
        if (nonzero[i] != 0) {
            x.new <- matrix(x[, which(model[i, ] == 1)], ncol = nonzero[i])
            LSL <- lm(y ~ x.new)
            rss <- sum(summary(LSL)$res^2)
            AIC.prior[i] <- n * log(rss/n) + nonzero[i] * 2 + psi * prior[i]
        }
    }
    
    m2 <- NULL
    AIC.prior2 <- NULL
    
    for (i in 1:m) {
        if (is.na(AIC[i]) == FALSE & abs(AIC[i]) != Inf) {
            m2 <- rbind(m2, model[i, ])
            AIC.prior2 <- c(AIC.prior2, AIC.prior[i])
        }
    }
    
    AIC.new.prior <- AIC.prior2 - min(AIC.prior2)
    
    weight.AIC.prior <- round(exp(-AIC.new.prior/2)/sum(exp(-AIC.new.prior/2)), 
        5)
    
    outlist <- list(weight.AIC.Prior = weight.AIC.prior, ending_candidate_model = m2)
}