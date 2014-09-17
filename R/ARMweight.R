ARMPweight <- function(x, y, n_rep = 100, psi = 1, candidate_model, n_train = ceiling(n/2), 
    prior = TRUE) {
    y <- drop(y)
    x <- as.matrix(x)
    p <- NCOL(x)
    n <- length(y)
    
    if (n != NROW(x)) 
        stop("x and y have different number of observations")
    if (n_train >= n) 
        stop("Training size must be less than the number of observations")
    
    if (missing(candidate_model)) 
        stop("Users must supply a candidate model.")
    if (is.matrix(candidate_model) != TRUE) 
        stop("Supplied model must be a matrix.")
    if (NCOL(candidate_model) != NCOL(x)) 
        stop("Number of variables in candidate model and x does not match.")
    
    model.ordered <- candidate_model[order(rowSums(candidate_model)), ]
    model <- model.ordered[rowSums(model.ordered) < n_train, ]
    n_mo <- NROW(model)
    d1 <- matrix(NA, n_rep, n_mo)
    s1 <- matrix(NA, n_rep, n_mo)
    
    for (i in 1:n_rep) {
        tindex <- sample(n, n_train, replace = F)
        if (any(model[1, ] == 1)) {
            for (j in 1:n_mo) {
                LSL <- lm(y[tindex] ~ x[tindex, model[j, ] == 1])
                d1[i, j] <- sum((y[-tindex] - cbind(1, x[-tindex, model[j, 
                  ] == 1]) %*% LSL$coef)^2)
                s1[i, j] <- summary(LSL)$sigma
            }
        } else {
            d1[i, 1] <- sum((y[-tindex] - mean(y[tindex]))^2)
            s1[i, 1] <- sd(y[tindex])
            for (j in 2:n_mo) {
                LSL <- lm(y[tindex] ~ x[tindex, model[j, ] == 1])
                d1[i, j] <- sum((y[-tindex] - cbind(1, x[-tindex, model[j, 
                  ] == 1]) %*% LSL$coef)^2)
                s1[i, j] <- summary(LSL)$sigma
            }
        }
    }
    
    tmp_omit <- na.omit(t(s1))
    index_omit <- attributes(tmp_omit)$na.action
    s1 <- s1[, -index_omit]
    d1 <- d1[, -index_omit]
    model <- model[-index_omit, ]
    sk <- rowSums(model)
    lw_num <- (-n/2) * log(s1) - (1/sqrt(s1)) * d1/2
    if (prior == TRUE) {
        n_mo_new <- NCOL(d1)
        ck <- rep(NA, n_mo_new)
        if (sk[1] == 0) {
            ck[1] <- 2 * log(sk[1] + 2)/choose(p, sk[1])
            ck[2:n_mo_new] <- sk[2:n_mo_new] * log(exp(1) * p/sk[2:n_mo_new]) + 
                2 * log(sk[2:n_mo_new] + 2)
        } else {
            ck <- sk * log(exp(1) * p/sk) + 2 * log(sk + 2)
        }
        lw_num <- sweep(lw_num, MARGIN = 2, psi * ck, "-")
    }
    lw_num <- sweep(lw_num, MARGIN = 1, apply(lw_num, 1, max), "-")
    w_num <- apply(lw_num, c(1, 2), function(x) ifelse(abs(x) > 700, 0, 
        exp(x)))
    w_den <- rowSums(w_num)
    weight <- apply(w_num/w_den, 2, mean)
    list(weight = weight, ending_candidate_model = model)
}