lsARMcore <- function(x, y, candidate_models, n_train, no_rep) {
	p <- NCOL(x)
    n <- length(y)
    n_mo <- NROW(candidate_models)
    sk <- rowSums(candidate_models)
    dk <- rep(NA, n_mo)
    sigma_k <- rep(NA, n_mo)
	tindex <- sample(n, n_train, replace = F)
        if (any(candidate_models[1, ] == 1)) {
            for (j in 1:n_mo) {
                LSL <- lm(y[tindex] ~ x[tindex, candidate_models[j, ] == 
                  1])
                sigma_k[j] <- summary(LSL)$sigma
				if(any(is.na(LSL$coef))){
					dk[j] <- Inf
					}else{
				        dk[j] <- sum((y[-tindex] - cbind(1, x[-tindex, candidate_models[j, 
		                  ] == 1]) %*% LSL$coef)^2)
					}
            }
        } else {
            dk[1] <- sum((y[-tindex] - mean(y[tindex]))^2)
            sigma_k[1] <- sd(y[tindex])
            for (j in 2:n_mo) {
                LSL <- lm(y[tindex] ~ x[tindex, candidate_models[j, ] == 
                  1])
	            sigma_k[j] <- summary(LSL)$sigma		
				if(any(is.na(LSL$coef))){
					dk[j] <- Inf
					}else{
				        dk[j] <- sum((y[-tindex] - cbind(1, x[-tindex, candidate_models[j, 
		                  ] == 1]) %*% LSL$coef)^2)
					}
            }
        }
	list(sigma_k=sigma_k,dk=dk)
}