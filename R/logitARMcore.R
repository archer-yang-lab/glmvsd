logitARMcore <- function(x, y, candidate_models, n_train, no_rep) {
	p <- NCOL(x)
    n <- length(y)
	n_mo <- NROW(candidate_models)
	tindex <- sample(n, n_train, replace = F)
	w_num <- rep(NA, n_mo)
	if (any(candidate_models[1, ] == 1)) {
	    for (j in 1:n_mo) {
	        glmfit <- glm(y[tindex] ~ x[tindex, candidate_models[j, 
	          ] == 1], family = "binomial", control = list(maxit = 1e7))
			if(any(is.na(glmfit$coef))){w_num[j] <- 0
				}else{
					gk <- cbind(1, x[-tindex, candidate_models[j, ] == 1]) %*% 
	                  glmfit$coef
	                fk <- exp(gk)/(exp(gk) + 1)
	                fk[gk > 700] <- 1
	                fk[gk < -700] <- 0
	                w_num[j] <- prod(fk^y[-tindex] * (1 - fk)^(1 - y[-tindex]))
			}  
	    }
	} else {
	    w_num[1] <- prod(mean(y[tindex])^y[-tindex] * (1 - mean(y[tindex]))^(1 - 
	        y[-tindex]))
	    for (j in 2:n_mo) {
	        glmfit <- glm(y[tindex] ~ x[tindex, candidate_models[j, 
	          ] == 1], family = "binomial", control = list(maxit = 1e7))
			if(any(is.na(glmfit$coef))){w_num[j] <- 0
				}else{
					gk <- cbind(1, x[-tindex, candidate_models[j, ] == 1]) %*% 
	                  glmfit$coef
	                fk <- exp(gk)/(exp(gk) + 1)
	                fk[gk > 700] <- 1
	                fk[gk < -700] <- 0
	                w_num[j] <- prod(fk^y[-tindex] * (1 - fk)^(1 - y[-tindex]))
			}
	    }
	}
	w_num
}
