lsARM <- function(x, y, candidate_models, 
					n_train, no_rep, psi, 
					prior = TRUE, parallel) {
    p <- NCOL(x)
    n <- length(y)
    n_mo <- NROW(candidate_models)
    sk <- rowSums(candidate_models)
    dk <- matrix(NA, no_rep, n_mo)
    sigma_k <- matrix(NA, no_rep, n_mo)
	if(parallel){
		outlist = foreach(i = seq(no_rep), .packages = c("glmvsd")) %dopar%{
			lsARMcore(x, y, candidate_models, n_train, no_rep)
		}
		res <- melt(outlist)
		sigma_k <- matrix(res$value[res$L2=="sigma_k"],no_rep,n_mo,byrow=TRUE)
		dk <- matrix(res$value[res$L2=="dk"],no_rep,n_mo,byrow=TRUE)
	}else{
		for (i in 1:no_rep) {
			res <- lsARMcore(x, y, candidate_models, n_train, no_rep)
	        sigma_k[i, ] <- res$sigma_k
			dk[i, ] <- res$dk
	    }
	}
    lw_num <- (-n/2) * log(sigma_k) - ((sigma_k)^(-2)) * dk/2
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