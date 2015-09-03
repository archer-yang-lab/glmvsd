logitARM <- function(x, y, candidate_models, 
					n_train, no_rep, psi, 
					prior = TRUE, parallel) {
    p <- NCOL(x)
    n <- length(y)
    n_mo <- NROW(candidate_models)
    sk <- rowSums(candidate_models)
    w_num <- matrix(NA, no_rep, n_mo)
	if(parallel){
		outlist = foreach(i = seq(no_rep), .packages = c("glmvsd")) %dopar%{
			logitARMcore(x, y, candidate_models, n_train, no_rep)
		} 
		w_num <- matrix(unlist(outlist),no_rep,n_mo,byrow=TRUE)
	}else{
		for (i in 1:no_rep) {
	        w_num[i, ] <- logitARMcore(x, y, candidate_models, n_train, no_rep)
	    }
	}
    if (prior == TRUE) {
        ck <- ck_compute(n_mo, sk, p)
        w_num <- sweep(w_num, MARGIN = 2, exp(-psi * ck), "*")
    }
    weight <- apply(w_num/rowSums(w_num), 2, mean)
    list(weight = weight)
}
