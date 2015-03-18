ck_compute <- function(n_mo, sk, p) {
	ck <- rep(NA, n_mo)
	if (sk[1] == 0) {
	    ck[1] <- 2 * log(sk[1] + 2)/choose(p, sk[1])
	    ck[2:n_mo] <- sk[2:n_mo] * log(exp(1) * p/sk[2:n_mo]) + 2 * 
	        log(sk[2:n_mo] + 2)
	} else {
	    ck <- sk * log(exp(1) * p/sk) + 2 * log(sk + 2)
	}
	ck
}


