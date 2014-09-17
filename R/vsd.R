vsd <- function(x, y, n_train = ceiling(n/2), n_rep = 100, k = 10, base_model, 
    psi = 1, candidate = c("union", "supplied"), candidate_model, weight_fun = c("ARM", 
        "ARM.Prior", "AIC", "AIC.Prior", "BIC", "BIC.Prior")) {
    # check data and parameter
    candidate <- match.arg(candidate)
    weight_fun <- match.arg(weight_fun)
    y <- drop(y)
    x <- as.matrix(x)
    p <- NCOL(x)
    n <- length(y)
    if (n != NROW(x)) 
        stop("x and y have different number of observations")
    if (n_train >= n) 
        stop("Training size must be less than the number of observations")
    if (missing(base_model)) 
        stop("User must provide a base model.")
    
	# use union option to compute candidate models
    if (candidate == "union") {
        lassofit <- glmnet(x = x, y = y, alpha = 1, maxit = 1e+08)
        scadfit <- ncvreg(X = x, y = y, family = "gaussian", penalty = "SCAD", 
            max.iter = 1e+07)
        mcpfit <- ncvreg(X = x, y = y, family = "gaussian", penalty = "MCP", 
            max.iter = 1e+07)
        lasso.path <- as.matrix(lassofit$beta)
        scad.path <- as.matrix(scadfit$beta[-1, ])
        mcp.path <- as.matrix(mcpfit$beta[-1, ])
        beta.path <- t(cbind(lasso.path, scad.path, mcp.path))
        ind.path <- (1 - (beta.path == 0))
        candidate_model <- unique(ind.path)
        rownames(candidate_model) <- NULL
    }
    if (candidate == "supplied") {
        if (missing(candidate_model)) 
            stop("Users must supply a candidate model.")
        if (is.matrix(candidate_model) != TRUE) 
            stop("Supplied model must be a matrix.")
        if (NCOL(candidate_model) != NCOL(x)) 
            stop("Number of variables in candidate model and x does not match.")
        ######################## add 0-1 check
    }
    # compute weights    
    if (weight_fun == "ARM") {
        fit <- ARMweight(x = x, y = y, n_rep = n_rep, candidate_model = candidate_model, 
            n_train = n_train, prior=FALSE)
        weight <- fit$weight
        ending_candidate_model <- fit$ending_candidate_model
    }
    
    if (weight_fun == "ARM.Prior") {
        fit <- ARMweight(x = x, y = y, n_rep = n_rep, candidate_model = candidate_model, 
            psi = psi, n_train = n_train, prior=TRUE)
        weight <- fit$weight
        ending_candidate_model <- fit$ending_candidate_model
    }
    
    if (weight_fun == "AIC") {
        fit <- AICweight(x = x, y = y, candidate_model = candidate_model)
        weight <- fit$weight
        ending_candidate_model <- fit$ending_candidate_model
    }
    
    if (weight_fun == "AIC.Prior") {
        fit <- AICPweight(x = x, y = y, candidate_model = candidate_model, 
            psi = psi)
        weight <- fit$weight
        ending_candidate_model <- fit$ending_candidate_model
    }
    
    if (weight_fun == "BIC") {
        fit <- BICweight(x = x, y = y, candidate_model = candidate_model)
        weight <- fit$weight
        ending_candidate_model <- fit$ending_candidate_model
    }
    
    if (weight_fun == "BIC.Prior") {
        fit <- BICPweight(x = x, y = y, candidate_model = candidate_model, 
            psi = psi)
        weight <- fit$weight
        ending_candidate_model <- fit$ending_candidate_model
    }
    TMP_matrix <- sweep(ending_candidate_model, MARGIN = 2, base_model, "-")
	DIFF <- rowSums(abs(TMP_matrix))
    DIFF_minus <- rowSums(TMP_matrix == -1)
    DIFF_plus <- rowSums(TMP_matrix == 1)
    VSD <- weight %*% DIFF  # vsd value
    VSD_minus <- weight %*% DIFF_minus  #false positive
    VSD_plus <- weight %*% DIFF_plus  #false negative    
    object <- list(VSD = VSD, VSD_minus = VSD_minus, VSD_plus = VSD_plus, 
        weight = weight, difference = DIFF)
    class(object) <- "vsd"
    object
}


