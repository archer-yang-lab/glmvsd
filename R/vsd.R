vsd <- function(x, y, n_train = ceiling(n/2), n_rep = 100, model_check, 
    psi = 1, candidate = c("union", "supplied"), candidate_models, weight_fun = c("ARM", 
        "ARM.Prior", "BIC", "BIC.Prior")) {
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
    if (missing(model_check)) 
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
    }
    if (candidate == "supplied") {
        if (missing(candidate_models)) 
            stop("Users must supply a candidate model.")
        if (is.matrix(candidate_models) != TRUE) 
            stop("Supplied model must be a matrix.")
        if (NCOL(candidate_models) != NCOL(x)) 
            stop("Number of variables in candidate model and x does not match.")
        ######################## add 0-1 check
    }
    candidate_models <- unique(ind.path)
    rownames(candidate_models) <- NULL
    candidate_models <- candidate_models[order(rowSums(candidate_models)), ]
    # compute weights    
    if (weight_fun == "ARM") {
	    candidate_models <- candidate_models[rowSums(candidate_models) < (n_train-2), ]
        fit <- ARMweight(x = x, y = y, candidate_models = candidate_models, 
            n_train = n_train, n_rep = n_rep, psi = psi, prior=FALSE)
        weight <- fit$weight
    }
    
    if (weight_fun == "ARM.Prior") {
	    candidate_models <- candidate_models[rowSums(candidate_models) < (n_train-2), ]
        fit <- ARMweight(x = x, y = y, candidate_models = candidate_models, 
            n_train = n_train, n_rep = n_rep, psi = psi, prior=TRUE)
        weight <- fit$weight
    }
    
    if (weight_fun == "BIC") {
	    candidate_models <- candidate_models[rowSums(candidate_models) < (n-2), ]
        fit <- BICweight(x = x, y = y, candidate_models = candidate_models, psi = psi, prior=FALSE)
        weight <- fit$weight
    }
    
    if (weight_fun == "BIC.Prior") {
	    candidate_models <- candidate_models[rowSums(candidate_models) < (n-2), ]
        fit <- BICweight(x = x, y = y, candidate_models = candidate_models, 
            psi = psi, prior=TRUE)
        weight <- fit$weight
    }
    TMP_matrix <- sweep(candidate_models, MARGIN = 2, model_check, "-")
	DIFF <- rowSums(abs(TMP_matrix))
    DIFF_minus <- rowSums(TMP_matrix == -1)
    DIFF_plus <- rowSums(TMP_matrix == 1)
    VSD <- weight %*% DIFF  # vsd value
    VSD_minus <- weight %*% DIFF_minus  #false positive
    VSD_plus <- weight %*% DIFF_plus  #false negative    
    object <- list(VSD = VSD, VSD_minus = VSD_minus, VSD_plus = VSD_plus, 
        weight = weight, difference = DIFF, candidate_models_cleaned = candidate_models)
    class(object) <- "vsd"
    object
}

