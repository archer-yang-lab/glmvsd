glmvsd <- function(x, y, n_train = ceiling(n/2), no_rep = 100, model_check, 
    psi = 1, family = c("gaussian", "binomial"), method = c("union", 
        "customize"), candidate_models, weight_function = c("ARM", "BIC"), 
    prior = TRUE) {
    # check data and parameter
    family <- match.arg(family)
    method <- match.arg(method)
    weight_function <- match.arg(weight_function)
    y <- drop(y)
    y <- as.numeric(y)
    x <- as.matrix(x)
    p <- NCOL(x)
    n <- length(y)
    if (family == "binomial") {
        if (!all(y %in% c(0, 1))) 
            stop("There can only be 0 or 1 in y when using binomial family")
    }
    # if (family == "tweedie") {
    #     if (any(y < 0)) 
    #         stop("y must be nonzero when using Tweedie family")
    # }
    if (n != NROW(x)) 
        stop("x and y have different number of observations")
    if (n_train >= n) 
        stop("Training size must be less than the number of observations")
    if (missing(model_check)) 
        stop("User must provide a base model.")
    # use union option to compute candidate models
    if (method == "union") {
        if (family == "gaussian") 
            candidate_models <- gaussianfit(x, y)
        if (family == "binomial") 
            candidate_models <- binomialfit(x, y)
        # if (family == "tweedie") 
        #     candidate_models <- tweediefit(x, y)
    }
    if (method == "customize") {
        if (missing(candidate_models)) 
            stop("Users must supply a candidate model.")
        if (is.matrix(candidate_models) != TRUE) 
            stop("Supplied model must be a matrix.")
        if (NCOL(candidate_models) != NCOL(x)) 
            stop("Number of variables in candidate model and x does not match.")
        if (!all(as.numeric(candidate_models) %in% c(0, 1))) 
            stop("There can only be 0 or 1 in candidate_models")
    }
    # clean the candidate models
    candidate_models <- unique(candidate_models)
    rownames(candidate_models) <- NULL
    candidate_models <- candidate_models[order(rowSums(candidate_models)), 
        ]
    candidate_models <- candidate_models[rowSums(candidate_models) < 
        (n_train - 2), ]
    # compute weights
    if (family == "gaussian") {
        if (weight_function == "ARM") {
            fit <- lsARM(x = x, y = y, candidate_models = candidate_models, 
                n_train = n_train, no_rep = no_rep, psi = psi, prior = prior)
        }
        if (weight_function == "BIC") {
            fit <- lsBIC(x = x, y = y, candidate_models = candidate_models, 
                psi = psi, prior = prior)
        }
    }
    if (family == "binomial") {
        if (weight_function == "ARM") {
            fit <- logitARM(x = x, y = y, candidate_models = candidate_models, 
                n_train = n_train, no_rep = no_rep, psi = psi, prior = prior)
        }
        if (weight_function == "BIC") {
            fit <- logitBIC(x = x, y = y, candidate_models = candidate_models, 
                psi = psi, prior = prior)
        }
    }
    weight <- fit$weight
	# compute VSD etc
    TMP_matrix <- sweep(candidate_models, MARGIN = 2, model_check, "-")
    DIFF <- rowSums(abs(TMP_matrix))
    DIFF_minus <- rowSums(TMP_matrix == -1)
    DIFF_plus <- rowSums(TMP_matrix == 1)
    VSD <- weight %*% DIFF  # glmvsd value
    VSD_minus <- weight %*% DIFF_minus  #false positive
    VSD_plus <- weight %*% DIFF_plus  #false negative    
    object <- list(VSD = VSD, VSD_minus = VSD_minus, VSD_plus = VSD_plus, 
        weight = weight, difference = DIFF, candidate_models_cleaned = candidate_models)
    class(object) <- "glmvsd"
    object
}