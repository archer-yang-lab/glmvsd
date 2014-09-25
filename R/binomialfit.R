binomialfit <- function(x, y) {
    lassofit <- glmnet(x = x, y = y, family = "binomial", alpha = 1, maxit = 1e+06)
    mcpfit <- ncvreg(X = x, y = y, family = "binomial", penalty = "MCP", 
        warn = FALSE, max.iter = 1e+06)
    lasso.path <- as.matrix(lassofit$beta)
    mcp.path <- as.matrix(mcpfit$beta[-1, ])
    beta.path <- t(cbind(lasso.path, mcp.path))
    candidate_models <- (1 - (beta.path == 0))
    candidate_models
}