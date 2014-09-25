gaussianfit <- function(x, y) {
    lassofit <- glmnet(x = x, y = y, alpha = 1, maxit = 1e+06)
    scadfit <- ncvreg(X = x, y = y, family = "gaussian", penalty = "SCAD", 
        max.iter = 1e+06)
    mcpfit <- ncvreg(X = x, y = y, family = "gaussian", penalty = "MCP", 
        max.iter = 1e+06)
    lasso.path <- as.matrix(lassofit$beta)
    scad.path <- as.matrix(scadfit$beta[-1, ])
    mcp.path <- as.matrix(mcpfit$beta[-1, ])
    beta.path <- t(cbind(lasso.path, scad.path, mcp.path))
    candidate_models <- (1 - (beta.path == 0))
    candidate_models
}