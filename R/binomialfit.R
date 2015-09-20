binomialfit <- function(x, y) {
  lassofit <- glmnet(x = x, y = y, family = "binomial", alpha = 1, maxit = 1e+06)
  # enetfit1 <- glmnet(x = x, y = y, family = "binomial", alpha = 0.1, maxit = 1e+06)
  # enetfit2 <- glmnet(x = x, y = y, family = "binomial", alpha = 0.3, maxit = 1e+06)
  # enetfit3 <- glmnet(x = x, y = y, family = "binomial", alpha = 0.9, maxit = 1e+06)
  scadfit <- ncvreg(X = x, y = y, family = "binomial", penalty = "SCAD", 
                    warn = FALSE, max.iter = 1e+04)
  mcpfit <- ncvreg(X = x, y = y, family = "binomial", penalty = "MCP", 
                   warn = FALSE, max.iter = 1e+04)
  lasso.path <- as.matrix(lassofit$beta)
  # enet.path1 <- as.matrix(enetfit1$beta)
  # enet.path2 <- as.matrix(enetfit2$beta)
  # enet.path3 <- as.matrix(enetfit3$beta)
  scad.path <- as.matrix(scadfit$beta[-1, ])
  mcp.path <- as.matrix(mcpfit$beta[-1, ])
  # beta.path <- t(cbind(lasso.path, enet.path1, enet.path2, enet.path3, scad.path, mcp.path))
  beta.path <- t(cbind(lasso.path, scad.path, mcp.path))
  candidate_models <- (1 - (beta.path == 0))
  candidate_models
}