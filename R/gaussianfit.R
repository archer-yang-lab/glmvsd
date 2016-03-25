gaussianfit <- function(x, y) {
  m1 <- glmnet(x = x, y = y, family = "gaussian", alpha = 1, maxit = 1e+06)
  m2 <- ncvreg(X = x, y = y, family = "gaussian", penalty = "SCAD",
                    warn = FALSE, max.iter = 1e+04)
  m3 <- ncvreg(X = x, y = y, family = "gaussian", penalty = "MCP",
                   warn = FALSE, max.iter = 1e+04)
  m1.path <- as.matrix(m1$beta)
  m2.path <- as.matrix(m2$beta[-1, ])
  m3.path <- as.matrix(m3$beta[-1, ])
  beta.path <- t(cbind(m1.path, m2.path, m3.path))
  candidate_models <- (1 - (beta.path == 0))
  candidate_models
}