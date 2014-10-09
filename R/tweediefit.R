tweediefit <- function(x, y) {
    # fit lasso model
    lassofit <- HDtweedie(x = x, y = y, alpha = 1, maxit = 1e+06)
    # fit elastic net model
    enetfit <- HDtweedie(x = x, y = y, alpha = 0.5, maxit = 1e+06)
    # expand x matrix using 5 basis splines for group lasso
    ex <- expnd(x)
    gp <- rep(1:NCOL(x), each = 5)
    # fit group lasso
    glassofit <- HDtweedie(x = ex, y = y, group = gp, alpha = 1, maxit = 1e+06)
    # fit group elastic net
    genetfit <- HDtweedie(x = ex, y = y, group = gp, alpha = 0.5, maxit = 1e+06)
    # fit tweedie GLM model
    tweediefit <- glm.fit(x = x, y = y, family = tweedie(var.power = 1.5, 
        link.power = 0))
    # stepwise selection both side
    bothfit <- step.tweedie(tweediefit, y, x, 1.5, grpscp = 1:NCOL(x), 
        bn = NCOL(x), ix = 1:NCOL(x), iy = 1:NCOL(x), direction = "both", 
        steps = 100, pvalf = 0.05, pvalb = 0.1, trace = 0)
    # backward selection
    backfit <- step.tweedie(tweediefit, y, x, 1.5, grpscp = 1:NCOL(x), 
        bn = NCOL(x), ix = 1:NCOL(x), iy = 1:NCOL(x), direction = "backward", 
        steps = 100, pvalf = 0.05, pvalb = 0.1, trace = 0)
    # forward selection
    forefit <- step.tweedie(tweediefit, y, x, 1.5, grpscp = 1:NCOL(x), 
        bn = NCOL(x), ix = 1:NCOL(x), iy = 1:NCOL(x), direction = "forward", 
        steps = 100, pvalf = 0.05, pvalb = 0.1, trace = 0)
    # gether lasso and enet result
    lasso.path <- as.matrix(lassofit$beta)
    enet.path <- as.matrix(enetfit$beta)
    # gether glasso and genet result
    gind <- seq(1, NCOL(x) * 5, by = 5)
    glasso.path <- as.matrix(glassofit$beta[gind, ])
    genet.path <- as.matrix(genetfit$beta[gind, ])
    beta.path <- t(cbind(lasso.path, enet.path, glasso.path, genet.path))
    candidate_models <- (1 - (beta.path == 0))
    # gether stepwise selection result
    zerov <- rep(0, NCOL(x))
    (zerov[order(bothfit$final.xindex)] <- 1)
    b1 <- zerov
    zerov <- rep(0, NCOL(x))
    (zerov[order(backfit$final.xindex)] <- 1)
    b2 <- zerov
    zerov <- rep(0, NCOL(x))
    (zerov[order(forefit$final.xindex)] <- 1)
    b3 <- zerov
    candidate_models <- rbind(candidate_models, b1, b2, b3)
    candidate_models
}