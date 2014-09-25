##############################################################################
### Wei Qian 6-2-14 Created
### Tweedie_stepwise.R
### This file contains the R function to perform stepwise selection for 
### Tweedie's compound Poisson model.

### The main function:
### step.tweedie(object, y, x, rho, grpscp, bn, ix, iy, 
### direction = c("both","backward","forward"), steps = 100, pvalf = 0.05, 
### pvalb = 0.1, trace = 0)

### object: an object used as the initial model in the stepwise search (see tweedie package for creation of tweedie GLM object)
### y: response variable
### x: matrix of predictors
### rho: the power used for variance-mean relation
### grpscp: index for groups in the intial model
### bn: number of groups
### ix: first index for each group
### iy: last index for each group
### direction: the mode of stepwise search
### steps: the maximum number of steps to be considered
### pvalf: term entering p-value
### pvalb: term removing p-value
### trace: if 1, the result of each step is printed during the running of step.tweedie

##############################################################################


deviance <- function(y,mu,rho) {
	2.0*sum(y*mu^(1.0-rho)/(rho-1.0)+mu^(2.0-rho)/(2.0-rho)-y^(2.0-rho)/(rho-1.0)/(2.0-rho))
}

pf_safe <- function(devdiff, df1, estscale, df2) {
	if (df1 <= 0) return(1)
	else {
		a <- devdiff/df1/estscale
		a <- pf(a,df1,df2,lower.tail=FALSE)
		return(a)
	}
}

add1.tweedie <- function(object, y, x, rho, grpscp, bn, ix, iy) {
	
	toadd <- setdiff(1:bn,grpscp)
	ns <- length(toadd)
	if (ns == 0) return (NULL) 
	
	addgrp <- estscale <- dfdiff <- devdiff <- pfval <- numeric(ns)
	
	cnt <- 0
	for (tt in toadd) {
		# create the index for x matrix
		xindex <- numeric(0)
		for (jj in union(grpscp,tt)) {
			xindex <- c(xindex,ix[jj]:iy[jj])
		}
		x1 <- x[,xindex]
		
		### see if the new design matrix is rank deficient
		mr <- as.numeric(rankMatrix(as.matrix(x1)))
		if (mr < length(xindex)) {
			cnt <- cnt+1
			pfval[cnt] <- Inf
			
		} else {
		result <- tryCatch(
		{
		z <- glm(y~x1,family=tweedie(link.power=0,var.power=rho))
		1 
		}, error = function(err) {
			print(paste("Myerr:",err,sep=" "))
			return(0)
		}		
		)
		if (result==1) {
		cnt <- cnt+1
		addgrp[cnt] <- tt
		estscale[cnt] <- sum((y-z$fitted.values)^2/z$fitted.values^rho)/z$df.residual
		dfdiff[cnt] <- object$df.residual-z$df.residual
		devdiff[cnt] <- object$deviance-z$deviance
		pfval[cnt] <- pf_safe(devdiff[cnt],dfdiff[cnt],estscale[cnt],z$df.residual)
		} else {
			cnt <- cnt+1
			pfval[cnt] <- Inf
		}	
		}
	}
	
	data.frame(addgrp = addgrp, estscale = estscale, dfdiff = dfdiff, devdiff = devdiff, pfval = pfval)
}

drop1.tweedie <- function(object, y, x, rho, grpscp, bn, ix, iy) {
	todrop <- grpscp
	ns <- length(todrop)
	if (ns == 0) return (NULL)
	
	estscale <- sum((y-object$fitted.values)^2/object$fitted.values^rho)/object$df.residual
	
	dropgrp <- dfdiff <- devdiff <- pfval <- numeric(ns)
	
	cnt <- 0
	for (tt in todrop) {
		# create the index for x matrix
		if (length(setdiff(grpscp,tt))==0) { # only one group in object
			z <- glm(y~1,family=tweedie(link.power=0,var.power=rho))
			result <- 1
		} else {
			xindex <- numeric(0)
			for (jj in setdiff(grpscp,tt)) {
				xindex <- c(xindex,ix[jj]:iy[jj])
			}
			x1 <- x[,xindex]
			result <- tryCatch({
			z <- glm(y~x1,family=tweedie(link.power=0,var.power=rho))
			1
			}, error = function(err) {
				print(paste("Myerr:",err,sep=" "))
				return(0)
			}
			)
		}
		if (result==1) {
		cnt <- cnt+1
		dropgrp[cnt] <- tt
		dfdiff[cnt] <- z$df.residual-object$df.residual
		devdiff[cnt] <- z$deviance-object$deviance
		pfval[cnt] <- pf_safe(devdiff[cnt],dfdiff[cnt],estscale,object$df.residual)	
		} else {
			cnt <- cnt+1
			pfval[cnt] <- 0
		}
	}
	
	data.frame(dropgrp = dropgrp, dfdiff = dfdiff, devdiff = devdiff, pfval = pfval)
	
}

step.tweedie <- function(object, y, x, rho, grpscp, bn, ix, iy, direction = c("both","backward","forward"), steps = 100, pvalf = 0.05, pvalb = 0.1, trace = 0) {
	library(Matrix)
	library(tweedie)
	direction <- match.arg(direction)
	backward <- direction == "both" | direction == "backward"
	forward <- direction == "both" | direction == "forward"
	fit <- object
	
	models <- vector("list", steps)
	nm <- steps
	models[[1]] <- list(grpscp = grpscp, change = "", step = nm-steps)
	if (trace) {
		print(paste("############## STEP", 0, sep=" "))
		print(models[[nm-steps+1]])
	}
	
	xindex <- numeric(0)
	if (length(grpscp)>0) {
		for (jj in grpscp) {
			xindex <- c(xindex,ix[jj]:iy[jj])
		}
	}
	
	while(steps > 0) {
		steps <- steps - 1
		
		## find the toadd and todrop groups
		toadd <- setdiff(1:bn,grpscp)
		todrop <- grpscp
		
		change <- NULL
		if (backward && length(todrop)) {
			aod <- drop1.tweedie(fit, y, x, rho, grpscp, bn, ix, iy)
			loc.max <- which.max(aod$pfval)
			## drop one term that is least signifant
			if (aod$pfval[loc.max] >= pvalb) {
				tt <- aod$dropgrp[loc.max]
				change <- paste("-",tt,sep="")
				
				# refit the model to update fit
				grpscp <- setdiff(grpscp,tt)
				if (length(grpscp)==0) { 
					xindex <- numeric(0)
					fit <- glm(y~1,family=tweedie(link.power=0,var.power=rho))
				} else {
					xindex <- numeric(0)
					for (jj in grpscp) {
						xindex <- c(xindex,ix[jj]:iy[jj])
					}
					x1 <- x[,xindex]
					fit <- glm(y~x1,family=tweedie(link.power=0,var.power=rho))
				}				
				models[[nm-steps+1]] <- list(grpscp = grpscp, change = change, step = nm-steps)
				if (trace) {
					print(paste("############## STEP", nm-steps, sep=" "))
					print(models[[nm-steps+1]])
				}			
			}	
		}
		if (is.null(change)) {
			if (forward && length(toadd)) {
				aod <- add1.tweedie(fit, y, x, rho, grpscp, bn, ix, iy)
				loc.min <- which.min(aod$pfval)
				## add one most significant term
				if (aod$pfval[loc.min] < pvalf) {
					tt <- aod$addgrp[loc.min]
					change <- paste("+",tt,sep="")
					
					# refit the model to update fit
					grpscp <- union(grpscp,tt)
					xindex <- numeric(0)
					for (jj in grpscp) {
						xindex <- c(xindex,ix[jj]:iy[jj])
					}
					x1 <- x[,xindex]
					fit <- glm(y~x1,family=tweedie(link.power=0,var.power=rho))
					models[[nm-steps+1]] <- list(grpscp = grpscp, change = change, step = nm-steps)
					if (trace) {
						print(paste("############## STEP", nm-steps, sep=" "))
						print(models[[nm-steps+1]])
					}	
				} 
			}
		}   
		
		## break out of the loop if nothing changes
		if (is.null(change)) break	
	}
	return(list(models = models[1:(nm-steps)], final.fit = fit, final.grpscp = grpscp, final.xindex = xindex, direction = direction))
}

