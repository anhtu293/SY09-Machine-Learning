library(MASS)
adq.app <- function(Xapp, zapp)
{
	n <- dim(Xapp)[1]
	p <- dim(Xapp)[2]
	g <- max(unique(zapp))

	param <- NULL
	param$MCov <- array(0, c(p,p,g))
	param$mean <- array(0, c(g,p))
	param$prop <- rep(0, g)

	for (k in 1:g)
	{
	        #index des elements de class k
		indk <- which(zapp==k)
                #calculer lq proportions
		prop <- length(indk) / n
		#calculer les vecteurs de moyennes (esperence)
		moyenne <- Xapp[indk,]
		moyenne <- apply(moyenne, 2, mean)
		#calculer la matrice de covariance
		covariance <- t(as.matrix(Xapp[indk,] - moyenne)) %*% (as.matrix(Xapp[indk,] - moyenne))
		covariance <- covariance / length(indk)
		covariance_corr <- covariance*(length(indk)/(length(indk) - 1))
		#affecter les resultats
		param$MCov[,,k] <- covariance_corr
		param$mean[k,] <- moyenne
		param$prop[k] <- prop
	}

	param
}

adl.app <- function(Xapp, zapp)
{
	n <- dim(Xapp)[1]
	p <- dim(Xapp)[2]
	g <- max(unique(zapp))

	param <- NULL
	MCov <- array(0, c(p,p))
	param$MCov <- array(0, c(p,p,g))
	param$mean <- array(0, c(g,p))
	param$prop <- rep(0, g)

	for (k in 1:g)
	{
	        #index des elements de classe k
		indk <- which(zapp==k)
                #calculer le proportion
		prop <- length(indk) / n
		#calculer les vecteurs de moyenne de classe k
		moyenne <- Xapp[indk,]
		moyenne <- apply(moyenne, 2, mean)
		#calculer la matrice de covariance
		covariance <- t(as.matrix(Xapp[indk,] - moyenne)) %*% (as.matrix(Xapp[indk,] - moyenne))
		covariance <- covariance / length(indk)
		covariance_corr <- covariance*(length(indk)/(length(indk) - 1))
		MCov <- MCov + covariance_corr*(length(indk) - 1)
		param$mean[k,] <- moyenne
		param$prop[k] <- prop
	}
	MCov <- MCov/(n - g)
	for (k in 1:g)
	{
		param$MCov[,,k] <- MCov
	}

	param
}

nba.app <- function(Xapp, zapp)
{
	n <- dim(Xapp)[1]
	p <- dim(Xapp)[2]
	g <- max(unique(zapp))

	param <- NULL
	param$MCov <- array(0, c(p,p,g))
	param$mean <- array(0, c(g,p))
	param$prop <- rep(0, g)

	for (k in 1:g)
	{
	        #index des elements de classe k
	        indk <- which(zapp==k)
	        #calculer le proportion
	        prop <- length(indk) / n
	        #calculer les vecteurs de moyenne de classe k
	        moyenne <- Xapp[indk,]
	        moyenne <- apply(moyenne, 2, mean)
	        #Calculer la matrice de covariance
	        covariance <- t(as.matrix(Xapp[indk,] - moyenne)) %*% (as.matrix(Xapp[indk,] - moyenne))
	        covariance <- covariance / length(indk)
		param$MCov[,,k] <- diag(diag(covariance))
		param$mean[k,] <- moyenne
		param$prop[k] <- prop
	}

	param
}

ad.val <- function(param, Xtst)
{
        # P(w_k|x )=f_k (x)/f(x) * pi_k
	n <- dim(Xtst)[1]
	p <- dim(Xtst)[2]
	g <- length(param$prop)

	out <- NULL

	prob <- matrix(0, nrow=n, ncol=g)
        f_x <- 0
        f_k <- 0
	for (k in 1:g)
	{
	        f_k <- mvdnorm(Xtst, param$mean[k,], param$MCov[,,k])
		prob[,k] <- f_k * param$prop[k]
		f_x <- f_x + param$prop[k] * f_k
	}
	prob <- prob / f_x
	pred <- max.col(prob)

	out$prob <- prob
	out$pred <- pred

	out
}
