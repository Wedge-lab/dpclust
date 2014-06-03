DensityEstimator <- function(clustered.thetas, thetas, density.smooth = 0.1, density.from = 0, x.max=NA, y.max=NA, main="") {
	# density.smooth is the smoothing factor used in R's density() function
		
	no.iters = ncol(clustered.thetas)
	post.ints <- matrix(NA, ncol=no.iters, nrow=512)	

	if(is.na(x.max)){
		x.max = ceiling(max(thetas,na.rm=T)*12)/10
	}
	if(nrow(clustered.thetas)==1){
		clustered.thetas = rbind(clustered.thetas,clustered.thetas)
	}
	xx <- density(rep(clustered.thetas[,1],2), adjust=density.smooth, from=density.from, to=x.max)$x

	for (i in 1:no.iters) {
		post.ints[,i] <- density(clustered.thetas[,i], adjust=density.smooth, from=density.from, to=x.max)$y
	}

	polygon.data = c(apply(post.ints, MARGIN=1, FUN=quantile, probs=0.975), rev(apply(post.ints, MARGIN=1, FUN=quantile, probs=0.025)))
	if(is.na(y.max)){
		y.max=ceiling(max(polygon.data)/10)*10
	}
	
	hist(thetas[thetas<=x.max], breaks=seq(-0.1, x.max, 0.025), col="lightgrey",freq=FALSE, xlab="subclonal fraction",main=main, ylim=c(0,y.max))

	polygon(c(xx, rev(xx)), polygon.data, border="plum4", col=cm.colors(1,alpha=0.3))
	yy = apply(post.ints, MARGIN=1, FUN=quantile, probs=0.5)
	lines(xx, yy, col="plum4", lwd=3)
	
	return(xx[which.max(yy)])
}