subclone.dirichlet.gibbs <- function(C=30, y, N, totalCopyNumber=array(1,length(y)),cellularity=1, normalCopyNumber=array(2,length(y)) , no.chrs.bearing.mut = array(1,length(y)),iter=1000) {
	# y is a vector of the number of reads reporting each variant
	# N is a vector of the number of reads in total across the base in question (in the same order as Y obviously!)
	# C is the maximum number of clusters in the Dirichlet process
	# iter is the number of iterations of the Gibbs sampler
	
	num.muts <- length(y)
	print(paste("num.muts=",num.muts,sep=""))
	
	# Hyperparameters for alpha
	A <- B <- 0.01
	
	# Set up data formats for recording iterations
	pi.h <- matrix(NA, nrow=iter, ncol=C)
	V.h <- matrix(1, nrow=iter, ncol=C)
	S.i <- matrix(NA, nrow=iter, ncol=num.muts)
	Pr.S <- matrix(NA, nrow=num.muts, ncol=C)
	alpha <- rep(NA, iter)	
	mutBurdens <- array(NA, c(iter, C,num.muts))

	mutCopyNum = mutationBurdenToMutationCopyNumber(y/N,totalCopyNumber,cellularity,normalCopyNumber) / no.chrs.bearing.mut
	lower=min(mutCopyNum)
	upper=max(mutCopyNum)
	difference=upper-lower
	lower=lower-difference/10
	upper=upper+difference/10
	# randomise starting positions of clusters
	pi.h[1,]=runif(C,lower,upper)
	for(c in 1:C){
		mutBurdens[1,c,]=mutationCopyNumberToMutationBurden(pi.h[1,c]*no.chrs.bearing.mut,totalCopyNumber,cellularity,normalCopyNumber)
	}	
	
	V.h[1,] <- c(rep(0.5,C-1), 1)
	S.i[1,] <- c(1, rep(0,num.muts-1))
	alpha[1] <- 1
	V.h[1:iter, C] <- rep(1, iter)
	
	for (m in 2:iter) {
		if (m / 100 == round(m/100)) {print(m)}
		
		# Update cluster allocation for each individual mutation
		for (k in 1:num.muts) {			
			#use log-space to avoid problems with very high counts
			Pr.S[k,1] <- log(V.h[m-1,1]) + y[k]*log(mutBurdens[m-1,1,k]) + (N[k]-y[k])*log(1-mutBurdens[m-1,1,k])
			Pr.S[k,2:C] <- sapply(2:C, function(V, obs.y, obs.N, pi, curr.k, j) {log(V[j]) + sum(log(1-V[1:(j-1)])) + obs.y[curr.k]*log(pi[j]) + (obs.N[curr.k] - obs.y[curr.k])*log(1-pi[j])}, V=V.h[m-1,], pi=mutBurdens[m-1,,k], obs.y=y, obs.N = N, curr.k=k)
			
			if(sum(is.na(Pr.S[k,]))>0){
				print("err1")
				print(Pr.S[k,])
				print(V.h[m-1,])
				print(pi.h[m-1,])
			}
			Pr.S[k,] = Pr.S[k,] - max(Pr.S[k,],na.rm=T)
			Pr.S[k,]=exp(Pr.S[k,])
			if(sum(is.na(Pr.S[k,]))>0){
				print("err2")
				print(Pr.S[k,])
			}
			Pr.S[k,] <- Pr.S[k,] / sum(Pr.S[k,],na.rm=T)
			if(sum(is.na(Pr.S[k,]))>0){
				print("err3")
				print(Pr.S[k,])
			}			
			Pr.S[k,is.na(Pr.S[k,])] = 0
		}

		S.i[m,] <- sapply(1:num.muts, function(Pr, k) {sum(rmultinom(1,1,Pr[k,]) * (1:length(Pr[k,])))}, Pr=Pr.S)
		
		# Update stick-breaking weights
		V.h[m,1:(C-1)]  <- sapply(1:(C-1), function(S, curr.m, curr.alpha, h) {rbeta(1, 1+sum(S[curr.m,] == h), curr.alpha+sum(S[curr.m,] > h))}, S=S.i, curr.m=m, curr.alpha=alpha[m-1])
		V.h[m,c(V.h[m,1:(C-1)] == 1,FALSE)] <- 0.999 # Need to prevent one stick from taking all the remaining weight
		
		# Update location of subclones / clusters
		countsPerCopyNum=N*mutationCopyNumberToMutationBurden(1,totalCopyNumber,cellularity,normalCopyNumber)*no.chrs.bearing.mut
		
		#190512 randomise unused pi.h
		pi.h[m,]=runif(C,lower,upper)
		
		mutBurdens[m,,]=mutBurdens[m-1,,]
		for(c in unique(S.i[m,])){
			pi.h[m,c] = rgamma(1,shape=sum(y[S.i[m,]==c]),rate=sum(countsPerCopyNum[S.i[m,]==c]))
			mutBurdens[m,c,]=mutationCopyNumberToMutationBurden(pi.h[m,c]*no.chrs.bearing.mut,totalCopyNumber,cellularity,normalCopyNumber)
		}
		
		# Update alpha
		alpha[m] <- rgamma(1, shape=C+A-1, rate=B-sum(log(1-V.h[m,1:(C-1)]))) 
	}
	return(list(S.i=S.i, V.h=V.h, pi.h=pi.h, alpha=alpha, y1=y, N1=N))
} 

Gibbs.subclone.density.est <- function(GS.data, pngFile, density.smooth = 0.1, post.burn.in.start = 3000, post.burn.in.stop = 10000, density.from = 0, y.max=5, mutationCopyNumber = NA,no.chrs.bearing.mut = NA) {
	print(paste("density.smooth=",density.smooth,sep=""))
	png(filename=pngFile,,width=1500,height=1000)
	# GS.data is the list output from the above function
	# density.smooth is the smoothing factor used in R's density() function
	# post.burn.in.start is the number of iterations to drop from the Gibbs sampler output to allow the estimates to equilibrate on the posterior
	
	xlabel = "mutation copy number"
	if(is.na(mutationCopyNumber)){
		print("No mutationCopyNumber. Using mutation burden")
		y <- GS.data$y1
		N <- GS.data$N1
		mutationCopyNumber = y/N
		xlabel = "mutation burden"
	}

	if(!is.na(no.chrs.bearing.mut)){
		mutationCopyNumber = mutationCopyNumber / no.chrs.bearing.mut
		xlabel = "fraction of tumour cells"
	}

	V.h.cols <- GS.data$V.h
	pi.h.cols <- GS.data$pi.h
	wts <- matrix(NA, nrow=dim(V.h.cols)[1], ncol=dim(V.h.cols)[2])
	wts[,1] <- V.h.cols[,1]
	wts[,2] <- V.h.cols[,2] * (1-V.h.cols[,1])
	for (i in 3:dim(wts)[2]) {wts[,i] <- apply(1-V.h.cols[,1:(i-1)], MARGIN=1, FUN=prod) * V.h.cols[,i]}

	post.ints <- matrix(NA, ncol=post.burn.in.stop - post.burn.in.start + 1, nrow=512)

	x.max = ceiling(max(mutationCopyNumber)*12)/10
	
	xx <- density(c(pi.h.cols[post.burn.in.start-1,]), weights=c(wts[post.burn.in.start,]) / sum(c(wts[post.burn.in.start,])), adjust=density.smooth, from=density.from, to=x.max)$x


	for (i in post.burn.in.start : post.burn.in.stop) {
		post.ints[,i - post.burn.in.start + 1] <- density(c(pi.h.cols[i-1,]), weights=c(wts[i,]) / sum(c(wts[i,])), adjust=density.smooth, from=density.from, to=x.max)$y
	}

	polygon.data = c(apply(post.ints, MARGIN=1, FUN=quantile, probs=0.975), rev(apply(post.ints, MARGIN=1, FUN=quantile, probs=0.025)))
	if(is.na(y.max)){
		y.max=ceiling(max(polygon.data)/10)*10
	}
	
	par(mar = c(5,6,4,1)+0.1)
	hist(mutationCopyNumber, breaks=seq(-0.1, x.max, 0.025), col="lightgrey",freq=FALSE, xlab=xlabel,main="", ylim=c(0,y.max),cex.axis=2,cex.lab=2)
	polygon(c(xx, rev(xx)), polygon.data, border="plum4", col=cm.colors(1,alpha=0.3))

	yy = apply(post.ints, MARGIN=1, FUN=quantile, probs=0.5)

	lines(xx, yy, col="plum4", lwd=3)
	
	dev.off()
	print(paste("highest density is at ",xx[which.max(yy)],sep=""))
	write.table(cbind(xx,yy),gsub(".png","density.txt",pngFile),sep="\t",col.names=c(gsub(" ",".",xlabel),"median.density"),row.names=F,quote=F)
	write.table(polygon.data,gsub(".png","polygonData.txt",pngFile),sep="\t",row.names=F,quote=F)
}