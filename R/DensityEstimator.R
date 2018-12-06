#
# This file contains various functions to estimate densities
#


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

	post.ints = sapply(1:no.iters,FUN=function(i,clustered.thetas1,adjust,from,to) { density(clustered.thetas[,i], adjust=adjust, from=from, to=to)$y },clustered.thetas1=clustered.thetas,adjust=density.smooth,from=density.from,to=x.max)

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

Gibbs.subclone.density.est.1d <- function(GS.data, pngFile, samplename, density.smooth=0.1, post.burn.in.start=3000, post.burn.in.stop=10000, density.from=0, x.max=2, y.max=5, mutationCopyNumber=NULL, no.chrs.bearing.mut=NULL) {
  
  # GS.data is the list output from the above function
  # density.smooth is the smoothing factor used in R's density() function
  # post.burn.in.start is the number of iterations to drop from the Gibbs sampler output to allow the estimates to equilibrate on the posterior
  xlabel = "mutation_copy_number" 
  if(is.null(mutationCopyNumber)){
    print("No mutationCopyNumber. Using mutation burden")
    y <- GS.data$y1
    N <- GS.data$N1
    mutationCopyNumber = y/N
    xlabel = "mutation_burden"
  }
  
  # Save the original MCN to give to the plotter below. If MCN and chroms.bearing.muts are both submitted MCN was
  # scaled twice (once just below and once in the plotter) giving rise to wrong figures.
  mutationCopyNumber.original = mutationCopyNumber
  
  if(!is.null(no.chrs.bearing.mut)){
    mutationCopyNumber = mutationCopyNumber / no.chrs.bearing.mut
    xlabel = "fraction_of_cells"
  }
  
  V.h.cols <- GS.data$V.h # weights
  pi.h.cols <- GS.data$pi.h[,,1] # discreteMutationCopyNumbers
  wts <- matrix(NA, nrow=dim(V.h.cols)[1], ncol=dim(V.h.cols)[2])
  wts[,1] <- V.h.cols[,1]
  wts[,2] <- V.h.cols[,2] * (1-V.h.cols[,1])
  for (i in 3:dim(wts)[2]) {wts[,i] <- apply(1-V.h.cols[,1:(i-1)], MARGIN=1, FUN=prod) * V.h.cols[,i]}
  
  post.ints <- matrix(NA, ncol=post.burn.in.stop - post.burn.in.start + 1, nrow=512)
  
  if (is.na(x.max)) {
    x.max = ceiling(max(mutationCopyNumber)*12)/10
  }
  
  xx <- density(c(pi.h.cols[post.burn.in.start-1,]), weights=c(wts[post.burn.in.start,]) / sum(c(wts[post.burn.in.start,])), adjust=density.smooth, from=density.from, to=x.max)$x
  
  for (i in post.burn.in.start : post.burn.in.stop) {
    post.ints[,i - post.burn.in.start + 1] <- density(c(pi.h.cols[i-1,]), weights=c(wts[i,]) / sum(c(wts[i,])), adjust=density.smooth, from=density.from, to=x.max)$y
  }
  
  polygon.data = c(apply(post.ints, MARGIN=1, FUN=quantile, probs=0.975), rev(apply(post.ints, MARGIN=1, FUN=quantile, probs=0.025)))
  if(is.na(y.max)){
    y.max=ceiling(max(polygon.data)/10)*10
  }
  
  yy = apply(post.ints, MARGIN=1, FUN=quantile, probs=0.5)
  density = cbind(xx, yy)
  
  # print(head(polygon.data))
  
  plot1D(density,
         polygon.data,
         pngFile=pngFile,
         density.from=0,
         y.max=y.max,
         x.max=x.max,
         mutationCopyNumber=mutationCopyNumber.original,
         no.chrs.bearing.mut=no.chrs.bearing.mut,
         samplename=samplename)
  
  write.table(density,gsub(".png","density.txt",pngFile),sep="\t",col.names=c(gsub(" ",".",xlabel),"median.density"),row.names=F,quote=F)
  write.table(polygon.data,gsub(".png","polygonData.txt",pngFile),sep="\t",row.names=F,quote=F)
  
  return(list(density=data.frame(fraction.of.tumour.cells=xx, median.density=yy), polygon.data=polygon.data))
}

Gibbs.subclone.density.est <- function(burden, GS.data, pngFile, density.smooth = 0.1, post.burn.in.start = 3000, post.burn.in.stop = 10000, samplenames = c("sample 1","sample 2"), indices = NA) {
  
  V.h.cols <- GS.data$V.h
  if("pi.h" %in% names(GS.data)){
    if(is.na(indices)){
      pi.h.cols <- GS.data$pi.h
    }else{
      pi.h.cols <- GS.data$pi.h[,,indices]
    }
  }else{
    pi.h.cols <- GS.data$theta$mu
  }
  
  wts <- matrix(NA, nrow=dim(V.h.cols)[1], ncol=dim(V.h.cols)[2])
  wts[,1] <- V.h.cols[,1]
  wts[,2] <- V.h.cols[,2] * (1-V.h.cols[,1])
  for (i in 3:dim(wts)[2]) {wts[,i] <- apply(1-V.h.cols[,1:(i-1)], MARGIN=1, FUN=prod) * V.h.cols[,i]}
  
  #2-D kernel smoother
  #library(ks)
  
  gridsize=c(64L,64L)
  num.timepoints = NCOL(burden)
  if(num.timepoints==2){
    range=list(c(floor(min(burden[,1])*10)-1,ceiling(max(burden[,1])*10)+1)/10,c(floor(min(burden[,2])*10)-1,ceiling(max(burden[,2])*10)+1)/10)
    
    #gridsize[1]=round(100*(range[[1]][2]-range[[1]][1]))+1
    #gridsize[2]=round(100*(range[[2]][2]-range[[2]][1]))+1
    if(range[[1]][2]-range[[1]][1]>20){
      gridsize[1]=round(20*(range[[1]][2]-range[[1]][1]))+1			
    }else if(range[[1]][2]-range[[1]][1]>8){
      gridsize[1]=round(50*(range[[1]][2]-range[[1]][1]))+1
    }else if(range[[1]][2]-range[[1]][1]>3){
      gridsize[1]=round(100*(range[[1]][2]-range[[1]][1]))+1
    }else{
      gridsize[1]=round(200*(range[[1]][2]-range[[1]][1]))+1
    }
    if(range[[2]][2]-range[[2]][1]>20){
      gridsize[2]=round(20*(range[[2]][2]-range[[2]][1]))+1
    }else if(range[[2]][2]-range[[2]][1]>8){
      gridsize[2]=round(50*(range[[2]][2]-range[[2]][1]))+1
    }else if(range[[2]][2]-range[[2]][1]>3){
      gridsize[2]=round(100*(range[[2]][2]-range[[2]][1]))+1
    }else{
      gridsize[2]=round(200*(range[[2]][2]-range[[2]][1]))+1
    }
    
    #if added 150512
    sampledIters = post.burn.in.start : post.burn.in.stop
    #don't use the intitial state
    sampledIters = sampledIters[sampledIters!=1]		
    if(length(sampledIters) > 1000){
      post.ints <- array(NA, c(gridsize[1], gridsize[2], 1000))
      sampledIters=floor(post.burn.in.start + (1:1000) * (post.burn.in.stop - post.burn.in.start)/1000)			
    }else{
      post.ints <- array(NA, c(gridsize[1], gridsize[2], length(sampledIters)))
    }
    
    
  }else{
    post.ints <- array(NA, c(512, post.burn.in.stop - post.burn.in.start + 1))
    xx=array(NA,c(1,512))
    range=c(floor(min(burden)*10)-1,ceiling(max(burden[,1])*10)+1)/10
  }
  
  #no.density.points=1000
  no.density.points=10000
  no.clusters=ncol(V.h.cols)
  
  for(i in 1:length(sampledIters)){
    if (i %% 100 ==0) { print(paste(i, "/", length(sampledIters))) }
    
    density.data = lapply(1:no.clusters, function(j) {
      t(array(pi.h.cols[sampledIters[i]-1,j,], c(num.timepoints,round(no.density.points*wts[sampledIters[i],j]))))
    })
    density.data = do.call(rbind, density.data)
    
    if(i==1){
      write.csv(density.data,gsub(".png",paste("_densityData",i,".csv",sep=""),pngFile))
    }
    if(num.timepoints==2){
      d=bkde2D(density.data,bandwidth=density.smooth,gridsize=gridsize,range.x=range)
      #if(i==post.burn.in.start){
      if(i==1){	
        xvals=d$x1
        yvals=d$x2
      }
      #post.ints[,,i - post.burn.in.start + 1]=d$fhat
      post.ints[,,i]=d$fhat	
      
    }else{
      d=bkde(density.data,bandwidth=density.smooth,gridsize=512L,range.x=range)
      #post.ints[,i - post.burn.in.start + 1]=d$y
      post.ints[,i]=d$y
    }
  }
  if(num.timepoints==2){
    median.density=apply(post.ints, MARGIN=c(1,2), FUN=median)
  }else{
    median.density=apply(post.ints, MARGIN=1, FUN=median)
  }
  
  # Create the plots both with and without mutations 
  plotnD(xvals=xvals, 
         yvals=yvals, 
         zvals=median.density, 
         subclonal.fraction_x=burden[,1], 
         subclonal.fraction_y=burden[,2], 
         pngFile=gsub(".png","_withoutMutations.png", pngFile), 
         samplename_x=samplenames[1], 
         samplename_y=samplenames[2], 
         max.plotted.value=NA)
  
  plotnD(xvals=xvals, 
         yvals=yvals, 
         zvals=median.density, 
         subclonal.fraction_x=burden[,1], 
         subclonal.fraction_y=burden[,2], 
         pngFile=gsub(".png","_withMutations.png", pngFile), 
         samplename_x=samplenames[1], 
         samplename_y=samplenames[2], 
         max.plotted.value=NA,
         plot_mutations=T)
  
  # Save the data of the plot to replot lateron for figure fine tuning
  write.csv(xvals,gsub(".png","_xvals.csv",pngFile))
  write.csv(yvals,gsub(".png","_yvals.csv",pngFile))	
  write.csv(median.density,gsub(".png","_zvals.csv",pngFile))
  
  return(list(fraction.of.tumour.cells=burden, median.density=median.density, xvals=xvals, yvals=yvals))
}


Gibbs.subclone.density.est.nd <- function(burden, GS.data, density.smooth = 0.1, post.burn.in.start = 200, post.burn.in.stop = 1000, max.burden=NA, resolution=NA) {
  #
  # Use this when assigning mutations to clusters by using multiDimensional clustering. 
  # It returns density values across the grid and a confidence interval that are used
  # by the multi D clustering for determining where a mutation should be assigned.
  #
  
  V.h.cols <- GS.data$V.h
  if("pi.h" %in% names(GS.data)){
    pi.h.cols <- GS.data$pi.h
  }else{
    pi.h.cols <- GS.data$theta$mu
  }
  C = dim(pi.h.cols)[2]
  
  print("Estimating density for all MCMC iterations...")
  wts <- matrix(NA, nrow=dim(V.h.cols)[1], ncol=dim(V.h.cols)[2])
  wts[,1] <- V.h.cols[,1]
  wts[,2] <- V.h.cols[,2] * (1-V.h.cols[,1])
  for (i in 3:dim(wts)[2]) {wts[,i] <- apply(1-V.h.cols[,1:(i-1)], MARGIN=1, FUN=prod) * V.h.cols[,i]}
  
  num.timepoints = NCOL(burden)
  gridsize=rep(64,num.timepoints)
  range = array(NA,c(num.timepoints,2))
  for(n in 1:num.timepoints){
    range[n,] = c(floor(min(burden[,1])*10)-2,ceiling(max(burden[,1])*10)+2)/10
    #don't calculate density for outliers (should only be applied when burden is subclonal fraction)
    if(!is.na(max.burden) && range[n,2] > max.burden){
      range[n,2] = max.burden
    }
    #avoid trying to create arrays that are too large and disallowed (and will be very slow)
    if(is.na(resolution)){
      if(num.timepoints>=6|num.timepoints==5){
        if(range[n,2]-range[n,1]>3){
          grid.resolution=5
        }else if(range[n,2]-range[n,1]>1.5){
          grid.resolution=10
        }else{
          grid.resolution=20
        }
      }else if(num.timepoints==4){
        if(range[n,2]-range[n,1]>8){
          grid.resolution=5
        }else if(range[n,2]-range[n,1]>3){
          grid.resolution=10
        }else if(range[n,2]-range[n,1]>1.5){
          grid.resolution=20
        }else{
          grid.resolution=40
        }                       
      }else{
        if(range[n,2]-range[n,1]>20){
          grid.resolution=5                       
        }else if(range[n,2]-range[n,1]>8){
          grid.resolution=10
        }else if(range[n,2]-range[n,1]>3){
          grid.resolution=25
        }else{
          grid.resolution=50
        }
      }
    }else{
      grid.resolution = resolution
    }
    
    gridsize[n]=round(grid.resolution*(range[n,2]-range[n,1]))+1
  }
  #if added 150512
  sampledIters = post.burn.in.start : post.burn.in.stop
  #don't use the intitial state
  sampledIters = sampledIters[sampledIters!=1]            
  if(length(sampledIters) > 1000){
    post.ints <- array(NA, c(gridsize, 1000))
    sampledIters=floor(post.burn.in.start + (1:1000) * (post.burn.in.stop - post.burn.in.start)/1000)                       
  }else{
    post.ints <- array(NA, c(gridsize, length(sampledIters)))
  }
  
  no.clusters=ncol(V.h.cols)
  no.eval.points = prod(gridsize)
  
  if(num.timepoints>=4){
    evaluation.points = array(NA,c(no.eval.points,num.timepoints))
    evaluation.points[,1] = rep(seq(range[1,1],range[1,2],(range[1,2]-range[1,1])/(gridsize[1]-1)),times = prod(gridsize[2:num.timepoints]))                
    for(i in 2:num.timepoints){
      if(i==num.timepoints){
        evaluation.points[,i] = rep(seq(range[i,1],range[i,2],(range[i,2]-range[i,1])/(gridsize[i]-1)),each = prod(gridsize[1:(i-1)]))
      }else{
        evaluation.points[,i] = rep(seq(range[i,1],range[i,2],(range[i,2]-range[i,1])/(gridsize[i]-1)),each = prod(gridsize[1:(i-1)]),times = prod(gridsize[(i+1):num.timepoints]))
      }
    }
    #evaluation.points = list()
    #for(i in 1:num.timepoints){
    #       evaluation.points[[i]] = seq(range[i,1],range[i,2],(range[i,2]-range[i,1])/(gridsize[i]-1))
    #}
  }
  
  for(i in 1:length(sampledIters)){
    if (i %% 100 == 0) { print(paste(i, "/", length(sampledIters), sep=" ")) }
    if(num.timepoints>=4){
      #use weights
      d=kde(pi.h.cols[sampledIters[i]-1,,],H=diag(num.timepoints)*density.smooth,eval.points = evaluation.points,w=C*wts[sampledIters[i],])                   
    }else{
      #use weights
      d=kde(pi.h.cols[sampledIters[i]-1,,],H=diag(num.timepoints)*density.smooth,gridsize=gridsize,xmin=range[,1],xmax=range[,2],w=C*wts[sampledIters[i],])
    }
    if(i==1){       
      eval.points = d$eval.points
    }
    post.ints[((i-1)*no.eval.points+1):(i*no.eval.points)]=d$estimate
  }
  
  print("Estimating density confidence intervals...")
  median.density = array(NA,gridsize)
  lower.CI = array(NA,gridsize)
  
  #median.density=apply(post.ints, MARGIN=1:num.timepoints, FUN=median)
  ##lower.CI=apply(post.ints, MARGIN=1:num.timepoints, FUN=quantile, probs = 0.025)
  ##upper.CI=apply(post.ints, MARGIN=1:num.timepoints, FUN=quantile, probs = 0.975)       
  ##return(list(range=range,gridsize=gridsize,median.density=median.density, lower.CI=lower.CI, upper.CI= upper.CI))
  ##just do one-tailed test, we're not interested in the upper threshold
  #lower.CI=apply(post.ints, MARGIN=1:num.timepoints, FUN=quantile, probs = 0.05)
  
  #try to reduce memory requirements
  #median.vector = vector(mode="numeric",length=prod(gridsize))
  #lower.CI.vector = vector(mode="numeric",length=prod(gridsize))
  no.samples = length(sampledIters)
  
  # Determine how often to print a status update
  num_iters = prod(gridsize)
  if (num_iters < 1000) {
    print_status_iters = 100
  } else if (num_iters < 10000) {
    print_status_iters = 1000
  } else if (num_iters < 100000) {
    print_status_iters = 10000
  } else {
    print_status_iters = 100000
  }
  
  for(i in 1:num_iters){
    if(i %% print_status_iters==0){ print(paste(i,"/", num_iters, sep=" ")) }
    indices = (i-1) %% gridsize[1] + 1
    for(j in 2:num.timepoints){
      new.index = (i-1) %/% prod(gridsize[1:(j-1)]) + 1
      new.index = (new.index - 1) %% gridsize[j] + 1
      indices = c(indices,new.index)
    }
    
    median.density[array(indices,c(1,num.timepoints))] = median(post.ints[cbind(array(rep(indices,each=no.samples),c(no.samples,num.timepoints)),1:no.samples)])
    lower.CI[array(indices,c(1,num.timepoints))] = quantile(post.ints[cbind(array(rep(indices,each=no.samples),c(no.samples,num.timepoints)),1:no.samples)],probs=0.05)
  }
  
  return(list(range=range, gridsize=gridsize, median.density=median.density, lower.CI=lower.CI))    
}
