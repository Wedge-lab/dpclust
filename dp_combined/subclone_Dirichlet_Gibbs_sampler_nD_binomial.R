library(lattice)
library(KernSmooth)
source("interconvertMutationBurdens.R")
source("PlotDensities.R")

subclone.dirichlet.gibbs <- function(mutCount, WTCount, totalCopyNumber=array(1,dim(mutCount)), normalCopyNumber=array(2,dim(mutCount)), copyNumberAdjustment = array(1,dim(mutCount)),C=30, cellularity=rep(1,ncol(mutCount)),iter=1000,conc_param=1,cluster_conc=10) {
  # y is a p-by-q matrix of the number of reads reporting each variant (p=number of mutations,q=number of timepoints / related samples)
  # N is a p-by-q matrix of the number of reads in total across the base in question (in the same order as Y obviously!, p=number of mutations,q=number of timepoints / related samples)
  # C is the maximum number of clusters in the Dirichlet process
  # iter is the number of iterations of the Gibbs sampler
  if (is.null(copyNumberAdjustment)) { copyNumberAdjustment = array(1,dim(mutCount)) }
  
  num.muts <- NROW(mutCount)
  num.timepoints  = NCOL(mutCount)
  
  # Hyperparameters for alpha
  #A <- B <- 0.01
  #strong prior on alpha
  A=1
  B=conc_param
  
  # Set up data formats for recording iterations
  ## Nr iterations dim not used - sd11
  pi.h <- array(NA, c(iter, C,num.timepoints))
  #mutBurdens <- array(NA, c(iter, C,num.timepoints,num.muts))
  #080214 reduce memory burden. It's not necessary to save all mutBurdens, because they can be derived from pi.h
  mutBurdens <- array(NA, c(C,num.timepoints,num.muts))
  
  V.h <- matrix(1, nrow=iter, ncol=C)
  S.i <- matrix(NA, nrow=iter, ncol=num.muts)
  Pr.S <- matrix(NA, nrow=num.muts, ncol=C)
  alpha <- rep(NA, iter)
  
  # Initialise (now randomised - see below)
  #pi.h[1,,]<-array(rep(seq(0.1, 1, 0.9/(C-1)),num.timepoints),c(C,num.timepoints))
  
  lower = array(NA,num.timepoints)
  upper = array(NA,num.timepoints)
  
  mutCopyNum = array(NA,c(num.muts,num.timepoints))
  for(t in 1:num.timepoints){  
    mutCopyNum[,t] = mutationBurdenToMutationCopyNumber(mutCount[,t]/(mutCount[,t]+WTCount[,t]),totalCopyNumber[,t] ,cellularity[t],normalCopyNumber[,t]) / copyNumberAdjustment[,t]
    lower[t]=min(mutCopyNum[,t])
    upper[t]=max(mutCopyNum[,t])
    difference = upper[t]-lower[t]
    lower[t]=lower[t]-difference/10
    upper[t]=upper[t]+difference/10
    # 190512- randomise starting positions of clusters
    pi.h[1,,t]=runif(C,lower[t],upper[t])
    for(c in 1:C){
      #mutBurdens[1,c,t,]=mutationCopyNumberToMutationBurden(pi.h[1,c,t],totalCopyNumber[,t],cellularity[t])
      #mutBurdens[1,c,t,]=mutationCopyNumberToMutationBurden(pi.h[1,c,t] * copyNumberAdjustment[,t], totalCopyNumber[,t], cellularity[t]) 
      #080214
      mutBurdens[c,t,]=mutationCopyNumberToMutationBurden(pi.h[1,c,t] * copyNumberAdjustment[,t], totalCopyNumber[,t], cellularity[t],normalCopyNumber[,t]) 
    }
    #mutBurdens[1,,t]=mutationCopyNumberToMutationBurden(pi.h[1,,t],totalCopyNumber[t,],cellularity[t])
  }	
  V.h[1,] <- c(rep(0.5,C-1), 1)
  S.i[1,] <- c(1, rep(0,num.muts-1))
  alpha[1] <- 1
  V.h[1:iter, C] <- rep(1, iter)
  
  for (m in 2:iter) {
    if(m %% 100 == 0){print(m)}
    
    # Update cluster allocation for each individual mutation
    for (k in 1:num.muts) {			
      #use log-space to avoid problems with very high counts
      Pr.S[k,1] <- log(V.h[m-1,1])
      Pr.S[k,2:C] = sapply(2:C, FUN=function(j, V) { log(V[j]) + sum(log(1-V[1:(j-1)])) }, V=V.h[m-1,])
      
      
      for(t in 1:num.timepoints){
        for(c in 1:C){
          #Pr.S[k,c] <- Pr.S[k,c] + mutCount[k,t]*log(mutBurdens[m-1,c,t,k]) + WTCount[k,t]*log(1-mutBurdens[m-1,c,t,k])
          #080214
          Pr.S[k,c] <- Pr.S[k,c] + mutCount[k,t]*log(mutBurdens[c,t,k]) + WTCount[k,t]*log(1-mutBurdens[c,t,k])
        }
        # It would be faster to use apply here, but it is returning values with small difference as compared to the above code. The scaling below 
        # blows these differences up to significant impact.
        #         Pr.S2[k,] <- Pr.S2[k,] + sapply(1:C, FUN=function(c, t, k, mutCount, mutBurdens, WTCount) { mutCount[k,t]*log(mutBurdens[c,t,k]) + WTCount[k,t]*log(1-mutBurdens[c,t,k]) }, t=t,k=k,mutCount=mutCount,mutBurdens=mutBurdens,WTCount=WTCount)
        
#         if(sum(is.na(Pr.S[k,]))>0){
#           print("err1")
#           print(Pr.S[k,])
#           print(V.h[m-1,])
#           print(pi.h[m-1,])
#         }
#         Pr.S[k,] = Pr.S[k,] - max(Pr.S[k,],na.rm=T)
#         Pr.S[k,]=exp(Pr.S[k,])
#         if(sum(is.na(Pr.S[k,]))>0){
#           print("err2")
#           print(Pr.S[k,])
#         }
#         Pr.S[k,] <- Pr.S[k,] / sum(Pr.S[k,],na.rm=T)
#         if(sum(is.na(Pr.S[k,]))>0){
#           print("err3")
#           print(Pr.S[k,])
#         }			
#         Pr.S[k,is.na(Pr.S[k,])] = 0
        
      }			
      Pr.S[k,is.na(Pr.S[k,])] = 0
      Pr.S[k,] = Pr.S[k,] - max(Pr.S[k,])
      Pr.S[k,]=exp(Pr.S[k,])
      Pr.S[k,] <- Pr.S[k,] / sum(Pr.S[k,])				
    }
    
    if(sum(is.na(Pr.S))>0){
      print(paste("Pr.S=",Pr.S))
    }
    
    ## Determine which cluster each mutation is assigned to
    # 		S.i[m,] <- sapply(1:num.muts, function(Pr, k) {sum(rmultinom(1,1,Pr[k,]) * (1:length(Pr[k,])))}, Pr=Pr.S)
    S.i[m,] = apply(Pr.S, 1, function(mut) { which(rmultinom(1,1,mut)==1) })
    
    # Update stick-breaking weights
    V.h[m,1:(C-1)]  <- sapply(1:(C-1), function(S, curr.m, curr.alpha, h) {rbeta(1, 1+sum(S[curr.m,] == h), curr.alpha+sum(S[curr.m,] > h))}, S=S.i, curr.m=m, curr.alpha=alpha[m-1])
    V.h[m,c(V.h[m,1:(C-1)] == 1,FALSE)] <- 0.999 # Need to prevent one stick from taking all the remaining weight
    
    #200512 - get expected number of mutant reads per mutation copy number
    countsPerCopyNum=array(NA,c(num.timepoints,num.muts))
    for(t in 1:num.timepoints){
      countsPerCopyNum[t,]=(mutCount[,t]+WTCount[,t])*mutationCopyNumberToMutationBurden(copyNumberAdjustment[,t],totalCopyNumber[,t],cellularity[t],normalCopyNumber[,t])
    }
    #pi.h[m,,]=pi.h[m-1,,]
    #190512 randomise unused pi.h
    for(t in 1:num.timepoints){
      pi.h[m,,t]=runif(C,lower[t],upper[t])
    }
    #080214 - no longer needed, because mutBurdens are not saved for every iter
    #mutBurdens[m,,,]=mutBurdens[m-1,,,]
    for(c in unique(S.i[m,])){
      for(t in 1:num.timepoints){
        #040213 - problem if sum(countsPerCopyNum[t,S.i[m,]==c])==0 fixed
        if(sum(countsPerCopyNum[t,S.i[m,]==c])==0){
          pi.h[m,c,t] = 0
        }else{
          pi.h[m,c,t] = rgamma(1,shape=sum(mutCount[S.i[m,]==c,t]),rate=sum(countsPerCopyNum[t,S.i[m,]==c]))
          #if(is.nan(pi.h[m,c,t])){
          #	stop()
          #}
        }
      }
    }		
    for(t in 1:num.timepoints){
      for(c in 1:C){
        #mutBurdens[m,c,t,]=mutationCopyNumberToMutationBurden(pi.h[m,c,t] * copyNumberAdjustment[,t],totalCopyNumber[,t],cellularity[t])
        #080214
        mutBurdens[c,t,]=mutationCopyNumberToMutationBurden(pi.h[m,c,t] * copyNumberAdjustment[,t],totalCopyNumber[,t],cellularity[t],normalCopyNumber[,t])
      }
    }
    
    if(sum(is.na(pi.h[m,,]))){
      print(paste("pi.h=",pi.h[m,,]))
      print(paste("m=",m,sep=""))
    }
    
    # Update alpha
    alpha[m] <- rgamma(1, shape=C+A-1, rate=B-sum(log(1-V.h[m,1:(C-1)]))) 
  }
  return(list(S.i=S.i, V.h=V.h, pi.h=pi.h, mutBurdens=mutBurdens, alpha=alpha, y1=mutCount, N1=mutCount+WTCount))
}


Gibbs.subclone.density.est <- function(burden, GS.data, pngFile, density.smooth = 0.1, post.burn.in.start = 3000, post.burn.in.stop = 10000, samplenames = c("sample 1","sample 2"), indices = NA) {
  print(paste("density.smooth=",density.smooth,sep=""))
  
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
  
  #print(paste("wts=",wts[1000,],sep=""))
  #2-D kernel smoother
  #library(ks)
  
  gridsize=c(64L,64L)
  
  num.timepoints = NCOL(burden)
  if(num.timepoints==2){
    range=list(c(floor(min(burden[,1])*10)-1,ceiling(max(burden[,1])*10)+1)/10,c(floor(min(burden[,2])*10)-1,ceiling(max(burden[,2])*10)+1)/10)
    print(range)
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
    print(paste("gridsize=",gridsize,sep=""))
    
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
  #for (i in post.burn.in.start : post.burn.in.stop) {
  for(i in 1:length(sampledIters)){
    density.data=array(NA,c(0,num.timepoints))
    for(j in 1:no.clusters){
      #density.data=rbind(density.data,t(array(pi.h.cols[i,j,],c(num.timepoints,round(no.density.points*wts[i,j])))))
      
      #density.data=rbind(density.data,t(array(pi.h.cols[sampledIters[i],j,],c(num.timepoints,round(no.density.points*wts[sampledIters[i],j])))))
      #180512 use pi.h from previous generation
      density.data=rbind(density.data,t(array(pi.h.cols[sampledIters[i]-1,j,],c(num.timepoints,round(no.density.points*wts[sampledIters[i],j])))))
    }
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
         pngFile=pngFile, 
         samplename_x=samplenames[1], 
         samplename_y=samplenames[2], 
         max.plotted.value=NA)
  
  # Save the data of the plot to replot lateron for figure fine tuning
  write.csv(xvals,gsub(".png","_xvals.csv",pngFile))
  write.csv(yvals,gsub(".png","_yvals.csv",pngFile))	
  write.csv(median.density,gsub(".png","_zvals.csv",pngFile))
  
  return(list(fraction.of.tumour.cells=burden, median.density=median.density, xvals=xvals, yvals=yvals))
}

# Gibbs.subclone.density.est.1d <- function(GS.data, pngFile, samplename, density.smooth=0.1, post.burn.in.start=3000, post.burn.in.stop=10000, density.from=0, x.max=2, y.max=5, mutationCopyNumber=NULL, no.chrs.bearing.mut=NULL) {
#   print(paste("density.smooth=",density.smooth,sep=""))
#   png(filename=pngFile,,width=1500,height=1000)
#   # GS.data is the list output from the above function
#   # density.smooth is the smoothing factor used in R's density() function
#   # post.burn.in.start is the number of iterations to drop from the Gibbs sampler output to allow the estimates to equilibrate on the posterior
#   
#   if(is.null(mutationCopyNumber)){
#     print("No mutationCopyNumber. Using mutation burden")
#     y <- GS.data$y1
#     N <- GS.data$N1
#     mutationCopyNumber = y/N
#   }
#   
#   if(!is.null(no.chrs.bearing.mut)){
#     mutationCopyNumber = mutationCopyNumber / no.chrs.bearing.mut
#   }
#   
#   V.h.cols <- GS.data$V.h # weights
#   pi.h.cols <- GS.data$pi.h[,,1] # discreteMutationCopyNumbers
#   wts <- matrix(NA, nrow=dim(V.h.cols)[1], ncol=dim(V.h.cols)[2])
#   wts[,1] <- V.h.cols[,1]
#   wts[,2] <- V.h.cols[,2] * (1-V.h.cols[,1])
#   for (i in 3:dim(wts)[2]) {wts[,i] <- apply(1-V.h.cols[,1:(i-1)], MARGIN=1, FUN=prod) * V.h.cols[,i]}
#   
#   post.ints <- matrix(NA, ncol=post.burn.in.stop - post.burn.in.start + 1, nrow=512)
#   
#   if (is.na(x.max)) {
#     x.max = ceiling(max(mutationCopyNumber)*12)/10
#   }
#   
#   xx <- density(c(pi.h.cols[post.burn.in.start-1,]), weights=c(wts[post.burn.in.start,]) / sum(c(wts[post.burn.in.start,])), adjust=density.smooth, from=density.from, to=x.max)$x
#  
#   for (i in post.burn.in.start : post.burn.in.stop) {
#     post.ints[,i - post.burn.in.start + 1] <- density(c(pi.h.cols[i-1,]), weights=c(wts[i,]) / sum(c(wts[i,])), adjust=density.smooth, from=density.from, to=x.max)$y
#   }
#   
#   polygon.data = c(apply(post.ints, MARGIN=1, FUN=quantile, probs=0.975), rev(apply(post.ints, MARGIN=1, FUN=quantile, probs=0.025)))
#   if(is.na(y.max)){
#     y.max=ceiling(max(polygon.data)/10)*10
#   }
#   
#   yy = apply(post.ints, MARGIN=1, FUN=quantile, probs=0.5)
#   density = cbind(xx, yy)
#   
#   plot1D(density, 
#          polygon.data, 
#          pngFile=pngFile, 
#          density.from=0, 
#          y.max=y.max, 
#          x.max=x.max, 
#          mutationCopyNumber=mutationCopyNumber, 
#          no.chrs.bearing.mut=no.chrs.bearing.mut,
#          samplename=samplename)
#   
#   
#   print(paste("highest density is at ",xx[which.max(yy)],sep=""))
#   write.table(density,gsub(".png","density.txt",pngFile),sep="\t",col.names=c(gsub(" ",".",xlabel),"median.density"),row.names=F,quote=F)
#   write.table(polygon.data,gsub(".png","polygonData.txt",pngFile),sep="\t",row.names=F,quote=F)
#   
#   return(data.frame(fraction.of.tumour.cells=xx, median.density=yy))
# }

Gibbs.subclone.density.est.1d <- function(GS.data, pngFile, samplename, density.smooth=0.1, post.burn.in.start=3000, post.burn.in.stop=10000, density.from=0, x.max=2, y.max=5, mutationCopyNumber=NULL, no.chrs.bearing.mut=NULL) {
  print(paste("density.smooth=",density.smooth,sep=""))
  png(filename=pngFile,,width=1500,height=1000)
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
  
  print(head(polygon.data))
  
  plot1D(density, 
         polygon.data, 
         pngFile=pngFile, 
         density.from=0, 
         y.max=y.max, 
         x.max=x.max, 
         mutationCopyNumber=mutationCopyNumber, 
         no.chrs.bearing.mut=no.chrs.bearing.mut,
         samplename=samplename)
  
  
  print(paste("highest density is at ",xx[which.max(yy)],sep=""))
  write.table(density,gsub(".png","density.txt",pngFile),sep="\t",col.names=c(gsub(" ",".",xlabel),"median.density"),row.names=F,quote=F)
  write.table(polygon.data,gsub(".png","polygonData.txt",pngFile),sep="\t",row.names=F,quote=F)
  
  return(list(density=data.frame(fraction.of.tumour.cells=xx, median.density=yy), polygon.data=polygon.data))
}

Gibbs.subclone.density.est.nd <- function(burden, GS.data, density.smooth = 0.1, post.burn.in.start = 200, post.burn.in.stop = 1000, max.burden=NA, resolution=NA) {
  #
  # Use this when assigning mutations to clusters by using multiDimensional clustering. 
  # It returns density values across the grid and a confidence interval that are used
  # by the multi D clustering for determining where a mutation should be assigned.
  #
  print(paste("density.smooth=",density.smooth,sep=""))
  
  V.h.cols <- GS.data$V.h
  if("pi.h" %in% names(GS.data)){
    pi.h.cols <- GS.data$pi.h
  }else{
    pi.h.cols <- GS.data$theta$mu
  }
  C = dim(pi.h.cols)[2]
  
  wts <- matrix(NA, nrow=dim(V.h.cols)[1], ncol=dim(V.h.cols)[2])
  wts[,1] <- V.h.cols[,1]
  wts[,2] <- V.h.cols[,2] * (1-V.h.cols[,1])
  for (i in 3:dim(wts)[2]) {wts[,i] <- apply(1-V.h.cols[,1:(i-1)], MARGIN=1, FUN=prod) * V.h.cols[,i]}
  
  #print(paste("wts=",wts[1000,],sep=""))
  #n-D kernel smoother
  library(ks)
  #library(KernSmooth)
  
  num.timepoints = NCOL(burden)
  print(paste("num.timepoints=",num.timepoints,sep=""))
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
    print("creating evaluation.points")
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
    print(sampledIters[i])
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
  #save(post.ints,file="DEBUG_post.ints.RData")
  
  
  print("calculating median and 95% CI")
  median.density = array(NA,gridsize)
  lower.CI = array(NA,gridsize)
  
  print("arrays created")
  
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
  for(i in 1:prod(gridsize)){
    indices = (i-1) %% gridsize[1] + 1
    for(j in 2:num.timepoints){
      new.index = (i-1) %/% prod(gridsize[1:(j-1)]) + 1
      new.index = (new.index - 1) %% gridsize[j] + 1
      indices = c(indices,new.index)
    }
    #print(indices)
    #print(dim(post.ints))
    #print(cbind(array(rep(indices,each=no.samples),c(no.samples,num.timepoints)),sampledIters))
    
    median.density[array(indices,c(1,num.timepoints))] = median(post.ints[cbind(array(rep(indices,each=no.samples),c(no.samples,num.timepoints)),1:no.samples)])
    lower.CI[array(indices,c(1,num.timepoints))] = quantile(post.ints[cbind(array(rep(indices,each=no.samples),c(no.samples,num.timepoints)),1:no.samples)],probs=0.05)
    if(i %% 1000 ==0){
      print(paste(i," of ",prod(gridsize),sep=""))
    }
  }
  print("finished calculating median and 95% CI")
  
  return(list(range=range,gridsize=gridsize,median.density=median.density, lower.CI=lower.CI))    
}