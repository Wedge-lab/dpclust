#
# DPClust core algorithm
#

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
  A = 1
  B = conc_param
  
  # Set up data formats for recording iterations
  pi.h = array(NA, c(iter, C,num.timepoints))
  mutBurdens = array(NA, c(C,num.timepoints,num.muts))
  
  V.h = matrix(1, nrow=iter, ncol=C)
  S.i = matrix(NA, nrow=iter, ncol=num.muts)
  Pr.S = matrix(NA, nrow=num.muts, ncol=C)
  alpha = rep(NA, iter)
  
  lower = array(NA,num.timepoints)
  upper = array(NA,num.timepoints)
  
  mutCopyNum = array(NA,c(num.muts,num.timepoints))
  for (t in 1:num.timepoints) {  
    mutCopyNum[,t] = mutationBurdenToMutationCopyNumber(mutCount[,t]/(mutCount[,t]+WTCount[,t]),totalCopyNumber[,t] ,cellularity[t],normalCopyNumber[,t]) / copyNumberAdjustment[,t]

    # DCW adjust to cope with countsPerCopyNum==0,
    # either because WTCount+mutCount==0 or copyNumberAdjustment==0
    mutCopyNum[copyNumberAdjustment[,t]==0,t] = 0

    lower[t] = min(mutCopyNum[,t])
    upper[t] = max(mutCopyNum[,t])
    difference = upper[t]-lower[t]
    lower[t] = lower[t]-difference/10
    upper[t] = upper[t]+difference/10
    # Randomise starting positions of clusters
    pi.h[1,,t]=runif(C,lower[t],upper[t])
    for(c in 1:C){
      mutBurdens[c,t,] = mutationCopyNumberToMutationBurden(pi.h[1,c,t] * copyNumberAdjustment[,t], totalCopyNumber[,t], cellularity[t],normalCopyNumber[,t]) 
      mutBurdens[c,t,mutBurdens[c,t,]==0] = 0.000001
    }
  }	
  V.h[1,] = c(rep(0.5,C-1), 1)
  S.i[1,] = c(1, rep(0,num.muts-1))
  alpha[1] = 1
  V.h[1:iter, C] = rep(1, iter)
  
  for (m in 2:iter) {
    if (m %% 100 == 0){ print(paste(m, " / ", iter, sep=" ")) }
    
    # Update cluster allocation for each individual mutation
    for (k in 1:num.muts) {			
      # use log-space to avoid problems with very high counts
      Pr.S[k,1] <- log(V.h[m-1,1])
      Pr.S[k,2:C] = sapply(2:C, FUN=function(j, V) { log(V[j]) + sum(log(1-V[1:(j-1)])) }, V=V.h[m-1,])
      
      for(t in 1:num.timepoints){
        for(c in 1:C){
          Pr.S[k,c] <- Pr.S[k,c] + mutCount[k,t]*log(mutBurdens[c,t,k]) + WTCount[k,t]*log(1-mutBurdens[c,t,k])
        }
        # It would be faster to use apply here, but it is returning values with small difference as compared to the above code. The scaling below 
        # blows these differences up to significant impact.
        #         Pr.S2[k,] <- Pr.S2[k,] + sapply(1:C, FUN=function(c, t, k, mutCount, mutBurdens, WTCount) { mutCount[k,t]*log(mutBurdens[c,t,k]) + WTCount[k,t]*log(1-mutBurdens[c,t,k]) }, t=t,k=k,mutCount=mutCount,mutBurdens=mutBurdens,WTCount=WTCount)
      }			
      Pr.S[k,is.na(Pr.S[k,])] = 0
      Pr.S[k,] = Pr.S[k,] - max(Pr.S[k,])
      Pr.S[k,] = exp(Pr.S[k,])
      Pr.S[k,] = Pr.S[k,] / sum(Pr.S[k,])				
    }
    
    # Determine which cluster each mutation is assigned to
    S.i[m,] = apply(Pr.S, 1, function(mut) { which(rmultinom(1,1,mut)==1) })
    
    # Update stick-breaking weights
    V.h[m,1:(C-1)] = sapply(1:(C-1), function(S, curr.m, curr.alpha, h) {rbeta(1, 1+sum(S[curr.m,] == h), curr.alpha+sum(S[curr.m,] > h))}, S=S.i, curr.m=m, curr.alpha=alpha[m-1])
    V.h[m,c(V.h[m,1:(C-1)] == 1,FALSE)] = 0.999 # Need to prevent one stick from taking all the remaining weight
    
    # Get expected number of mutant reads per mutation copy number
    countsPerCopyNum = array(NA,c(num.timepoints,num.muts))
    for (t in 1:num.timepoints) {
      countsPerCopyNum[t,] = (mutCount[,t]+WTCount[,t])*mutationCopyNumberToMutationBurden(copyNumberAdjustment[,t],totalCopyNumber[,t],cellularity[t],normalCopyNumber[,t])
    }
    
    # randomise unused pi.h
    for (t in 1:num.timepoints) {
      pi.h[m,,t] = runif(C,lower[t],upper[t])
    }

    for(c in unique(S.i[m,])){
      for(t in 1:num.timepoints){
        #040213 - problem if sum(countsPerCopyNum[t,S.i[m,]==c])==0 fixed
        if(sum(countsPerCopyNum[t,S.i[m,]==c])==0){
          pi.h[m,c,t] = 0.000001
        }else{
          pi.h[m,c,t] = rgamma(1,shape=sum(mutCount[S.i[m,]==c,t]),rate=sum(countsPerCopyNum[t,S.i[m,]==c]))
	  if (is.nan(pi.h[m,c,t])) {
	    print(paste("shape",sum(mutCount[S.i[m,]==c,t]),"rate",sum(countsPerCopyNum[t,S.i[m,]==c]),sep=","))
	    stop("pi.h is NaN")
	  }
        }
      }
    }		
    
    for(t in 1:num.timepoints){
      for(c in 1:C){
        mutBurdens[c,t,]=mutationCopyNumberToMutationBurden(pi.h[m,c,t] * copyNumberAdjustment[,t],totalCopyNumber[,t],cellularity[t],normalCopyNumber[,t])
      }
    }
    
    # Update alpha
    alpha[m] <- rgamma(1, shape=C+A-1, rate=B-sum(log(1-V.h[m,1:(C-1)]))) 
  }
  return(list(S.i=S.i, V.h=V.h, pi.h=pi.h, mutBurdens=mutBurdens, alpha=alpha, y1=mutCount, N1=mutCount+WTCount))
}
