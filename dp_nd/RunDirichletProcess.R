library(MASS)
library(MCMCpack)
library(mvtnorm)
source('interconvertMutationBurdens.R')
source('subclone_Dirichlet_Gibbs_sampler_nD_binomial.R')
source("OneDimensionalClustering.R")

RunDirichletProcess <- function(mutCount, WTCount, totalCopyNumber, copyNumberAdjustment, mutation.copy.number, cellularity, output_folder, noiters, burn.in.iters, subsamplesrun, samplename, conc_param, cluster_conc) {
#   print(head(mutation.copy.number))
#   print(sum(mutation.copy.number==0))
#   print(head(copyNumberAdjustment))  
#   print(sum(copyNumberAdjustment==0))
  if(!file.exists(output_folder)){
    dir.create(output_folder)
  }

  GS.data<-subclone.dirichlet.gibbs(mutCount=mutCount,WTCount=WTCount,totalCopyNumber=totalCopyNumber, copyNumberAdjustment=copyNumberAdjustment, cellularity=cellularity,iter=noiters,conc_param=conc_param,cluster_conc=cluster_conc)
#   
#   finalStates=GS.data$S.i[noiters,]
#   finalMus=GS.data$pi.h[noiters,,]
#   finalFittedMus=finalMus[finalStates,]
#   
  write.csv(GS.data$S.i,paste(output_folder,"/",samplename,"_2D_iters",noiters,"_concParam",conc_param,"_clusterWidth",1/cluster_conc,"_states.csv",sep=""))
  write.csv(GS.data$V.h,paste(output_folder,"/",samplename,"_2D_iters",noiters,"_concParam",conc_param,"_clusterWidth",1/cluster_conc,"_stickbreaking_weights.csv",sep=""))
  write.csv(GS.data$pi.h,paste(output_folder,"/",samplename,"_2D_iters",noiters,"_concParam",conc_param,"_clusterWidth",1/cluster_conc,"_discreteMutationCopyNumbers.csv",sep=""))
  write.csv(GS.data$alpha,paste(output_folder,"/",samplename,"_2D_iters",noiters,"_concParam",conc_param,"_clusterWidth",1/cluster_conc,"_alpha.csv",sep=""))
  
  
#   GS.data = list()
#   GS.data$S.i = read.csv(paste(output_folder,"/",samplename,"_2D_iters",noiters,"_concParam",conc_param,"_clusterWidth",1/cluster_conc,"_states.csv",sep=""),row.names=1)
#   GS.data$V.h = read.csv(paste(output_folder,"/",samplename,"_2D_iters",noiters,"_concParam",conc_param,"_clusterWidth",1/cluster_conc,"_stickbreaking_weights.csv",sep=""),row.names=1)
#   GS.data$pi.h = read.csv(paste(output_folder,"/",samplename,"_2D_iters",noiters,"_concParam",conc_param,"_clusterWidth",1/cluster_conc,"_discreteMutationCopyNumbers.csv",sep=""),row.names=1)
#   GS.data$alpha = read.csv(paste(output_folder,"/",samplename,"_2D_iters",noiters,"_concParam",conc_param,"_clusterWidth",1/cluster_conc,"_alpha.csv",sep=""),row.names=1)
  # nD dataset, plot sample versus sample
  if (ncol(mutCount) > 1) {
    for(i in 1:(length(subsamplesrun)-1)){
      for(j in (i+1):length(subsamplesrun)){
        imageFile = paste(output_folder,"/",samplename,subsamplesrun[i],subsamplesrun[j],"_iters",noiters,"_concParam",conc_param,"_clusterWidth",1/cluster_conc,"_2D_binomial.png",sep="")
        Gibbs.subclone.density.est(mutation.copy.number[,c(i,j)]/copyNumberAdjustment[,c(i,j)],GS.data,imageFile, post.burn.in.start = burn.in.iters, post.burn.in.stop = noiters, samplenames = paste(samplename,subsamplesrun[c(i,j)],sep=""),indices=c(i,j))		
      }
    }
  # 1D dataset, plot just the single density
  } else {
    density = Gibbs.subclone.density.est.1d(GS.data, 
                                            paste(samplename,"_DirichletProcessplot.png", sep=''), 
                                            post.burn.in.start=no.iters.burn.in, 
                                            post.burn.in.stop=no.iters, # TODO: no.iters or no.iters-no.iters.burn.in ??
                                            y.max=50, 
                                            mutationCopyNumber=mutation.copy.number, 
                                            no.chrs.bearing.mut=copyNumberAdjustment)
    subclonal.fraction = mutation.copy.number / copyNumberAdjustment
    subclonal.fraction[is.nan(subclonal.fraction)] = 0
    oneDimensionalClustering(samplename, subclonal.fraction, GS.data, density, no.iters, no.iters.burn.in)
  }

}
