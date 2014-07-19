library(MASS)
library(MCMCpack)
library(mvtnorm)
source('interconvertMutationBurdens.R')
source('subclone_Dirichlet_Gibbs_sampler_nD_binomial.R')

RunDirichletProcess <- function(mutCount, WTCount, totalCopyNumber, copyNumberAdjustment, cellularity, output_folder, noiters, burn.in.iters, conc_param, cluster_conc) {
  
  if(!file.exists(output_folder)){
    dir.create(output_folder)
  }
  
  GS.data<-subclone.dirichlet.gibbs(mutCount=mutCount,WTCount=WTCount,totalCopyNumber=totalCopyNumber, copyNumberAdjustment=copyNumberAdjustment, cellularity=cellularity,iter=noiters,conc_param=conc_param,cluster_conc=cluster_conc)
  
  finalStates=GS.data$S.i[noiters,]
  finalMus=GS.data$pi.h[noiters,,]
  finalFittedMus=finalMus[finalStates,]
  
  write.csv(GS.data$S.i,paste(output_folder,"/",samplename,"_2D_iters",noiters,"_concParam",conc_param,"_clusterWidth",1/cluster_conc,"_states.csv",sep=""))
  write.csv(GS.data$V.h,paste(output_folder,"/",samplename,"_2D_iters",noiters,"_concParam",conc_param,"_clusterWidth",1/cluster_conc,"_stickbreaking_weights.csv",sep=""))
  write.csv(GS.data$pi.h,paste(output_folder,"/",samplename,"_2D_iters",noiters,"_concParam",conc_param,"_clusterWidth",1/cluster_conc,"_discreteMutationCopyNumbers.csv",sep=""))
  write.csv(GS.data$alpha,paste(output_folder,"/",samplename,"_2D_iters",noiters,"_concParam",conc_param,"_clusterWidth",1/cluster_conc,"_alpha.csv",sep=""))
  
  for(i in 1:(length(subsamples[[run]])-1)){
    for(j in (i+1):length(subsamples[[run]])){
      imageFile = paste(output_folder,"/",samplenames[run],subsamples[[run]][i],subsamples[[run]][j],"_iters",noiters,"_concParam",conc_param,"_clusterWidth",1/cluster_conc,"_2D_binomial.png",sep="")
      Gibbs.subclone.density.est(mutation.copy.number[,c(i,j)]/copyNumberAdjustment[,c(i,j)],GS.data,imageFile, post.burn.in.start = burn.in.iters, post.burn.in.stop = noiters, samplenames = paste(samplename,subsamples[[run]][c(i,j)],sep=""),indices=c(i,j))		
    }
  }

}