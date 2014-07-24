source("interconvertMutationBurdens.R")
source("subclone_Dirichlet_Gibbs_sampler_binomial.R")
source("OneDimensionalClustering.R")

RunDirichlet_1D <- function(mutCount, WTCount, no.iters, no.iters.burn.in, cellularity, totalCopyNumber, copyNumberAdjustment, samplename, outdir) {
  GS.data<-subclone.dirichlet.gibbs(y=mutCount,N=mutCount+WTCount, iter=no.iters, cellularity=cellularity, totalCopyNumber=totalCopyNumber, no.chrs.bearing.mut=copyNumberAdjustment)
  density = Gibbs.subclone.density.est(GS.data, paste(samplename,"_DirichletProcessplot.png", sep=''), post.burn.in.start=no.iters.burn.in, post.burn.in.stop=no.iters, y.max=50, mutationCopyNumber=totalCopyNumber, no.chrs.bearing.mut=copyNumberAdjustment)
  
  write.csv(GS.data$S.i,paste(outdir,"/",samplename,"_iters",no.iters,"_states.csv",sep=""))
  write.csv(GS.data$V.h,paste(outdir,"/",samplename,"_iters",no.iters,"_stickbreaking_weights.csv",sep=""))
  write.csv(GS.data$pi.h,paste(outdir,"/",samplename,"_iters",no.iters,"_discreteMutationCopyNumbers.csv",sep=""))
  write.csv(GS.data$alpha,paste(outdir,"/",samplename,"_iters",no.iters,"_alpha.csv",sep=""))
  
  subclonal.fraction = totalCopyNumber / copyNumberAdjustment
  subclonal.fraction[is.nan(subclonal.fraction)] = 0
  oneDimensionalClustering(samplename, subclonal.fraction, GS.data, density, no.iters, no.iters.burn.in)
}
