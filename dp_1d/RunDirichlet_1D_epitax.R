outdir = getwd()
# datpath = "/lustre/scratch110/sanger/sd11/epitax/dirichlet_input/"
setwd("~/repo/dirichlet/dp_nd")
source("LoadData.R")
# source("RunDirichlet_1D.R")
source("interconvertMutationBurdens.R")
source("subclone_Dirichlet_Gibbs_sampler_nD_binomial.R")
# source("OneDimensionalClustering.R")
setwd(outdir)

# args=commandArgs(TRUE)
# run = as.integer(args[1])
# no.iters = as.integer(args[2]) # 1300
# burn.in.fraction = as.double(args[3])

# no.iters.burn.in = floor(no.iters*burn.in.fraction) # 300
# set.seed(123)
# 
# samplenames = c("PD7404a","PD7405a","PD7407a","PD7408a","PD7409a","PD7410a","PD7414a","PD7416a","PD7418a","PD7422a","PD7426a","PD7427a","PD7428a","PD7431a","PD7433a")
# cellularities = c(0.35224,0.8586,0.27594,0.34848,0.26,0.28252,0.24216,0.987,0.44247,0.1528,0.5566,0.191,0.5782,0.34918,0.81844)

# samplename = samplenames[run]
# cellularity = cellularities[run]


getSimpleSubclones <- function(n, size, prob) {
  mutCount = c()
  WTCount = c()
  no.chrs.bearing.mut = c()
  totalCopyNumber = c()
  for (i in 1:length(n)) {
    new.muts = rbinom(n[i], size[i], prob[i])
    mutCount = c(mutCount,new.muts)
    WTCount = c(WTCount, rep(size[i],n[i])-new.muts)
    no.chrs.bearing.mut = c(no.chrs.bearing.mut, rep(1,n[i]))
    totalCopyNumber = c(totalCopyNumber, rep(2,n[i]))
  }
  mutation.copy.number = mutCount/(mutCount+WTCount) * totalCopyNumber #no.chrs.bearing.mut#
  return(list(mutCount=matrix(mutCount), WTCount=matrix(WTCount), copyNumberAdjustment=matrix(no.chrs.bearing.mut), totalCopyNumber=matrix(totalCopyNumber), mutation.copy.number=matrix(mutation.copy.number)))
}

samplename = "simulated"
cellularity = 1
no.iters = 1500
no.iters.burn.in = 300

res = getSimpleSubclones(c(25,30,25), c(100,100,100), c(0.1,0.5,0.9))

C=30
y=res$mutCount
N=res$mutCount+res$WTCount
iter=no.iters
cellularity=cellularity
totalCopyNumber=res$totalCopyNumber
no.chrs.bearing.mut=res$copyNumberAdjustment
normalCopyNumber=array(2,length(y))
iter = no.iters

C=30
mutCount=res$mutCount
WTCount=res$WTCount
iter=no.iters
cellularity=cellularity
totalCopyNumber=res$totalCopyNumber
copyNumberAdjustment=res$copyNumberAdjustment
normalCopyNumber=matrix(array(2,length(mutCount)))
iter = no.iters


# RunDirichlet_1D(mutCount=res$mutCount, WTCount=res$WTCount, no.iters=no.iters, no.iters.burn.in=no.iters.burn.in, cellularity=cellularity, totalCopyNumber=res$totalCopyNumber, 
#                 mutation.copy.number=res$mutation.copy.number, copyNumberAdjustment=res$copyNumberAdjustment, samplename=samplename, outdir=outdir)

# GS.data<-subclone.dirichlet.gibbs(y=res$mutCount,
#                                   N=res$mutCount+res$WTCount, 
#                                   iter=no.iters, 
#                                   cellularity=cellularity, 
#                                   totalCopyNumber=res$totalCopyNumber, 
#                                   no.chrs.bearing.mut=res$copyNumberAdjustment)
GS.data<-subclone.dirichlet.gibbs(mutCount=res$mutCount,
                                  WTCount=res$WTCount, 
                                  iter=no.iters, 
                                  cellularity=cellularity, 
                                  totalCopyNumber=res$totalCopyNumber, 
                                  copyNumberAdjustment=res$copyNumberAdjustment)
density = Gibbs.subclone.density.est.1d(GS.data, 
                                     paste(samplename,"_DirichletProcessplot.png", sep=''), 
                                     post.burn.in.start=no.iters.burn.in, 
                                     post.burn.in.stop=no.iters, # TODO: no.iters or no.iters-no.iters.burn.in ??
                                     y.max=50, 
                                     mutationCopyNumber=res$mutation.copy.number, 
                                     no.chrs.bearing.mut=res$copyNumberAdjustment)

write.csv(GS.data$S.i,paste(outdir,"/",samplename,"_iters",no.iters,"_states.csv",sep=""))
write.csv(GS.data$V.h,paste(outdir,"/",samplename,"_iters",no.iters,"_stickbreaking_weights.csv",sep=""))
write.csv(GS.data$pi.h,paste(outdir,"/",samplename,"_iters",no.iters,"_discreteMutationCopyNumbers.csv",sep=""))
write.csv(GS.data$alpha,paste(outdir,"/",samplename,"_iters",no.iters,"_alpha.csv",sep=""))

subclonal.fraction = res$mutation.copy.number / res$copyNumberAdjustment
subclonal.fraction[is.nan(subclonal.fraction)] = 0
oneDimensionalClustering(samplename, subclonal.fraction, GS.data, density, no.iters, no.iters.burn.in)



# res = load.data(datpath,samplename,list(c("")), cellularity=cellularity, Chromosome="chr", WT.count="WT.count", mut.count="mut.count", subclonal.CN="subclonal.CN", no.chrs.bearing.mut="no.chrs.bearing.mut", data_file_suffix="_allDirichletProcessInfo.txt")
# RunDirichlet_1D(mutCount=res$mutCount, WTCount=res$WTCount, no.iters=no.iters, no.iters.burn.in=no.iters.burn.in, cellularity=cellularity, 
#                 totalCopyNumber=res$totalCopyNumber, mutation.copy.number=res$mutation.copy.number, copyNumberAdjustment=res$copyNumberAdjustment, 
#                 samplename=samplename, outdir=outdir)



