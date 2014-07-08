datpath = "/lustre/scratch110/sanger/sd11/dirichlet/dp_tree_based/Data/"
setwd("/lustre/scratch110/sanger/sd11/dirichlet/dp_tree_based/")
source("RunTreeBasedDP.R")
source('LoadData.R')

args=commandArgs(TRUE)
run = as.integer(args[1])
bin.size = as.double(args[2])
no.iters = as.integer(args[3])
burn.in.fraction = as.double(args[4])
if (length(args) < 5) {
  phase = NA
} else {
  phase = args[5]
}

parallel = TRUE
resort.mutations = T
no.iters.burn.in = floor(no.iters*burn.in.fraction)
set.seed(123)


samplenames=c("H3122","C32","ES7","ES8")
subsamples = list()
subsamples[[1]] = c("","_TR_1","_TR_7","_TR_9")
subsamples[[2]] = c("","_AZDR_1","_AZDR_9","_AZDR_11","-PLAZ-1","-PLAZ-2","-PLAZ-3","-PLAZ-4","_PLXR_1")
subsamples[[3]] = c("-OLAR-1-2uM","-OLAR-3-2uM")
subsamples[[4]] = c("-OLAR-2-5uM","-OLAR-3-5uM","-OLAR-4-5uM","-OLAR-5-5uM","-OLAR-6-5uM")
samplename = samplenames[run]

no.subsamples = length(subsamples[[run]])
cellularity = rep(1,no.subsamples)


res = load.data(datpath, samplename, subsamples[[run]], 'chr', 'WTCount', 'mutCount', 'subclonal.CN', 'no.chrs.bearing.mut', "_muts_withAllSubclonalCNinfoAndWithoutMutsOnDeletedChrs_June2013.txt")
WTCount = res$WTCount
mutCount = res$mutCount
totalCopyNumber = res$totalCopyNumber
copyNumberAdjustment = res$copyNumberAdjustment
non.deleted.muts = res$non.deleted.muts
kappa = res$kappa


if(is.na(bin.size)){
	outdir = paste(samplename,"_treeBasedDirichletProcessOutputs_noIters",no.iters,"_burnin",no.iters.burn.in,sep="")
	RunTreeBasedDP(mutCount,WTCount,kappa = kappa, samplename = samplename, subsamplenames = subsamples[[run]], no.iters=no.iters,no.iters.burn.in=no.iters.burn.in,bin.size = bin.size, resort.mutations = resort.mutations, outdir = outdir, parallel=parallel, phase=phase)
  
  }else{
	outdir = paste(samplename,"_treeBasedDirichletProcessOutputs_noIters",no.iters,"_binSize",bin.size,"_burnin",no.iters.burn.in,sep="")
	RunTreeBasedDP(mutCount,WTCount,kappa = kappa, samplename = samplename, subsamplenames = subsamples[[run]], no.iters=no.iters,no.iters.burn.in=no.iters.burn.in,bin.size = bin.size, resort.mutations = resort.mutations, outdir = outdir, parallel=parallel, phase=phase)
}

print(warnings())

q(save="no")
