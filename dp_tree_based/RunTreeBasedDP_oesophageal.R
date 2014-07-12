datpath = "/lustre/scratch110/sanger/sd11/dirichlet/dp_nd/"
setwd("~/repo/dirichlet/dp_tree_based/")
source("RunTreeBasedDP.R")
source("LoadData.R")
setwd("/lustre/scratch110/sanger/sd11/dirichlet/dp_tree_based/")

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
resort.mutations = TRUE
no.iters.burn.in = floor(no.iters*burn.in.fraction)
set.seed(123)

samplenames=c("Case1","Case2","Case3","Case4","Case1","Case1_1050","Case1_3300","Case1_1050","Case1_3300")
subsamples = list()
subsamples[[1]] = c("Primary1","Primary2","Metastasis1","Metastasis2")
subsamples[[2]] = c("Primary1","Primary2","Primary3","Metastasis1","Metastasis2")
subsamples[[3]] = c("Primary1","Metastasis1")
subsamples[[4]] = c("Primary1","Metastasis1")
subsamples[[5]] = c("Primary1","Metastasis1")
subsamples[[6]] = c("_Primary1","_Primary2","_Metastasis1","_Metastasis2")
subsamples[[7]] = c("_Primary1","_Primary2","_Metastasis1","_Metastasis2")
subsamples[[8]] = c("_Primary1","_Metastasis1")
subsamples[[9]] = c("_Primary1","_Metastasis1")

cellularities = list(
  c(0.5857628,0.4327561,0.4124731,0.5564892),
  c(0.6810764,0.2532192,0.1817953,0.06062584,0.102386),
  c(0.7727005,0.6570972),
  c(0.8886442,0.7912862),
  c(0.5857628,0.4124731),
  c(0.5857628,0.4327561,0.4124731,0.5564892),
  c(0.5857628,0.4327561,0.4124731,0.5564892),
  c(0.5857628,0.4124731),
  c(0.5857628,0.4124731)
)
samplename = samplenames[run]
no.subsamples = length(subsamples[[run]])
cellularity = cellularities[[run]]


res = load.data(datpath, samplename, subsamples[[run]], 'Chromosome', 'WT.count', 'mut.count', 'subclonal.CN', 'no.chrs.bearing.mut', "_allDirichletProcessInfo.txt")
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
