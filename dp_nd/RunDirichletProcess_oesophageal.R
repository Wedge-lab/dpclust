datpath = "/lustre/scratch110/sanger/sd11/dirichlet/dp_nd/"
wd = getwd()
setwd('~/repo/dirichlet/dp_nd')
source('LoadData.R')
source('RunDirichletProcess.R')
setwd(wd)

args=commandArgs(TRUE)
run = as.integer(args[1])
no.iters = as.integer(args[2])
burn.in.fraction = as.double(args[3])

conc_param=1
cluster_conc = 5
#no.iters=500
#no.iters=1000
#test
#noiters=5
# no.iters=20
no.iters.burn.in = floor(no.iters*burn.in.fraction)

samplenames = c("Case1","Case2","Case3","Case4")
subsamples = list(
  c("Primary1","Primary2","Metastasis1","Metastasis2"),
  c("Primary1","Primary2","Primary3","Metastasis1","Metastasis2"),
  c("Primary1","Metastasis1"),
  c("Primary1","Metastasis1")
)
cellularities = list(
  c(0.5857628,0.4327561,0.4124731,0.5564892),
  c(0.6810764,0.2532192,0.1817953,0.06062584,0.102386),
  c(0.7727005,0.6570972),
  c(0.8886442,0.7912862)
)

samplename = samplenames[run]
no.subsamples = length(subsamples[[run]])
cellularity = cellularities[[run]]
outdir = paste(samplename,"_nD_DirichletProcessOutputs_noIters",no.iters,"_burnin",no.iters.burn.in,sep="")

res = load.data(datpath, samplename, subsamples[[run]], cellularity, 'Chromosome', 'WT.count', 'mut.count', 'subclonal.CN', 'no.chrs.bearing.mut', "_allDirichletProcessInfo.txt")

RunDirichletProcess(mutCount=res$mutCount, WTCount=res$WTCount, totalCopyNumber=res$totalCopyNumber, copyNumberAdjustment=res$copyNumberAdjustment, mutation.copy.number=res$mutation.copy.number, cellularity=cellularity, output_folder=outdir, noiters=no.iters, burn.in.iters=no.iters.burn.in, subsamplesrun=subsamples[[run]], samplename=samplename, conc_param=conc_param, cluster_conc=cluster_conc)

