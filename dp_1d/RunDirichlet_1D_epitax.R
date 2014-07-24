datpath = "/lustre/scratch110/sanger/sd11/epitax/dirichlet_input/"
setwd("~/repo/dirichlet/dp_1d")
source("LoadData.R")
source("RunDirichlet_1D.R")
setwd(datpath)

args=commandArgs(TRUE)
run = as.integer(args[1])
no.iters = as.integer(args[2]) # 1300
burn.in.fraction = as.double(args[3])

no.iters.burn.in = floor(no.iters*burn.in.fraction) # 300
set.seed(123)

samplenames = c("PD7404a","PD7405a","PD7407a","PD7408a","PD7409a","PD7410a","PD7414a","PD7416a","PD7418a","PD7422a","PD7426a","PD7427a","PD7428a","PD7431a","PD7433a")
cellularities = c(0.35224,0.8586,0.27594,0.34848,0.26,0.28252,0.24216,0.987,0.44247,0.1528,0.5566,0.191,0.5782,0.34918,0.81844)

samplename = samplenames[run]
cellularity = cellularities[run]

res = load.data(datpath,samplename,list(c("")), cellularity=cellularity, Chromosome="chr", WT.count="WT.count", mut.count="mut.count", subclonal.CN="subclonal.CN", no.chrs.bearing.mut="no.chrs.bearing.mut", data_file_suffix="_allDirichletProcessInfo.txt")
RunDirichlet_1D(mutCount=res$mutCount, WTCount=res$WTCount, no.iters=no.iters, no.iters.burn.in=no.iters.burn.in, cellularity=cellularity, 
                totalCopyNumber=res$totalCopyNumber, mutation.copy.number=res$mutation.copy.number, copyNumberAdjustment=res$copyNumberAdjustment, 
                samplename=samplename, outdir=outdir)
