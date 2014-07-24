outdir = getwd()
args=commandArgs(TRUE)
run = as.integer(args[1]) # Sample id
no.iters = as.integer(args[2]) # Number of iters
no.iters.burn.in = as.double(args[3]) # Number of iters used for burnin
datpath = toString(args[4]) # Where are input files stored
purity_file = toString(args[5]) # A file containing samplenames and purity

setwd("~/repo/dirichlet/dp_1d")
source("LoadData.R")
source("RunDirichlet_1D.R")
setwd(outdir)

set.seed(123)

sample2purity = read.table(purity_file, header=T)
samplename = sample2purity[run,1]
cellularity = sample2purity[run,2]

print(samplename)

res = load.data(datpath,samplename,list(c("")), cellularity=cellularity, Chromosome="chr", WT.count="WT.count", mut.count="mut.count", subclonal.CN="subclonal.CN", no.chrs.bearing.mut="no.chrs.bearing.mut", data_file_suffix="_allDirichletProcessInfo.txt")
RunDirichlet_1D(mutCount=res$mutCount, WTCount=res$WTCount, no.iters=no.iters, no.iters.burn.in=no.iters.burn.in, cellularity=cellularity, 
                totalCopyNumber=res$totalCopyNumber, mutation.copy.number=res$mutation.copy.number, copyNumberAdjustment=res$copyNumberAdjustment, 
                samplename=samplename, outdir=outdir)
