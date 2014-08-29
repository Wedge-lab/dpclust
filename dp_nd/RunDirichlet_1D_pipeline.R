outdir = getwd()
args=commandArgs(TRUE)
run = as.integer(args[1]) # Sample id
no.iters = as.integer(args[2]) # Number of iters
no.iters.burn.in = as.double(args[3]) # Number of iters used for burnin
datpath = toString(args[4]) # Where are input files stored
purity_file = toString(args[5]) # A file containing samplenames and purity

conc_param=0.01
cluster_conc = 5

setwd("~/repo/dirichlet/dp_nd")
source("LoadData.R")
source("RunDirichletProcess.R")
setwd(outdir)

set.seed(123)

sample2purity = read.table(purity_file, header=T, stringsAsFactors=F)
samplename = sample2purity[run,1]
cellularity = sample2purity[run,2]

print(samplename)

res = load.data(datpath,samplename,list(c("")), cellularity=cellularity, Chromosome="chr", WT.count="WT.count", mut.count="mut.count", subclonal.CN="subclonal.CN", no.chrs.bearing.mut="no.chrs.bearing.mut", data_file_suffix="_allDirichletProcessInfo.txt")
RunDirichletProcess(mutCount=matrix(res$mutCount), WTCount=matrix(res$WTCount), no.iters=no.iters, no.iters.burn.in=no.iters.burn.in, cellularity=cellularity, 
                totalCopyNumber=matrix(res$totalCopyNumber), mutation.copy.number=matrix(res$mutation.copy.number), copyNumberAdjustment=matrix(res$copyNumberAdjustment), 
                samplename=samplename, output_folder=outdir, conc_param=conc_param, cluster_conc=cluster_conc)
