#setwd("/lustre/scratch110/sanger/dw9/TreeBasedDirichletProcess/Release")
setwd("/lustre/scratch110/sanger/sd11/dirichlet/dp_tree_based/")
source("RunTreeBasedDP.R")

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

#bin.size = 0.025
#increased for ES8
#bin.size = 0.05
#bin.size = NA

resort.mutations = T

#burn.in.fraction = 0.2
#TEST
#no.iters=20
#no.iters=125
#no.iters=1250
#no.iters=10000

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

data=list()
for(s in 1:length(subsamples[[run]])){
	#data[[s]] = read.table(paste("Data/",samplename,subsamples[[run]][s],"_muts_withAllSubclonalCNinfo_June2013.txt",sep=""),header=T)
	data[[s]] = read.table(paste("Data/",samplename,subsamples[[run]][s],"_muts_withAllSubclonalCNinfoAndWithoutMutsOnDeletedChrs_June2013.txt",sep=""),header=T)
}

WTCount = array(0,c(nrow(data[[1]]),no.subsamples))
mutCount = array(0,c(nrow(data[[1]]),no.subsamples))
totalCopyNumber = array(0,dim(WTCount))
copyNumberAdjustment = array(0,dim(WTCount))
non.deleted.muts = vector(mode="logical",length=nrow(data[[1]]))
for(s in 1:length(subsamples[[run]])){
	WTCount[,s] = as.numeric(data[[s]]$WTCount)
	mutCount[,s] = as.numeric(data[[s]]$mutCount)
	totalCopyNumber[,s] = as.numeric(data[[s]]$subclonal.CN)
	copyNumberAdjustment[,s] = as.numeric(data[[s]]$no.chrs.bearing.mut)
	non.deleted.muts[data[[s]]$no.chrs.bearing.mut>0]=T
}


non.zero = which(rowSums(mutCount)>0 & !is.na(rowSums(totalCopyNumber)) & non.deleted.muts)
mutCount = mutCount[non.zero,]
WTCount = WTCount[non.zero,]
totalCopyNumber = totalCopyNumber[non.zero,]
copyNumberAdjustment = copyNumberAdjustment[non.zero,]
for(i in 1:length(subsamples[[run]])){
	data[[i]] = data[[i]][non.zero,]
}

kappa = array(1,dim(mutCount))
for(i in 1:length(subsamples[[run]])){
	#multiply by no.chrs.bearing.mut, so that kappa is the fraction of reads required for fully clonal mutations, rather than mutation at MCN = 1
  kappa[,i] = mutationCopyNumberToMutationBurden(1,data[[i]]$subclonal.CN,cellularity[i]) * data[[i]]$no.chrs.bearing.mut
}

start_time = Sys.time()
if(is.na(bin.size)){
	outdir = paste(samplename,"_treeBasedDirichletProcessOutputs_noIters",no.iters,sep="")
	RunTreeBasedDP(mutCount,WTCount,kappa = kappa, samplename = samplename, subsamplenames = subsamples[[run]], no.iters=no.iters,no.iters.burn.in=no.iters.burn.in,bin.size = bin.size, resort.mutations = resort.mutations, outdir = outdir, parallel=parallel, phase=phase)
  
  }else{
	outdir = paste(samplename,"_treeBasedDirichletProcessOutputs_noIters",no.iters,"_binSize",bin.size,sep="")
	RunTreeBasedDP(mutCount,WTCount,kappa = kappa, samplename = samplename, subsamplenames = subsamples[[run]], no.iters=no.iters,no.iters.burn.in=no.iters.burn.in,bin.size = bin.size, resort.mutations = resort.mutations, outdir = outdir, parallel=parallel, phase=phase)
}
end_time = Sys.time()
# working dir has changed, therefore write this file directly to current dir
write.table(data.frame(diff=c(difftime(end_time, start_time, units='sec')), unit=c('seconds')), file='runtime.txt', quote=F, row.names=F)
print(end_time-start_time)

print(warnings())

q(save="no")
