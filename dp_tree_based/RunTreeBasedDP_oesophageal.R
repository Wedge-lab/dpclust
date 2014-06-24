datpath = "/lustre/scratch110/sanger/sd11/dirichlet/dp_nd/"
setwd("/lustre/scratch110/sanger/sd11/dirichlet/dp_tree_based/")
source("RunTreeBasedDP.R")

args=commandArgs(TRUE)
run = as.integer(args[1])
bin.size = as.double(args[2])
no.iters = as.integer(args[3])
burn.in.fraction = as.double(args[4])

parallel = TRUE
resort.mutations = TRUE
no.iters.burn.in = floor(no.iters*burn.in.fraction)
set.seed(456)

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

data=list()
for(s in 1:length(subsamples[[run]])){
  data[[s]] = read.table(paste(datpath,samplename,subsamples[[run]][s],"_allDirichletProcessInfo.txt",sep=""),header=T,sep="\t")
  data[[s]] = data[[s]][data[[s]]$Chromosome %in% 1:22,]
}

# TODO: Add in code that makes sure all mutations are mentioned in all files? See nD Run script.

WTCount = array(0,c(nrow(data[[1]]),no.subsamples))
mutCount = array(0,c(nrow(data[[1]]),no.subsamples))
totalCopyNumber = array(0,dim(WTCount))
copyNumberAdjustment = array(0,dim(WTCount))
non.deleted.muts = vector(mode="logical",length=nrow(data[[1]]))
for(s in 1:length(subsamples[[run]])){
  WTCount[,s] = as.numeric(data[[s]]$WT.count)
  mutCount[,s] = as.numeric(data[[s]]$mut.count)
  totalCopyNumber[,s] = as.numeric(data[[s]]$subclonal.CN)
  copyNumberAdjustment[,s] = as.numeric(data[[s]]$no.chrs.bearing.mut)
  non.deleted.muts[data[[s]]$no.chrs.bearing.mut>0]=T
}

kappa = array(1,dim(mutCount))
for(i in 1:length(subsamples[[run]])){
  #multiply by no.chrs.bearing.mut, so that kappa is the fraction of reads required for fully clonal mutations, rather than mutation at MCN = 1
  kappa[,i] = mutationCopyNumberToMutationBurden(1,data[[i]]$subclonal.CN,cellularity[i]) * data[[i]]$no.chrs.bearing.mut
}

# Remove those mutations that have missing values
no.muts = nrow(totalCopyNumber)
not.there.wt = sapply(1:no.muts, FUN=function(i, dat1) { sum(is.na(dat1[i,]))>0}, dat1=WTCount)
not.there.mut = sapply(1:no.muts, FUN=function(i, dat1) { sum(is.na(dat1[i,]))>0}, dat1=mutCount)
not.there.cn = sapply(1:no.muts, FUN=function(i, dat1) { sum(is.na(dat1[i,]))>0}, dat1=totalCopyNumber)
not.there.cna = sapply(1:no.muts, FUN=function(i, dat1) { sum(is.na(dat1[i,]))>0}, dat1=copyNumberAdjustment)
not.there.kappa = sapply(1:no.muts, FUN=function(i, dat1) { sum(is.na(dat1[i,]))>0}, dat1=kappa)
select = !(not.there.wt | not.there.mut | not.there.cn | not.there.cna | not.there.kappa)
WTCount = WTCount[select,]
mutCount = mutCount[select,]
totalCopyNumber = totalCopyNumber[select,]
copyNumberAdjustment = copyNumberAdjustment[select,]
non.deleted.muts = non.deleted.muts[select]
kappa = kappa[select,]
print(paste("Removed ",no.muts-nrow(WTCount)))

start_time = Sys.time()
if(is.na(bin.size)){
  outdir = paste(samplename,"_treeBasedDirichletProcessOutputs_noIters",no.iters,sep="")
  RunTreeBasedDP(mutCount,WTCount,kappa = kappa, samplename = samplename, subsamplenames = subsamples[[run]], no.iters=no.iters,no.iters.burn.in=no.iters.burn.in,bin.size = bin.size, resort.mutations = resort.mutations, outdir = outdir, parallel=parallel)
  
}else{
  outdir = paste(samplename,"_treeBasedDirichletProcessOutputs_noIters",no.iters,"_binSize",bin.size,sep="")
  RunTreeBasedDP(mutCount,WTCount,kappa = kappa, samplename = samplename, subsamplenames = subsamples[[run]], no.iters=no.iters,no.iters.burn.in=no.iters.burn.in,bin.size = bin.size, resort.mutations = resort.mutations, outdir = outdir, parallel=parallel)
}
end_time = Sys.time()
# working dir has changed, therefore write this file directly to current dir
write.table(data.frame(diff=c(difftime(end_time, start_time, units='sec')), unit=c('seconds')), file='runtime.txt', quote=F, row.names=F)
print(end_time-start_time)

print(warnings())

q(save="no")
