#args=commandArgs(TRUE)
#run = as.integer(args[1])

run = 1

#noiters=500
#noiters=1000
#test
#noiters=5
noiters=20

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
cellularity = cellularities[[run]]

samplename = samplenames[run]
cellularity = cellularities[[run]]

#basedir = "/lustre/scratch109/sanger/dw9/2D_Dirichlet_Process/"
datadir = '/lustre/scratch110/sanger/sd11/dirichlet/dp_nd/'
#setwd(paste(basedir,"oesophageal",sep=""))
setwd(datadir)
codedir = '/nfs/users/nfs_s/sd11/repo/dirichlet/dp_nd/'
library(MASS)
library(MCMCpack)
library(mvtnorm)
source(paste(codedir,"subclone_Dirichlet_Gibbs_sampler_nD_binomial.R",sep=''))
source(paste(codedir,"interconvertMutationBurdens.R",sep=''))

output_folder = "oesophageal_nD_DirichletProcess_adjusted"
if(!file.exists(output_folder)){
	dir.create(output_folder)
}

data=list()
for(s in 1:length(subsamples[[run]])){
	data[[s]] = read.table(paste(samplename,subsamples[[run]][s],"_allDirichletProcessInfo.txt",sep=""),header=T,sep="\t")
	#don't use X chromosome
	data[[s]] = data[[s]][data[[s]]$Chromosome %in% 1:22,]
}
matched.data = list()
for(s in 1:length(subsamples[[run]])){
	matched.data[[s]] = data[[s]][0,]
}
matched.loci = list()
for(chr in 1:22){
	for(s in 1:length(subsamples[[run]])){
		if(s==1){
			matched.loci[[chr]] = data[[s]]$Position[data[[s]]$Chromosome==chr]
		}else{
			matched.loci[[chr]] = matched.loci[[chr]][matched.loci[[chr]] %in% data[[s]]$Position[data[[s]]$Chromosome==chr]]
		}
	}
	for(s in 1:length(subsamples[[run]])){
		matched.data[[s]] = rbind(matched.data[[s]],data[[s]][data[[s]]$Chromosome==chr & data[[s]]$Position %in% matched.loci[[chr]],])
	}
}
data = matched.data
print(paste("dim data[[1]]=",dim(data[[1]]),sep=""))

conc_param=1
cluster_conc = 5

## Store nmuts and nsubsamples in separate variable - sd11
WTCount = array(0,c(nrow(data[[1]]),length(subsamples[[run]])))
mutCount = array(0,c(nrow(data[[1]]),length(subsamples[[run]])))
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
totalCount = WTCount+mutCount
no.reads = sapply(1:nrow(totalCount),function(i){0 %in% totalCount[i,]}) | sapply(1:nrow(copyNumberAdjustment),function(i){0 %in% copyNumberAdjustment[i,]})

non.zero = which(rowSums(mutCount)>0 & !is.na(rowSums(totalCopyNumber)) & non.deleted.muts & !no.reads)
mutCount = mutCount[non.zero,]
WTCount = WTCount[non.zero,]
totalCopyNumber = totalCopyNumber[non.zero,]
copyNumberAdjustment = copyNumberAdjustment[non.zero,]
for(i in 1:length(subsamples[[run]])){
	data[[i]] = data[[i]][non.zero,]
}

## This does not remove anything - sd11
in.seg = which(!is.na(rowSums(totalCopyNumber)))
mutCount = mutCount[in.seg,]
WTCount = WTCount[in.seg,]
totalCopyNumber = totalCopyNumber[in.seg,]
copyNumberAdjustment = copyNumberAdjustment[in.seg,]
for(i in 1:length(subsamples[[run]])){
	data[[i]] = data[[i]][in.seg,]
}

## Can be merged with for loop that creates WTCount - sd11
mutation.copy.number = array(NA,dim(totalCopyNumber))
for(s in 1:length(subsamples[[run]])){
	mutation.copy.number[,s] = as.numeric(data[[s]]$mutation.copy.number)
}


GS.data<-subclone.dirichlet.gibbs(mutCount=mutCount,WTCount=WTCount,totalCopyNumber=totalCopyNumber, copyNumberAdjustment=copyNumberAdjustment, cellularity=cellularity,iter=noiters,conc_param=conc_param,cluster_conc=cluster_conc)

finalStates=GS.data$S.i[noiters,]
finalMus=GS.data$pi.h[noiters,,]
finalFittedMus=finalMus[finalStates,]

write.csv(GS.data$S.i,paste(output_folder,"/",samplename,"_2D_iters",noiters,"_concParam",conc_param,"_clusterWidth",1/cluster_conc,"_states.csv",sep=""))
write.csv(GS.data$V.h,paste(output_folder,"/",samplename,"_2D_iters",noiters,"_concParam",conc_param,"_clusterWidth",1/cluster_conc,"_stickbreaking_weights.csv",sep=""))
write.csv(GS.data$pi.h,paste(output_folder,"/",samplename,"_2D_iters",noiters,"_concParam",conc_param,"_clusterWidth",1/cluster_conc,"_discreteMutationCopyNumbers.csv",sep=""))
write.csv(GS.data$alpha,paste(output_folder,"/",samplename,"_2D_iters",noiters,"_concParam",conc_param,"_clusterWidth",1/cluster_conc,"_alpha.csv",sep=""))

for(i in 1:(length(subsamples[[run]])-1)){
	for(j in (i+1):length(subsamples[[run]])){
		imageFile = paste(output_folder,"/",samplenames[run],subsamples[[run]][i],subsamples[[run]][j],"_iters",noiters,"_concParam",conc_param,"_clusterWidth",1/cluster_conc,"_2D_binomial.png",sep="")
		Gibbs.subclone.density.est(mutation.copy.number[,c(i,j)]/copyNumberAdjustment[,c(i,j)],GS.data,imageFile, post.burn.in.start = noiters*0.2, post.burn.in.stop = noiters, samplenames = paste(samplename,subsamples[[run]][c(i,j)],sep=""),indices=c(i,j))		
	}
}
