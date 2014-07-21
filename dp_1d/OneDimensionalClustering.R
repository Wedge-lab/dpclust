args=(commandArgs(TRUE))
run = 1
if(length(args)>0){
	run = as.integer(args[1])
}

source("/lustre/scratch109/sanger/dw9/2D_Dirichlet_Process/interconvertMutationBurdens.R")

setwd("/lustre/scratch109/sanger/dw9/2D_Dirichlet_Process/melanoma_cell_lines")
samplenames=c("CP50-MEL-B","CP66-MEL","LB2518-MEL","LB373-MEL-D","MZ7-mel")
samplename = samplenames[run]

output_folder="melanoma_cell_lines_DP_output"

no.iters = 1000
burn.in = 200

log.f.of.y <- function(y1, n1, kappa1, x) {
	#x=1 and kappa=1 causes problems
	x[x>0.999 & kappa1==1] = 0.999
	#allow kappa = 0, for mutations on deleted chromosomes
	non.zero.inds = which(kappa1!=0)
	if(length(non.zero.inds)>0){
		return(sum(lchoose(n1[non.zero.inds], y1[non.zero.inds]) + y1[non.zero.inds] * log(kappa1[non.zero.inds]*x[non.zero.inds]) + (n1[non.zero.inds]-y1[non.zero.inds]) * log(1-kappa1[non.zero.inds]*x[non.zero.inds])) * length(y1)/length(non.zero.inds))
	}else{
		print("WARNING. ALL KAPPAS ARE ZERO. MUTATION ABSENT FROM ALL SAMPLES")
		print(y1)
		print(n1)
		print(kappa1)
		print(x)
		#return(NA)
		return(NaN)
	}
}

mutdata=read.table(paste(samplename,"_allDirichletProcessInfo.txt",sep=""),sep="\t",header=T,stringsAsFactors=F)
mutdata = mutdata[!is.na(mutdata$subclonal.fraction),]

mutdata = mutdata[mutdata$CHR %in% 1:22,]
mutCount = mutdata$mut.count
WTCount = mutdata$WT.count
totalReads = mutCount + WTCount

no.muts = length(mutCount)
total.copy.number = mutdata$subclonal.CN
normal.copy.number = rep(2,no.muts)
mutation.copy.number = mutdata$mutation.copy.number
no.chrs.bearing.mut = mutdata$no.chrs.bearing.mut
subclonal.fraction = mutation.copy.number / no.chrs.bearing.mut

node.assignments = read.table(paste("melanoma_cell_lines_DP_output/",samplename,"_iters",no.iters,"_states.csv",sep=""),sep=",",row.names=1, header=T)
no.muts = ncol(node.assignments)

library(MASS)
library(MCMCpack)
library(mvtnorm)
source("/lustre/scratch109/sanger/dw9/2D_Dirichlet_Process/multidimensionalDensityEstimator.R")
source("/lustre/scratch109/sanger/dw9/2D_Dirichlet_Process/interconvertMutationBurdens.R")


S.i = read.csv(paste(output_folder,"/",samplename,"_iters",no.iters,"_states.csv",sep=""),row.names=1)
V.h = read.csv(paste(output_folder,"/",samplename,"_iters",no.iters,"_stickbreaking_weights.csv",sep=""),row.names=1)
pi.h = read.csv(paste(output_folder,"/",samplename,"_iters",no.iters,"_discreteMutationCopyNumbers.csv",sep=""),row.names=1)

density = read.table(paste(output_folder,"/",samplename,"_iters",no.iters,"density.txt",sep=""),row.names=NULL,header=T,sep="\t")

hypercube.size = 20
localOptima = NULL
peak.indices = NULL
for(i in (1+hypercube.size):(nrow(density)-hypercube.size)){
	if(density$median.density[i] == max(density$median.density[(i-hypercube.size):(i+hypercube.size)])){
		localOptima = c(localOptima,density$MCN[i])
		peak.indices = c(peak.indices,i)
	}
}

print("localOptima")
print(localOptima)

write.table(localOptima,paste(samplename,"_localOptima.txt",sep=""),quote=F,sep="\t")


#ASSIGN mutations to clusters
no.optima = length(localOptima)
if(no.optima>1){	
	boundary = array(NA,no.optima-1)
	mutation.preferences = array(0,c(no.muts,no.optima))
	for(i in 1:(no.optima-1)){
		min.density = min(density$median.density[(peak.indices[i]+1):(peak.indices[i+1]-1)])
		min.indices = intersect(which(density$median.density == min.density),(peak.indices[i]+1):(peak.indices[i+1]-1))
		
		#what distance along the line between a pair of optima do we have to go to reach the minimum density
		boundary[i] = (density$MCN[max(min.indices)] + density$MCN[min(min.indices)])/2
	}

	sampledIters = (burn.in + 1) : no.iters
	#don't use the intitial state
	sampledIters = sampledIters[sampledIters!=1]		
	if(length(sampledIters) > 1000){
		sampledIters=floor(post.burn.in.start + (1:1000) * (no.iters - burn.in)/1000)			
	}

	S.i = data.matrix(S.i)
	for(s in sampledIters){
		temp.preferences = array(0,c(no.muts,no.optima))
		for(c in unique(S.i[s,])){
			bestOptimum = sum(pi.h[s,c]>boundary)+1
			temp.preferences[S.i[s,]==c,bestOptimum] = temp.preferences[S.i[s,]==c,bestOptimum] + 1		
		}
		iter.preferences = t(sapply(1:no.muts,function(p,i){as.integer(p[i,]==max(p[i,])) / sum(p[i,]==max(p[i,]))},p=temp.preferences))
		mutation.preferences = mutation.preferences + iter.preferences
	}
	mutation.preferences = mutation.preferences / length(sampledIters)
	most.likely.cluster = sapply(1:no.muts,function(m,i){which.max(m[i,])},m=mutation.preferences)
	
	out = cbind(mutdata,mutation.preferences,most.likely.cluster)
	names(out)[(ncol(out)-no.optima):ncol(out)] = c(paste("prob.cluster",1:no.optima,sep=""),"most.likely.cluster")
	write.table(cbind(1:no.optima,localOptima,colSums(mutation.preferences)),paste(samplename,"_optimaInfo.txt",sep=""),col.names = c("cluster.no","MCN","no.of.mutations"),row.names=F,sep="\t",quote=F)		
	write.table(out,paste(samplename,"_DP_and cluster_info.txt",sep=""),sep="\t",row.names=F,quote=F)

}else{
	most.likely.cluster = rep(1,no.muts)
}

q(save="no")
