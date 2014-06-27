source("Tree_based_DP_Gibbs_sampler.R")
source("GetConsensusTrees.R")
source("PlotTreeWithIgraph.R")
source("interconvertMutationBurdens.R")
source("AnnotateTree.R")
source("InformationCriterions.R")

library(compiler)
mutationCopyNumberToMutationBurden = cmpfun(mutationCopyNumberToMutationBurden)

library(foreach)
library(doParallel)
library(doRNG)
# library(snowfall)
library(snow)

RunTreeBasedDP<-function(mutCount, WTCount, cellularity = rep(1,ncol(mutCount)), kappa = array(0.5,dim(mutCount)), samplename = "sample", subsamplenames = 1:ncol(mutCount), annotation = vector(mode="character",length=nrow(mutCount)), no.iters = 1250, no.iters.burn.in = 250, bin.size = NA, resort.mutations = T, outdir = paste(samplename,"_treeBasedDirichletProcessOutputs",sep=""), init.alpha = 0.01, shrinkage.threshold = 0.1, remove.node.frequency = NA, remove.branch.frequency = NA, parallel=FALSE, phase=NA){

  start = Sys.time()

  # Setup threads for parallel computations
  if(parallel) {
    clp = makeCluster(4)
    registerDoParallel(clp)
  } else {
    clp = NA
  }
  
  ###################################
  # Read in the data files
  ###################################
	no.subsamples = ncol(mutCount)

	#aggregate mutations that have similar allele burden and the same kappa
	if(!is.na(bin.size)){
		unique.kappa = unique(kappa)
		allele.fractions = mutCount / (mutCount + WTCount)
		rounded.allele.fractions = floor(allele.fractions / bin.size) * bin.size
    
		binned.mutCount = array(NA,c(0,no.subsamples))
		binned.WTCount = array(NA,c(0,no.subsamples))
		binned.kappa = array(NA,c(0,no.subsamples))
		bin.indices = list()
    
		for(u in 1:nrow(unique.kappa)){
      if(parallel) {
        inds = which(parSapply(cl=clp,X=1:nrow(kappa),FUN=function(k,i,u1){all(k[i,]==unique.kappa[u1,])},k=kappa, u1=u, simplify=TRUE, USE.NAMES=TRUE))
      } else {
        inds = which(sapply(1:nrow(kappa),function(k,i,u1){all(k[i,]==unique.kappa[u,])},k=kappa,u1=u))
      }
			unique.AF = unique(array(rounded.allele.fractions[inds,],c(length(inds),no.subsamples)))
			for(v in 1:nrow(unique.AF)){
        if(parallel && length(inds) > 1) {
          inds2 = inds[parSapply(cl=clp,X=inds,function(r,i,v1){all(r[i,]==unique.AF[v1,])},r=rounded.allele.fractions, v1=v, simplify=TRUE, USE.NAMES=TRUE)]
        } else {
          inds2 = inds[sapply(inds,function(r,i,v1){all(r[i,]==unique.AF[v1,])},r=rounded.allele.fractions, v1=v)]
        }
				bin.indices[[length(bin.indices)+1]] = inds2
				if(length(inds2)>1){
					binned.mutCount = rbind(binned.mutCount,colSums(mutCount[inds2,]))
					binned.WTCount = rbind(binned.WTCount,colSums(WTCount[inds2,]))
					binned.kappa = rbind(binned.kappa,unique.kappa[u,])
				}else{
					binned.mutCount = rbind(binned.mutCount,mutCount[inds2,])
					binned.WTCount = rbind(binned.WTCount,WTCount[inds2,])
					binned.kappa = rbind(binned.kappa,unique.kappa[u,])			
				}
			}
		}
    
		no.muts = nrow(binned.mutCount)
	}else{
		no.muts = nrow(mutCount)
	}

	if(!file.exists(outdir)){
		dir.create(outdir)
	}

  ###################################
  # Run the tree phase of the method
  ###################################
  if (is.na(phase) | phase == 'tree') {
  
  	if(is.na(bin.size)){
  		temp.list = tree.struct.dirichlet.gibbs(y=mutCount,n=WTCount+mutCount,kappa=kappa,iter=no.iters,shrinkage.threshold=shrinkage.threshold,init.alpha=init.alpha, remove.node.frequency = NA, remove.branch.frequency = NA, parallel=parallel, cluster=clp)
  	}else{
  		save(bin.indices,file=paste(outdir,"/",samplename,"_",no.iters,"iters_burnin",no.iters.burn.in,"_binnedIndices.RData",sep=""))
  		temp.list = tree.struct.dirichlet.gibbs(y=binned.mutCount,n=binned.WTCount+binned.mutCount,kappa=binned.kappa,iter=no.iters,shrinkage.threshold=shrinkage.threshold,init.alpha=init.alpha, remove.node.frequency = NA, remove.branch.frequency = NA, parallel=parallel, cluster=clp)	
  	}
  
    print(paste("Finished tree.struct.dirichlet in", as.numeric(Sys.time()-start,units="secs"), "seconds"))

    setwd(outdir)

  	trees = temp.list[[1]]
  	if(is.na(bin.size)){
  		node.assignments = temp.list[[2]]
  	}else{
  		binned.node.assignments = temp.list[[2]]
  		node.assignments = array(NA,c(nrow(mutCount),ncol(binned.node.assignments)))
  		for(i in 1:length(bin.indices)){
  			node.assignments[bin.indices[[i]],] = binned.node.assignments[i,]
  		}
  		write.table(binned.node.assignments,paste("aggregated_node_assignments_",samplename,"_",no.iters,"iters.txt",sep=""),sep="\t",row.names=F,quote=F)
  	}
  	alphas = temp.list[[3]]
  	lambdas = temp.list[[4]]
  	gammas = temp.list[[5]]
  	likelihoods = temp.list[[6]]
  	BIC = temp.list[[7]]
    AIC = temp.list[[8]]
    DIC = temp.list[[9]]
    if(parallel) {
      tree.sizes = parSapply(cl=clp,X=1:length(trees),function(t,i){nrow(t[[i]])},t=trees, simplify=TRUE, USE.NAMES=TRUE)
    } else {
      tree.sizes = sapply(1:length(trees),function(t,i){nrow(t[[i]])},t=trees)
    }
  	best.index = which.max(likelihoods)
  	best.BIC.index = which.min(BIC)
    best.AIC.index = which.min(AIC)
    best.DIC.index = which.min(DIC)
  
  	pdf(paste("all_trees_",samplename,"_",no.iters,"iters.pdf",sep=""),height=4,width=3*no.subsamples)
  	for(iter in 1:no.iters){
  		tree = trees[[iter]]
  		#tree$annotation = NA
  		tree = annotateTree(tree,node.assignments[,iter],annotation)
  		title = paste(samplename,"iter",iter)
  		if(iter==best.index){
  			title = paste(title,"(most likely tree)")
  		}
  		if(iter==best.BIC.index){
  			title = paste(title,"(best BIC tree)")
  		}
      if(iter==best.AIC.index) {
        title = paste(title,"(best AIC tree)")
      }
  		if(iter==best.DIC.index) {
  		  title = paste(title,"(best DIC tree)")
  		}
  		plotTree(tree,title)
  	}
  	dev.off()
  	
  	png(paste("log_likelihood_plot_",samplename,"_",no.iters,"iters.png",sep=""),width=1000)
  	par(cex.lab=3,cex.axis=3,lwd=3,mar=c(7,7,5,2))
  	plot(1:no.iters,likelihoods,type="l",col="red",xlab="MCMC iteration", ylab="log-likelihood",main=samplename)
  	dev.off()
  	png(paste("BIC_plot_",samplename,"_",no.iters,"iters.png",sep=""),width=1000)
  	par(cex.lab=3,cex.axis=3,lwd=3,mar=c(7,7,5,2))
  	plot(1:no.iters,BIC,type="l",col="red",xlab="MCMC iteration", ylab="BIC",main=samplename)
  	dev.off()
  	png(paste("AIC_plot_",samplename,"_",no.iters,"iters.png",sep=""),width=1000)
  	par(cex.lab=3,cex.axis=3,lwd=3,mar=c(7,7,5,2))
  	plot(1:no.iters,AIC,type="l",col="red",xlab="MCMC iteration", ylab="AIC",main=samplename)
  	dev.off()
  	png(paste("DIC_plot_",samplename,"_",no.iters,"iters.png",sep=""),width=1000)
  	par(cex.lab=3,cex.axis=3,lwd=3,mar=c(7,7,5,2))
  	plot(1:no.iters,DIC,type="l",col="red",xlab="MCMC iteration", ylab="DIC",main=samplename)
  	dev.off()
  	
  	write.table(node.assignments,paste("node_assignments_",samplename,"_",no.iters,"iters.txt",sep=""),sep="\t",row.names=F,quote=F)
  	write.table(alphas,paste("alphas_",samplename,"_",no.iters,"iters.txt",sep=""),sep="\t",row.names=F,quote=F)
  	write.table(lambdas,paste("lambdas_",samplename,"_",no.iters,"iters.txt",sep=""),sep="\t",row.names=F,quote=F)
  	write.table(gammas,paste("gammas_",samplename,"_",no.iters,"iters.txt",sep=""),sep="\t",row.names=F,quote=F)
  	write.table(likelihoods,paste("likelihoods_",samplename,"_",no.iters,"iters.txt",sep=""),sep="\t",row.names=F,quote=F)
  	#write.table(best.tree,paste("best_tree_",samplename,"_",no.iters,"iters.txt",sep=""),sep="\t",row.names=F,quote=F)
  	save(trees,file=paste("",samplename,"_trees_iters",no.iters,".Rdata",sep=""))
  
  } else {
    ## Load the data
    load(file=paste("",samplename,"_trees_iters",no.iters,".Rdata",sep=""))
    node.assignments = read.table(paste("node_assignments_",samplename,"_",no.iters,"iters.txt",sep=""),sep="\t",row.names=F,quote=F)
    alphas = read.table(paste("alphas_",samplename,"_",no.iters,"iters.txt",sep=""),sep="\t",row.names=F,quote=F)
    lambdas = read.table(lambdas,paste("lambdas_",samplename,"_",no.iters,"iters.txt",sep=""),sep="\t",row.names=F,quote=F)
    gammas = read.table(gammas,paste("gammas_",samplename,"_",no.iters,"iters.txt",sep=""),sep="\t",row.names=F,quote=F)
    likelihoods = read.table(paste("likelihoods_",samplename,"_",no.iters,"iters.txt",sep=""),sep="\t",row.names=F,quote=F)
  } 
	
  if (is.na(phase) | phase == 'cons') {
  	if(is.na(bin.size)){
  		temp.list = GetConsensusTrees(trees, node.assignments, mutCount, WTCount, kappa = kappa, samplename = samplename, subsamplenames = subsamplenames, no.iters = no.iters, no.iters.burn.in = no.iters.burn.in, resort.mutations = resort.mutations)
  	}else{
  		temp.list = GetConsensusTrees(trees, binned.node.assignments, binned.mutCount, binned.WTCount, kappa = binned.kappa, samplename = samplename, subsamplenames = subsamplenames, no.iters = no.iters, no.iters.burn.in = no.iters.burn.in, resort.mutations = resort.mutations,bin.indices = bin.indices)
  	}
  	print("Finished GetConsensusTrees")
  	print(Sys.time()-start)
  	print(names(temp.list))
  	
  	all.consensus.trees = temp.list$all.consensus.trees
  	save(all.consensus.trees,file=paste(samplename,"_",no.iters,"iters_",no.iters.burn.in,"_allConsensusTrees.RData",sep=""))
  	if(!is.na(bin.size)){
  		all.consensus.assignments = temp.list$all.consensus.assignments
  		all.disaggregated.consensus.assignments = temp.list$all.disaggregated.consensus.assignments
  		save(all.consensus.assignments,file=paste(samplename,"_",no.iters,"iters_",no.iters.burn.in,"_allBinnedConsensusAssignments.RData",sep=""))
  		save(all.disaggregated.consensus.assignments,file=paste(samplename,"_",no.iters,"iters_",no.iters.burn.in,"_allConsensusAssignments.RData",sep=""))		
  	}else{
  		all.consensus.assignments = temp.list$all.consensus.assignments
  		save(all.consensus.assignments,file=paste(samplename,"_",no.iters,"iters_",no.iters.burn.in,"_allConsensusAssignments.RData",sep=""))
  	}
  	likelihoods = temp.list$likelihoods
  	BIC = temp.list$BIC
    AIC = temp.list$AIC
    DIC = temp.list$DIC
  	plotScores("log_likelihood", samplename, no,iters, likelihoods, "log-likelihood")
  	plotScores("BIC", samplename, no,iters, BIC, "BIC")
  	plotScores("AIC", samplename, no,iters, AIC, "AIC")
  	plotScores("DIC", samplename, no,iters, BIC, "DIC")
#   	png(paste("log_likelihood_plot_",samplename,"_consensus_",no.iters,"iters.png",sep=""),width=1000)
#   	par(cex.lab=3,cex.axis=3,lwd=3,mar=c(7,7,5,2))
#   	plot(1:length(likelihoods),likelihoods,type="l",col="red",xlab="posterior tree index", ylab="log-likelihood",main=samplename)
#   	dev.off()
#   	png(paste("BIC_plot_",samplename,"_consensus_",no.iters,"iters.png",sep=""),width=1000)
#   	par(cex.lab=3,cex.axis=3,lwd=3,mar=c(7,7,5,2))
#   	plot(1:length(BIC),BIC,type="l",col="red",xlab="posterior tree index", ylab="BIC",main=samplename)
#   	dev.off()
    writeScoresTable(likelihoods, "likelihood", samplename, no,iters)
    writeScoresTable(BIC, "BIC", samplename, no,iters)
    writeScoresTable(AIC, "AIC", samplename, no,iters)
    writeScoresTable(DIC, "DIC", samplename, no,iters)
#   	write.table(data.frame(tree.index = 1:length(likelihoods),likelihood=likelihoods),paste("consensus_likelihoods_",samplename,"_",no.iters,"iters.txt",sep=""),sep="\t",row.names=F,quote=F)
#   	write.table(data.frame(tree.index = 1:length(BIC),BIC=BIC),paste("consensus_BICs_",samplename,"_",no.iters,"iters.txt",sep=""),sep="\t",row.names=F,quote=F)
  	
  	best.index = which.min(BIC)
  	best.tree = all.consensus.trees[[best.index]]
  	if(!is.na(bin.size)){
  		best.node.assignments = all.disaggregated.consensus.assignments[[best.index]]
  	}else{
  		best.node.assignments = all.consensus.assignments[[best.index]]
  	}
  	
  	#png(paste("bestConsensusTree_",samplename,"_",no.iters,"iters.png",sep=""))
  	png(paste("bestConsensusTree_",samplename,"_",no.iters,"iters.png",sep=""),width=no.subsamples*1000,height=1000)
  	if(!is.na(bin.size)){
  		plotTree(annotateTree(best.tree,all.disaggregated.consensus.assignments,annotation),font.size=2,main=paste(samplename,subsamplenames,sep=""))
  		
  	}else{
  		plotTree(annotateTree(best.tree,all.consensus.assignments,annotation),font.size=2,main=paste(samplename,subsamplenames,sep=""))
  	}
  	dev.off()
  	write.table(best.node.assignments,paste(samplename,"_",no.iters,"iters_",no.iters.burn.in,"_bestConsensusAssignments.txt",sep=""),sep="\t",quote=F,row.names=F,col.names=F)
  	
  	if(no.subsamples>1){
  	  plotBestScatter(mutCount, WTCount, kappa, best.tree, samplename, no.iters)
#   		subclonal.fraction = array(NA,dim(mutCount))
#   		for(i in 1:no.subsamples){
#   			subclonal.fraction[,i] = mutCount[,i] / ((mutCount[,i]+WTCount[,i])*kappa[,i])
#   			subclonal.fraction[kappa[,i]==0,i] = NA
#   		}	
#   		
#   		#its hard to distinguish more than 8 different colours
#   		max.cols = 8
#   		no.nodes = nrow(best.tree)
#   		unique.nodes = rownames(best.tree)
#   		cols = rainbow(min(max.cols,no.nodes))
#   		png(paste("bestScatterPlots_",samplename,"_",no.iters,"iters.png",sep=""),width=no.subsamples*(no.subsamples-1)*500,height=1000)
#   		par(mfrow = c(1,no.subsamples*(no.subsamples-1)/2),cex=2)
#   		plot.data = subclonal.fraction
#   		plot.data[is.na(plot.data)]=0
#   		for(i in 1:(no.subsamples-1)){
#   			for(j in (i+1):no.subsamples){
#   				plot(plot.data[,i],plot.data[,j],type = "n",xlab = paste(samplename,subsamplenames[i]," subclonal fraction",sep=""), ylab = paste(samplename,subsamplenames[j]," subclonal fraction",sep=""),xlim = c(0,max(plot.data[,i])*1.3))
#   				for(n in 1:no.nodes){
#   					points(plot.data[,i][best.node.assignments==unique.nodes[n]],plot.data[,j][best.node.assignments==unique.nodes[n]],col=cols[(n-1) %% max.cols + 1],pch=20 + floor((n-1)/max.cols))
#   				}
#   				#legend(max(plot.data[,i]),max(plot.data[,j]),legend = unique.nodes[n],col=cols[(n-1) %% max.cols + 1],pch=20 + floor((n-1)/max.cols))
#   				legend(max(plot.data[,i])*1.1,max(plot.data[,j]),legend = unique.nodes,col=cols[(0:(no.nodes-1)) %% max.cols + 1],pch=20 + floor((0:(no.nodes-1))/max.cols))
#   			}
#   		}
#   		dev.off()
  	}	
  }

  if(parallel) {
    stopCluster(clp)
  }
}

plotScores <- function(file_prefix, samplename, no,iters, scores, ylab) {
#   png(paste("log_likelihood_plot_",samplename,"_consensus_",no.iters,"iters.png",sep=""),width=1000)
  png(paste(file_prefix,"_plot_",samplename,"_consensus_",no.iters,"iters.png",sep=""),width=1000)
  par(cex.lab=3,cex.axis=3,lwd=3,mar=c(7,7,5,2))
#   plot(1:length(likelihoods),likelihoods,type="l",col="red",xlab="posterior tree index", ylab="log-likelihood",main=samplename)
  plot(1:length(scores),scores,type="l",col="red",xlab="posterior tree index", ylab=ylab,main=samplename)
  dev.off()
}

writeScoresTable <- function(scores, score_name, samplename, no,iters) {
#   write.table(data.frame(tree.index = 1:length(likelihoods),likelihood=likelihoods),paste("consensus_likelihoods_",samplename,"_",no.iters,"iters.txt",sep=""),sep="\t",row.names=F,quote=F)
  data = data.frame(tree.index = 1:length(scores),likelihood=scores)
  colnames(data) = c("tree.index", score_name)
  write.table(data,paste("consensus_",score_name,"_",samplename,"_",no.iters,"iters.txt",sep=""),sep="\t",row.names=F,quote=F)
}

plotBestScatter <- function(mutCount, WTCount, kappa, best.tree, samplename, no.iters) {
  subclonal.fraction = array(NA,dim(mutCount))
  for(i in 1:no.subsamples){
    subclonal.fraction[,i] = mutCount[,i] / ((mutCount[,i]+WTCount[,i])*kappa[,i])
    subclonal.fraction[kappa[,i]==0,i] = NA
  }	
  
  #its hard to distinguish more than 8 different colours
  max.cols = 8
  no.nodes = nrow(best.tree)
  unique.nodes = rownames(best.tree)
  cols = rainbow(min(max.cols,no.nodes))
  png(paste("bestScatterPlots_",samplename,"_",no.iters,"iters.png",sep=""),width=no.subsamples*(no.subsamples-1)*500,height=1000)
  par(mfrow = c(1,no.subsamples*(no.subsamples-1)/2),cex=2)
  plot.data = subclonal.fraction
  plot.data[is.na(plot.data)]=0
  for(i in 1:(no.subsamples-1)){
    for(j in (i+1):no.subsamples){
      plot(plot.data[,i],plot.data[,j],type = "n",xlab = paste(samplename,subsamplenames[i]," subclonal fraction",sep=""), ylab = paste(samplename,subsamplenames[j]," subclonal fraction",sep=""),xlim = c(0,max(plot.data[,i])*1.3))
      for(n in 1:no.nodes){
        points(plot.data[,i][best.node.assignments==unique.nodes[n]],plot.data[,j][best.node.assignments==unique.nodes[n]],col=cols[(n-1) %% max.cols + 1],pch=20 + floor((n-1)/max.cols))
      }
      #legend(max(plot.data[,i]),max(plot.data[,j]),legend = unique.nodes[n],col=cols[(n-1) %% max.cols + 1],pch=20 + floor((n-1)/max.cols))
      legend(max(plot.data[,i])*1.1,max(plot.data[,j]),legend = unique.nodes,col=cols[(0:(no.nodes-1)) %% max.cols + 1],pch=20 + floor((0:(no.nodes-1))/max.cols))
    }
  }
  dev.off()
}

