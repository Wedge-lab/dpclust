source("Tree_based_DP_Gibbs_sampler.R")
source("GetConsensusTrees.R")
source("PlotTreeWithIgraph.R")
source("interconvertMutationBurdens.R")
source("AnnotateTree.R")
source("InformationCriterions.R")
source("RunTreeBasedDPConsensus.R")
source("RunTreeBasedDPMCMC.R")

library(compiler)
mutationCopyNumberToMutationBurden = cmpfun(mutationCopyNumberToMutationBurden)

library(foreach)
library(doParallel)
library(doRNG)
# library(snowfall)
library(snow)

RunTreeBasedDP<-function(mutCount, WTCount, cellularity = rep(1,ncol(mutCount)), kappa = array(0.5,dim(mutCount)), samplename = "sample", subsamplenames = 1:ncol(mutCount), annotation = vector(mode="character",length=nrow(mutCount)), no.iters = 1250, no.iters.burn.in = 250, bin.size = NA, resort.mutations = T, outdir = paste(samplename,"_treeBasedDirichletProcessOutputs",sep=""), init.alpha = 0.01, shrinkage.threshold = 0.1, remove.node.frequency = NA, remove.branch.frequency = NA, parallel=FALSE, phase=NA, blockid=1, no.of.blocks=NULL){


  seeds = c(123,321,213,231,456,654,465,645,789,987,978,798)
  set.seed(seeds[blockid])


  start_time = Sys.time()
  
  if(!file.exists(outdir)){
    dir.create(outdir)
  }

  ###################################
  # Perform binning
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
      inds = which(apply(kappa, 1, function(k,uk){all(k==uk)},uk=unique.kappa[u,]))
			unique.AF = unique(array(rounded.allele.fractions[inds,],c(length(inds),no.subsamples)))
			for(v in 1:nrow(unique.AF)){
        inds2 = inds[apply(matrix(rounded.allele.fractions[inds,], nrow=length(inds)), 1, function(r,ua){all(r==ua)},ua=unique.AF[v,])]
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
		save(bin.indices,file=paste(outdir,"/",samplename,"_",no.iters,"iters_burnin",no.iters.burn.in,"_binnedIndices.RData",sep=""))
    
		no.muts = nrow(binned.mutCount)
	}else{
		no.muts = nrow(mutCount)
	}

  ###################################
  # Run the tree phase of the method
  ###################################
  if (is.na(phase) | phase == 'tree') {
    # Setup threads for parallel computations
    if(parallel) {
      clp = makeCluster(4)
      registerDoParallel(clp)
    } else {
      clp = NA
    }    
    
    trees_start_time = Sys.time()
    if(!is.na(bin.size)) {
      RunTreeBasedDPMCMC(binned.mutCount, binned.WTCount, binned.kappa, nrow(mutCount), annotation, no.iters, shrinkage.threshold, init.alpha, outdir, parallel, clp, bin.indices=bin.indices, blockid=blockid)
    } else {
      RunTreeBasedDPMCMC(mutCount, WTCount, kappa, nrow(mutCount), annotation, no.iters, shrinkage.threshold, init.alpha, outdir, parallel, clp, bin.indices=NULL, blockid=blockid)
    }
    trees_end_time = Sys.time()
    print(paste("Finished tree.struct.dirichlet in", as.numeric(trees_end_time-trees_start_time,units="secs"), "seconds"))
    write.table(data.frame(diff=c(difftime(trees_end_time, trees_start_time, units='sec')), unit=c('seconds')), file=paste(outdir,'runtime_trees_',blockid,'.txt', sep=''), quote=F, row.names=F)
  
    if(parallel) { stopCluster(clp) }
    
  } else {
    # Change working dir for the cons phase
    setwd(outdir)
  }
  
  ###################################
  # Run the cons phase of the method
  ###################################
  if (is.na(phase) | phase == 'cons') {
    ###################################
    # Load data
    ###################################
    trees_all <- vector("list", length=0)
    node.assignments <- matrix(NA, nrow=nrow(mutCount), ncol=0)
    likelihoods = vector(mode='numeric',length=0)
    alphas = vector(mode='numeric',length=0)
    lambdas = vector(mode='numeric',length=0)
    gammas = vector(mode='numeric',length=0)
    if(!is.na(bin.size)) {
      binned.node.assignments = matrix(NA, nrow=nrow(binned.mutCount), ncol=0)
    }
    
    ## TODO: Work out how to merge block-parallel data
    for(i in 1:no.of.blocks) {
    
      ## list, dataframe per iter
      load(file=paste("",samplename,"_trees_iters",no.iters,"_block",blockid,".Rdata",sep=""))
      trees_all = append(trees_all, trees)
      print("assignments")
      ## Col per iter
      node.assignments_block = read.table(paste("node_assignments_",samplename,"_",no.iters,"iters_block",blockid,".txt",sep=""),sep="\t",header=T)
      node.assignments = cbind(node.assignments,node.assignments_block)

      ## Four below: List with item per iter
      alphas_block = read.table(paste("alphas_",samplename,"_",no.iters,"iters_block",blockid,".txt",sep=""),sep="\t")
      alphas = c(alphas, alphas_block)

      lambdas_block = read.table(paste("lambdas_",samplename,"_",no.iters,"iters_block",blockid,".txt",sep=""),sep="\t")
      lambdas = c(lambdas, lambdas_block)

      gammas_block = read.table(paste("gammas_",samplename,"_",no.iters,"iters_block",blockid,".txt",sep=""),sep="\t")
      gammas = c(gammas, gammas_block)

      likelihoods_block = read.table(paste("likelihoods_",samplename,"_",no.iters,"iters_block",blockid,".txt",sep=""),sep="\t")
      likelihoods = c(likelihoods, likelihoods_block)

      if(!is.na(bin.size)) {
        ## Col per iter
        binned.node.assignments_block = read.table(paste("aggregated_node_assignments_",samplename,"_",no.iters,"iters_block",blockid,".txt",sep=""),sep="\t",header=T)
        binned.node.assignments = cbind(binned.node.assignments, binned.node.assignments_block)
      }
    }
	
    ###################################
    # Start Cons
    ###################################
    cons_start_time = Sys.time()
    if(!is.na(bin.size)) { 
      RunTreeBasedDPConsensus(trees_all, binned.node.assignments, binned.mutCount, binned.WTCount, binned.kappa, samplename, subsamplenames, annotation, no.iters, no.iters.burn.in, resort.mutations, bin.indices=bin.indices)
    } else {
      RunTreeBasedDPConsensus(trees_all, node.assignments, mutCount, WTCount, kappa, samplename, subsamplenames, annotation, no.iters, no.iters.burn.in, resort.mutations, bin.indices=NULL)
    }
    cons_end_time = Sys.time()
    print(paste("Finished RunTreeBasedDPConsensus in", as.numeric(cons_end_time-cons_start_time,units="secs"), "seconds"))
    write.table(data.frame(diff=c(difftime(cons_end_time, cons_start_time, units='sec')), unit=c('seconds')), file='runtime_cons.txt', quote=F, row.names=F)
  }

  if (is.na(phase)) {
    end_time = Sys.time()
    # This file should always be written to the outdir
    write.table(data.frame(diff=c(difftime(end_time, start_time, units='sec')), unit=c('seconds')), file='runtime.txt', quote=F, row.names=F)
    print(paste("Finished in", as.numeric(end_time-start_time,units="secs"), "seconds"))
  }
}
