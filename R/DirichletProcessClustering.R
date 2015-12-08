# nD method includes
# source('subclone_Dirichlet_Gibbs_sampler_nD_binomial.R')
# source("OneDimensionalClustering.R")
# 
# # Tree based method includes
# source("Tree_based_DP_Gibbs_sampler.R")
# source("GetConsensusTrees.R")
# source("PlotTreeWithIgraph.R")
# source("AnnotateTree.R")
# source("InformationCriterions.R")
# source("RunTreeBasedDPConsensus.R")
# source("RunTreeBasedDPMCMC.R")
# 
# # Shared includes
# source("interconvertMutationBurdens.R")
# 
# library(doParallel)
# library(doRNG)

TreeBasedDP<-function(mutCount, WTCount, removed_indices=c(), cellularity=rep(1,ncol(mutCount)), kappa=array(0.5,dim(mutCount)), samplename="sample", subsamplenames=1:ncol(mutCount), annotation=vector(mode="character",length=nrow(mutCount)), no.iters=1250, no.iters.burn.in=250, bin.size=NA, resort.mutations=T, outdir=paste(samplename,"_treeBasedDirichletProcessOutputs",sep=""), init.alpha=0.01, shrinkage.threshold=0.1, remove.node.frequency=NA, remove.branch.frequency=NA, parallel=FALSE, phase=NA, blockid=1, no.of.blocks=NULL, conflict.array=NA){
  #
  # Tree based method that will yield a tree that contains the estimated clone/subclone structure for 
  # the samples given as input.
  #
  # It consists of two phases:
  # * tree : A random walk that creates a series of random trees
  # * cons : The random trees are summarised in a new series of consensus trees. The best consensus tree is selected as the resulting tree. Cons can handle input from multiple parallel trees calls.
  #
  
  # Set the seed for each block
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
  if (is.na(phase) | phase == 'tree_dp' | phase == 'tree') {
    # Setup threads for parallel computations
    if(parallel) {
      clp = makeCluster(4)
      registerDoParallel(clp)
    } else {
      clp = NA
    }    
    
    trees_start_time = Sys.time()
    if(!is.na(bin.size)) {
      mcmcResults = RunTreeBasedDPMCMC(mutCount=binned.mutCount, 
                         WTCount=binned.WTCount, 
                         kappa=binned.kappa, 
                         no.muts.input=nrow(mutCount), 
                         annotation=annotation,
                         samplename=samplename,
                         no.iters=no.iters, 
                         no.iters.burn.in=no.iters.burn.in,
                         shrinkage.threshold=shrinkage.threshold, 
                         init.alpha=init.alpha, 
                         outdir=outdir, 
                         parallel=parallel, 
                         clp=clp, 
                         bin.indices=bin.indices, 
                         blockid=blockid,
                         remove.node.frequency=remove.node.frequency,
                         remove.branch.frequency=remove.branch.frequency,
                         conflict.array=conflict.array)
    } else {
      mcmcResults = RunTreeBasedDPMCMC(mutCount=mutCount, 
                         WTCount=WTCount, 
                         kappa=kappa, 
                         no.muts.input=nrow(mutCount), 
                         annotation=annotation, 
                         samplename=samplename,
                         no.iters=no.iters, 
                         no.iters.burn.in=no.iters.burn.in,
                         shrinkage.threshold=shrinkage.threshold, 
                         init.alpha=init.alpha, 
                         outdir=outdir, 
                         parallel=parallel, 
                         clp=clp, 
                         bin.indices=NULL, 
                         blockid=blockid,
                         remove.node.frequency=remove.node.frequency,
                         remove.branch.frequency=remove.branch.frequency,
                         conflict.array=conflict.array)
    }
    
    # Set these variables in case both trees and cons need to be run consecutively in one R call
    trees_interleaved = mcmcResults$trees
    node.assignments_interleaved = mcmcResults$node.assignments
    binned.node.assignments_interleaved = mcmcResults$binned.node.assignments
    strengths = mcmcResults$strengths
    
    trees_end_time = Sys.time()
    print(paste("Finished tree.struct.dirichlet in", as.numeric(trees_end_time-trees_start_time,units="secs"), "seconds"))
    write.table(data.frame(diff=c(difftime(trees_end_time, trees_start_time, units='sec')), unit=c('seconds')), file=paste(outdir,'/runtime_trees_',blockid,'.txt', sep=''), quote=F, row.names=F)
    
    if(parallel) { stopCluster(clp) }
    
  } else {
    ###################################
    # Load data
    ###################################
    # Change working dir for the cons phase
    setwd(outdir)
    
    trees_all = vector("list", length=no.of.blocks)
    node.assignments_all = vector("list", length=no.of.blocks)
    ancestor.strengths_all = vector("list", length=no.of.blocks)
    sibling.strengths_all = vector("list", length=no.of.blocks)
    identity.strengths_all = vector("list", length=no.of.blocks)
    parent.child.strengths_all = vector("list", length=no.of.blocks)
    child.parent.strengths_all = vector("list", length=no.of.blocks)
    if(!is.na(bin.size)) {
      binned.node.assignments_all = vector("list", length=no.of.blocks)
    }
    
    for(i in 1:no.of.blocks) {
      load(file=paste("",samplename,"_trees_iters",no.iters,"_block",i,".Rdata",sep=""))
      ## list, dataframe per iter
      trees_all[[i]] = trees
      ## Col per iter
      node.assignments_all[[i]] = read.table(paste("node_assignments_",samplename,"_",no.iters,"iters_block",i,".txt",sep=""),sep="\t",header=T, stringsAsFactors=F)
      ancestor.strengths_all[[i]] = read.table(paste("ancestor.strengths_",samplename,"_",no.iters,"iters_block",i,".txt",sep=""),sep="\t",header=F)
      sibling.strengths_all[[i]] = read.table(paste("sibling.strengths_",samplename,"_",no.iters,"iters_block",i,".txt",sep=""),sep="\t",header=F)
      identity.strengths_all[[i]] = read.table(paste("identity.strengths_",samplename,"_",no.iters,"iters_block",i,".txt",sep=""),sep="\t",header=F)
      parent.child.strengths_all[[i]] = read.table(paste("parent.child.strengths_",samplename,"_",no.iters,"iters_block",i,".txt",sep=""),sep="\t",header=F)
      child.parent.strengths_all[[i]] = read.table(paste("child.parent.strengths_",samplename,"_",no.iters,"iters_block",i,".txt",sep=""),sep="\t",header=F)
      if(!is.na(bin.size)) {
        ## Col per iter
        binned.node.assignments_all[[i]] = read.table(paste("aggregated_node_assignments_",samplename,"_",no.iters,"iters_block",i,".txt",sep=""),sep="\t",header=T)
      }
    }
    
    # Storage for interleaving the data
    trees_appended = vector("list", length=0)
    node.assignments_appended = vector("list", length=0)
    if(!is.na(bin.size)) {
      binned.node.assignments_appended = vector("list", length=0)
    }

    # Append all data into a single list
    for (i in 1:no.of.blocks) {
      trees_appended = append(trees_appended, trees_all[[i]])
      node.assignments_appended = append(node.assignments_appended, node.assignments_all[[i]])
      if(!is.na(bin.size)) {
        binned.node.assignments_appended = append(binned.node.assignments_appended, binned.node.assignments_all[[i]])
      }
    }
    # Determine the order in which items should be
    idx <- order(c(unlist(lapply(trees_all,seq_along))))

    # Select the items in the right order and make sure the data remain strings
    trees_interleaved = trees_appended[idx]
    node.assignments_interleaved = do.call(cbind.data.frame,node.assignments_appended[idx])
    node.assignments_interleaved = data.frame(lapply(node.assignments_interleaved, as.character), stringsAsFactors=FALSE)
    
    if(!is.na(bin.size)) {
      binned.node.assignments_interleaved = do.call(cbind.data.frame,binned.node.assignments_appended[idx])
      binned.node.assignments_interleaved = data.frame(lapply(binned.node.assignments_interleaved, as.character), stringsAsFactors=FALSE)
    }

    # Sum the strengths for the different categories
    ancestor.strengths = as.matrix(Reduce('+', ancestor.strengths_all))
    sibling.strengths = as.matrix(Reduce('+', sibling.strengths_all))
    identity.strengths = as.matrix(Reduce('+', identity.strengths_all))
    parent.child.strengths = as.matrix(Reduce('+', parent.child.strengths_all))
    child.parent.strengths = as.matrix(Reduce('+', child.parent.strengths_all))
    strengths = list(ancestor.strengths=ancestor.strengths, sibling.strengths=sibling.strengths, identity.strengths=identity.strengths)

    # Save the matrices
    write.strengths.table(ancestor.strengths, removed_indices, paste("ancestor.strengths_",samplename,"_",no.iters,"iters_combined.txt",sep=""))
    write.strengths.table(sibling.strengths, removed_indices, paste("sibling.strengths_",samplename,"_",no.iters,"iters_combined.txt",sep=""))
    write.strengths.table(identity.strengths, removed_indices, paste("identity.strengths_",samplename,"_",no.iters,"iters_combined.txt",sep=""))
    write.strengths.table(parent.child.strengths, removed_indices, paste("parent.child.strengths_",samplename,"_",no.iters,"iters_combined.txt",sep=""))
    write.strengths.table(child.parent.strengths, removed_indices, paste("child.parent.strengths_",samplename,"_",no.iters,"iters_combined.txt",sep=""))

    # Adjust the no.iters and no.iters.burn.in to take into account the combined data
    no.iters.burn.in = no.iters.burn.in*no.of.blocks
    no.iters = no.iters*no.of.blocks

    # Print the combined trees
    pdf(paste("all_trees_",samplename,"_",no.iters,"iters.pdf",sep=""),height=4,width=4*no.subsamples)#ncol(mutCount)
    for(iter in 1:no.iters){
      tree = trees_interleaved[[iter]]
      tree$annotation = NA
      tree = annotateTree(tree,node.assignments_interleaved[,iter],annotation)
      title = paste(samplename,"iter",iter)
      plotTree(tree,title)
    }
    dev.off()
    
    # Clean up unused variables
    rm(trees_appended, node.assignments_appended, trees_all, node.assignments_all)
    rm(ancestor.strengths_all, sibling.strengths_all, identity.strengths_all, parent.child.strengths_all, child.parent.strengths_all)
    gc()
  }
  
  ###################################
  # Run the cons phase of the method
  ###################################
  if (is.na(phase) | phase == 'tree_dp' | phase == 'cons') {
    ###################################
    # Start Cons
    ###################################
    cons_start_time = Sys.time()
    if(!is.na(bin.size)) { 
      consResults = RunTreeBasedDPConsensus(trees=trees_interleaved, 
                                            node.assignments=binned.node.assignments_interleaved, 
                                            mutCount=binned.mutCount, 
                                            WTCount=binned.WTCount, 
                                            kappa=binned.kappa, 
                                            samplename=samplename, 
                                            subsamplenames=subsamplenames, 
                                            annotation=annotation, 
                                            no.iters=no.iters, 
                                            no.iters.burn.in=no.iters.burn.in, 
                                            resort.mutations=resort.mutations, 
                                            bin.indices=bin.indices,
                                            strengths=strengths)
    } else {
      consResults = RunTreeBasedDPConsensus(trees=trees_interleaved, 
                                            node.assignments=node.assignments_interleaved, 
                                            mutCount=mutCount, 
                                            WTCount=WTCount, 
                                            kappa=kappa, 
                                            samplename=samplename, 
                                            subsamplenames=subsamplenames, 
                                            annotation=annotation, 
                                            no.iters=no.iters, 
                                            no.iters.burn.in=no.iters.burn.in, 
                                            resort.mutations=resort.mutations, 
                                            bin.indices=NULL,
                                            strengths=strengths)
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

  # Return the consensus assignments when available (not after running the tree phase only)
  if (is.na(phase) | phase == 'tree_dp' | phase == 'cons') {
    return(consResults)
  } else {
    return(list())
  }

}


DirichletProcessClustering <- function(mutCount, WTCount, totalCopyNumber, copyNumberAdjustment, mutation.copy.number, cellularity, output_folder, no.iters, no.iters.burn.in, subsamplesrun, samplename, conc_param, cluster_conc, mut.assignment.type, most.similar.mut, mutationTypes) {
  #
  # Run the regular Dirichlet Process based method. Will perform clustering using the given data. The method
  # decides automatically whether the 1D or nD method is run based on the number of samples given at the input.
  # The number of samples is determined through the number of columns of the input.
  #
  # The nD method will yield a series of figures in which each sample is plotted against each other sample. 
  # The 1D method yields a density plot for just the single sample.
  #
  set.seed(123)
  
  print(dim(mutCount))
  if(!file.exists(output_folder)){
    dir.create(output_folder)
  }
  GS.data = subclone.dirichlet.gibbs(mutCount=mutCount,
                                    WTCount=WTCount,
                                    totalCopyNumber=totalCopyNumber, 
                                    copyNumberAdjustment=copyNumberAdjustment, 
                                    cellularity=cellularity,
                                    iter=no.iters,
                                    conc_param=conc_param,
                                    cluster_conc=cluster_conc)
  
  write.csv(GS.data$S.i,paste(output_folder,"/",samplename,"_2D_iters",no.iters,"_concParam",conc_param,"_clusterWidth",1/cluster_conc,"_states.csv",sep=""))
  write.csv(GS.data$V.h,paste(output_folder,"/",samplename,"_2D_iters",no.iters,"_concParam",conc_param,"_clusterWidth",1/cluster_conc,"_stickbreaking_weights.csv",sep=""))
  write.csv(GS.data$pi.h,paste(output_folder,"/",samplename,"_2D_iters",no.iters,"_concParam",conc_param,"_clusterWidth",1/cluster_conc,"_discreteMutationCopyNumbers.csv",sep=""))
  write.csv(GS.data$alpha,paste(output_folder,"/",samplename,"_2D_iters",no.iters,"_concParam",conc_param,"_clusterWidth",1/cluster_conc,"_alpha.csv",sep=""))
  
#      GS.data = list()
#      GS.data$S.i = as.matrix(read.csv(paste(output_folder,"/",samplename,"_2D_iters",no.iters,"_concParam",conc_param,"_clusterWidth",1/cluster_conc,"_states.csv",sep=""),row.names=1))
#      GS.data$V.h = as.matrix(read.csv(paste(output_folder,"/",samplename,"_2D_iters",no.iters,"_concParam",conc_param,"_clusterWidth",1/cluster_conc,"_stickbreaking_weights.csv",sep=""),row.names=1))
#      GS.data$pi.h = as.matrix(read.csv(paste(output_folder,"/",samplename,"_2D_iters",no.iters,"_concParam",conc_param,"_clusterWidth",1/cluster_conc,"_discreteMutationCopyNumbers.csv",sep=""),row.names=1))
#      GS.data$alpha = read.csv(paste(output_folder,"/",samplename,"_2D_iters",no.iters,"_concParam",conc_param,"_clusterWidth",1/cluster_conc,"_alpha.csv",sep=""),row.names=1)
# nD dataset, plot sample versus sample
if (ncol(mutCount) > 1) {
  ########################
  # Plot density and Assign mutations to clusters - nD
  ########################
  for(i in 1:(length(subsamplesrun)-1)){
    for(j in (i+1):length(subsamplesrun)){
      imageFile = paste(output_folder,"/",samplename,subsamplesrun[i],subsamplesrun[j],"_iters",no.iters,"_concParam",conc_param,"_clusterWidth",1/cluster_conc,"_2D_binomial.png",sep="")
      density = Gibbs.subclone.density.est(mutation.copy.number[,c(i,j)]/copyNumberAdjustment[,c(i,j)],
                                           GS.data,
                                           imageFile, 
                                           post.burn.in.start = no.iters.burn.in, 
                                           post.burn.in.stop = no.iters, 
                                           samplenames = paste(samplename,subsamplesrun[c(i,j)],sep=""),
                                           indices=c(i,j))  
      save(file=paste(output_folder,"/",samplename, subsamplesrun[i], subsamplesrun[j], "_densityoutput.RData", sep=""), GS.data, density)
    }
  }
  
  # Assign mutations to clusters using one of the different assignment methods
  opts = list(samplename=samplename, subsamplenames=subsamplesrun, no.iters=no.iters, no.iters.burn.in=no.iters.burn.in, no.iters.post.burn.in=no.iters-no.iters.burn.in, outdir=output_folder)
  print("Assigning mutations to clusters")
  if (mut.assignment.type == 1) {
    consClustering = multiDimensionalClustering(mutation.copy.number=mutation.copy.number, 
                                                copyNumberAdjustment=copyNumberAdjustment, 
                                                GS.data=GS.data, 
                                                density.smooth=0.01, 
                                                opts=opts)
  } else if (mut.assignment.type == 2) {
    consClustering = mutation_assignment_em(mutCount=mutCount, WTCount=WTCount, node.assignments=GS.data$S.i, opts=opts)
    
  } else if (mut.assignment.type == 3) {  
    warning("binom mut assignment not implemented for multiple timepoints")
    q(save="no", status=1)
    
  } else {
    warning(paste("Unknown mutation assignment type", mut.assignment.type, sep=" "))
    q(save="no", status=1)
  }
  return(consClustering)
    
    # 1D dataset, plot just the single density
  } else {
    ########################
    # Plot density and Assign mutations to clusters - 1D
    ########################
    # 1D dataset, plot just the single density
    wd = getwd()
    setwd(output_folder)
    res = Gibbs.subclone.density.est.1d(GS.data, 
                                            paste(samplename,"_DirichletProcessplot.png", sep=''), 
                                            samplename=samplename,
                                            post.burn.in.start=no.iters.burn.in, 
                                            post.burn.in.stop=no.iters,
                                            y.max=15,
					                                  x.max=NA, 
                                            mutationCopyNumber=mutation.copy.number, 
                                            no.chrs.bearing.mut=copyNumberAdjustment)
    density = res$density
    polygon.data = res$polygon.data
    
    # Assign mutations to clusters using one of the different assignment methods
    opts = list(samplename=samplename, subsamplenames=subsamplesrun, no.iters=no.iters, no.iters.burn.in=no.iters.burn.in, no.iters.post.burn.in=no.iters-no.iters.burn.in, outdir=output_folder)
    print("Assigning mutations to clusters")
    if (mut.assignment.type == 1) {
      subclonal.fraction = mutation.copy.number / copyNumberAdjustment
      subclonal.fraction[is.nan(subclonal.fraction)] = 0
      consClustering = oneDimensionalClustering(samplename, subclonal.fraction, GS.data, density, no.iters, no.iters.burn.in)
     
      setwd(wd) # Go back to original work directory 
      # Replot the data with cluster locations
      library(ggplot2)
      plot1D_2(density=density, 
             polygon.data=polygon.data, 
             pngFile=paste(output_folder, "/", samplename, "_DirichletProcessplot_with_cluster_locations.png", sep=""), 
             density.from=0, 
             x.max=1.5, 
             mutationCopyNumber=mutation.copy.number, 
             no.chrs.bearing.mut=copyNumberAdjustment,
             samplename=samplename,
             cluster.locations=consClustering$cluster.locations,
             mutation.assignments=consClustering$best.node.assignments,
             mutationTypes=mutationTypes)
      
      library(gridExtra)
      # Plot a table with the assignments
      plotAssignmentTable(cluster_locations=consClustering$cluster.locations, 
                          pngFile=paste(output_folder, "/", samplename, "_mutation_assignments.png", sep=""))
      
    } else if (mut.assignment.type == 2) {
      setwd(wd) # set the wd back earlier. The oneD clustering and gibbs sampler do not play nice yet and need the switch, the em assignment doesnt
      consClustering = mutation_assignment_em(mutCount=mutCount, WTCount=WTCount, node.assignments=GS.data$S.i, opts=opts)
      
    } else if (mut.assignment.type == 3) {  
      consClustering = mutation_assignment_binom(clustering_density=density,
                                                 mutCount=mutCount, 
                                                 WTCount=WTCount, 
                                                 copyNumberAdjustment=copyNumberAdjustment, 
                                                 tumourCopyNumber=totalCopyNumber,
                                                 normalCopyNumber=array(2, dim(mutCount)),
                                                 cellularity=cellularity)
      
      all.likelihoods = consClustering$all.likelihoods
      colnames(all.likelihoods) = paste("prob.cluster", 1:ncol(all.likelihoods))
      write.table(all.likelihoods, file=paste(output_folder, "/", samplename, "_DP_and_cluster_info.txt", sep=""), quote=F, row.names=F, sep="\t")
      
      setwd(wd) # Go back to original work directory 
      # Replot the data with cluster locations
      plot1D(density=density, 
             polygon.data=polygon.data, 
             pngFile=paste(output_folder, "/", samplename, "_DirichletProcessplot_with_cluster_locations.png", sep=""), 
             density.from=0, 
             x.max=1.5, 
             mutationCopyNumber=mutation.copy.number, 
             no.chrs.bearing.mut=copyNumberAdjustment,
             samplename=samplename,
             cluster.locations=consClustering$cluster.locations,
             mutation.assignments=consClustering$best.node.assignments)
      
    } else {
      warning(paste("Unknown mutation assignment type", mut.assignment.type, sep=" "))
      q(save="no", status=1)
    }
    setwd(wd)
    
    return(consClustering)
  }
  
}

write.strengths.table = function(dat, removed_indices, filename) {
  #
  # Adds in an empty column/row for each of the removed_indices and writes
  # the subsequent matrix to disk
  #
  write.table(add.muts.back.in(dat, removed_indices),filename,sep="\t",row.names=F,quote=F,col.names=F)
}

add.muts.back.in = function(dat, removed_indices, def.value=0) {
  #
  # Adds in empty columns and rows for mutations that were removed.
  # Mutations are added sequentially, so this method expects indices
  # of removed mutations in the original (full) matrix. The empty
  # mutations will be added in the place where they were removed,
  # keeping the order in tact.
  #
  for (i in removed_indices) { 
    if (i >= ncol(dat)) {
      	dat = cbind(dat, rep(def.value, nrow(dat)))
      	dat = rbind(dat, rep(def.value, ncol(dat)))
    } else {
        dat = cbind(dat[,1:(i-1)], rep(def.value, nrow(dat)), dat[,i:ncol(dat)])
        dat = rbind(dat[1:(i-1),], rep(def.value, ncol(dat)), dat[i:nrow(dat),])
    }
  }
  return(dat)
}

