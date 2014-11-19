source("DirichletProcessClustering.R")

RunDP <- function(analysis_type, dataset, samplename, subsamples, no.iters, no.iters.burn.in, outdir, conc_param, cluster_conc, resort.mutations, parallel, blockid, no.of.blocks, annotation=vector(mode="character",length=nrow(dataset$mutCount)), init.alpha=0.01, shrinkage.threshold=0.1, remove.node.frequency=NA, remove.branch.frequency=NA, bin.size=NA) {
  # Pick the analysis to run
  if (analysis_type == 'nd_dp') {
    clustering = DirichletProcessClustering(mutCount=dataset$mutCount, 
                                            WTCount=dataset$WTCount, 
                                            no.iters=no.iters, 
                                            no.iters.burn.in=no.iters.burn.in, 
                                            cellularity=cellularity, 
                                            totalCopyNumber=dataset$totalCopyNumber, 
                                            mutation.copy.number=dataset$mutation.copy.number,
                                            copyNumberAdjustment=dataset$copyNumberAdjustment, 
                                            samplename=samplename, 
                    			                  subsamplesrun=subsamples,
                                            output_folder=outdir, 
                                            conc_param=conc_param, 
                                            cluster_conc=cluster_conc)
  } else if(analysis_type == "tree_dp" | analysis_type == 'tree' | analysis_type == 'cons') {
    if(is.na(bin.size)){
      outdir = paste(samplename,"_treeBasedDirichletProcessOutputs_noIters",no.iters,"_burnin",no.iters.burn.in,sep="")
      
    }else{
      outdir = paste(samplename,"_treeBasedDirichletProcessOutputs_noIters",no.iters,"_binSize",bin.size,"_burnin",no.iters.burn.in,sep="")
    }
    clustering = TreeBasedDP(mutCount=dataset$mutCount,
                             WTCount=dataset$WTCount,
                             kappa=dataset$kappa, 
                             samplename=samplename, 
                             subsamplenames=subsamples,
                             no.iters=no.iters,
                             no.iters.burn.in=no.iters.burn.in,
                             cellularity=cellularity,
                             bin.size=bin.size, 
                             resort.mutations=resort.mutations, 
                             outdir=outdir, 
                             parallel=parallel, 
                             phase=analysis_type, 
                             blockid=blockid, 
                             no.of.blocks=no.of.blocks,
                             remove.node.frequency=remove.node.frequency, 
                             remove.branch.frequency=remove.branch.frequency,
                             annotation=annotation,
                             init.alpha=init.alpha, 
                             shrinkage.threshold=shrinkage.threshold)
  } else {
    print(paste("Unknown type of analysis",analysis_type))
    q(save="no", status=1)
  }

  if (analysis_type != 'tree') {
    # Write final output
    outfiles.prefix = paste(samplename, "_", no.iters, "iters_", no.iters.burn.in, sep="")
    output = cbind(dataset$chromosome[,1], dataset$position[,1]-1, dataset$position[,1], clustering$best.node.assignments, clustering$best.assignment.likelihoods)
    
    # Add the removed mutations back in
    for (i in dataset$removed_indices) {
      output = rbind(output[1:(i-1),], c(dataset$chromosome.not.filtered[i], dataset$mut.position.not.filtered[i]-1, dataset$mut.position.not.filtered[i], NA, NA), output[i:nrow(output),])
    }
    
    # Save the indices of the mutations that were not used during the analysis
    write.table(data.frame(mut.index=dataset$removed_indices), file=paste(outfiles.prefix,"removedMutationsIndex.txt", sep=""))

    # Save the consensus mutation assignments
    save(file=paste(outfiles.prefix, "burnin_bestConsensusResults.RData", sep=""), output, clustering, samplename, outdir, no.iters, no.iters.burn.in)
    colnames(output) = c("chr", "start", "end", "cluster", "likelihood")
    write.table(output, file=paste(outfiles.prefix, "burnin_bestConsensusAssignments.bed", sep=""), quote=F, row.names=F, sep="\t")

    # If tree based analysis, also save the tree
    if (analysis_type != 'nd_dp') {
      write.table(clustering$best.tree, file=paste(outfiles.prefix, "burnin_bestConsensusTree.txt", sep=""), quote=F, row.names=F, sep="\t")
    }
  }
}
