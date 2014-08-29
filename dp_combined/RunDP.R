source("DirichletProcessClustering.R")

RunDP <- function(analysis_type, dataset, samplename, subsamples, no.iters, no.iters.burn.in, outdir, conc_param, cluster_conc, resort.mutations, parallel, blockid, no.of.blocks, annotation=vector(mode="character",length=nrow(dataset$mutCount)), init.alpha=0.01, shrinkage.threshold=0.1, remove.node.frequency=NA, remove.branch.frequency=NA)
  # Pick the analysis to run
  if (analysis_type == 'nd_dp') {
    DirichletProcessClustering(mutCount=dataset$mutCount, 
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
    TreeBasedDP(mutCount=dataset$mutCount,
                   WTCount=dataset$WTCount,
                   kappa=dataset$kappa, 
                   samplename = samplename, 
                   subsamplenames = subsamples,
                   no.iters=no.iters,
                   no.iters.burn.in=no.iters.burn.in,
                   cellularity=cellularity,
                   bin.size = bin.size, 
                   resort.mutations = resort.mutations, 
                   outdir = outdir, 
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
  }
