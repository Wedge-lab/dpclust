# source("DirichletProcessClustering.R")
# source("PlotDensities.R")

RunDP <- function(analysis_type, dataset, samplename, subsamples, no.iters, no.iters.burn.in, outdir, conc_param, cluster_conc, resort.mutations, parallel, blockid, no.of.blocks, mut.assignment.type, annotation=vector(mode="character",length=nrow(dataset$mutCount)), init.alpha=0.01, shrinkage.threshold=0.1, remove.node.frequency=NA, remove.branch.frequency=NA, bin.size=NA, muts.sampled=F) {
  # Obtain the mutations that were not sampled, as these must be assigned to clusters separately
  if (muts.sampled) {
    most.similar.mut = dataset$most.similar.mut
  } else {
    most.similar.mut = NA
  }
  
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
                                            cluster_conc=cluster_conc,
                    			                  mut.assignment.type=mut.assignment.type,
							                              most.similar.mut=most.similar.mut)
    
  } else if (analysis_type == "tree_dp" | analysis_type == 'tree' | analysis_type == 'cons') {
	# REMOVE temp CNA branching testing
# 	mutCount = dataset$mutCount
#   	WTCount = dataset$WTCount
# 	kappa = dataset$kappa
# 	conflict_indices = c(which(dataset$position==20734478), which(dataset$position==24377093), which(dataset$position==32866944))

    clustering = TreeBasedDP(mutCount=dataset$mutCount,
                             WTCount=dataset$WTCount,
                             removed_indices=dataset$removed_indices,
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
                             shrinkage.threshold=shrinkage.threshold,
			     conflict.array=array(1,c(nrow(dataset$mutCount),nrow(dataset$mutCount))))

  } else if (analysis_type == "replot_1d") {
    ##############################
    # Replot 1D clustering
    ##############################
    density = read.table(paste(outdir, "/", samplename, "_DPoutput_", no.iters, "iters_", no.iters.burn.in, "burnin", "/", samplename, "_DirichletProcessplotdensity.txt", sep=""), header=T)
    polygon.data = read.table(paste(outdir, "/", samplename, "_DPoutput_", no.iters, "iters_", no.iters.burn.in, "burnin", "/", samplename, "_DirichletProcessplotpolygonData.txt", sep=""), header=T)
    pngFile = paste(outdir, "/", samplename, "_DPoutput_", no.iters, "iters_", no.iters.burn.in, "burnin", "/", samplename, "_DirichletProcessplot_replot.png", sep="")
    
    plot1D(density, 
           polygon.data[,1], 
           pngFile=pngFile, 
           density.from=0, 
           #y.max=6, 
           x.max=1.5, #2.7
           mutationCopyNumber=dataset$mutation.copy.number, 
           no.chrs.bearing.mut=dataset$copyNumberAdjustment,
           samplename=samplename)
    
  } else if (analysis_type == "replot_nd") {  
    ##############################
    # Replot nD clustering
    ##############################
    for(i in 1:(length(subsamples)-1)){
      for(j in (i+1):length(subsamples)){
        filename_prefix = paste(outdir,"/",samplename, "_DPoutput_", no.iters, "iters_", no.iters.burn.in, "burnin/", samplename,subsamples[i],subsamples[j],"_iters",no.iters,"_concParam",conc_param,"_clusterWidth",1/cluster_conc,sep="")
        pngFile = paste(filename_prefix,"_2D_binomial_limitedRange_replot.png",sep="")
        xvals = read.table(paste(filename_prefix,"_2D_binomial_xvals.csv",sep=""),header=T,sep=",",row.names=1)
        yvals = read.table(paste(filename_prefix,"_2D_binomial_yvals.csv",sep=""),header=T,sep=",",row.names=1)
        zvals = read.table(paste(filename_prefix,"_2D_binomial_zvals.csv",sep=""),header=T,sep=",",row.names=1)
        
        #
        # TODO
        # subclonal.fraction directly is not right for replotting the Myeloma samples :: check!
        #
        
        plotnD(xvals=xvals, 
               yvals=yvals, 
               zvals=zvals, 
               subclonal.fraction_x=dataset$subclonal.fraction[,i], 
               subclonal.fraction_y=dataset$subclonal.fraction[,j], 
               pngFile=pngFile, 
               samplename_x=paste(samplename,subsamples[i], sep=""), 
               samplename_y=paste(samplename,subsamples[j], sep=""), 
               max.plotted.value=1.4)
      }
    }
    
  } else if (analysis_type == "reassign_muts_1d") {
    ##############################
    # Reassign mutations to clusters using a previous 1D clustering run
    ##############################
    GS.data = list()
    full_outdir = paste(outdir, "/", samplename, "_DPoutput_", no.iters, "iters_", no.iters.burn.in, "burnin", "/", sep="")
    filename_prefix = paste(full_outdir, samplename, "_2D_iters", no.iters, "_concParam", conc_param, "_clusterWidth", 1/cluster_conc, sep="")
    GS.data$S.i = as.matrix(read.csv(paste(filename_prefix, "_states.csv", sep=""), row.names=1))
    GS.data$V.h = as.matrix(read.csv(paste(filename_prefix, "_stickbreaking_weights.csv", sep=""), row.names=1))
    GS.data$pi.h = as.matrix(read.csv(paste(filename_prefix, "_discreteMutationCopyNumbers.csv", sep=""), row.names=1))
    GS.data$alpha = read.csv(paste(filename_prefix, "_alpha.csv", sep=""), row.names=1)
    
    # Fix the dimensions of these arrays
    S.i = array(NA, c(nrow(GS.data$S.i), ncol(GS.data$S.i), 1))
    S.i[1:nrow(GS.data$S.i), 1:ncol(GS.data$S.i), 1] = GS.data$S.i
    GS.data$S.i = S.i
    
    pi.h = array(NA, c(nrow(GS.data$pi.h), ncol(GS.data$pi.h), 1))
    pi.h[1:nrow(GS.data$pi.h), 1:ncol(GS.data$pi.h), 1] = GS.data$pi.h
    GS.data$pi.h = pi.h
    
    # TODO: REMOVE - temp overwrite of no.iters.burn.in
#     no.iters.burn.in.supplied = no.iters.burn.in
#     no.iters.burn.in = no.iters-1000 # Taking only the latest 500 iterations
#     print(paste("BURNIN", no.iters.burn.in))
    
    opts = list(samplename=samplename, subsamplenames=subsamples, no.iters=no.iters, no.iters.burn.in=no.iters.burn.in, no.iters.post.burn.in=no.iters-no.iters.burn.in, outdir=outdir)
    
    #clustering = mutation_assignment_em(mutCount=dataset$mutCount, WTCount=dataset$WTCount, node.assignments=GS.data$S.i, opts=opts)
    
    if (mut.assignment.type == 1) {
      wd = getwd()
      setwd(paste(full_outdir, sep=""))
      res = Gibbs.subclone.density.est.1d(GS.data, 
                                          paste(samplename,"_DirichletProcessplot.png", sep=''), 
                                          samplename=samplename,
                                          post.burn.in.start=no.iters.burn.in, 
                                          post.burn.in.stop=no.iters,
                                          y.max=15,
                                          x.max=NA, 
                                          mutationCopyNumber=dataset$mutation.copy.number, 
                                          no.chrs.bearing.mut=dataset$copyNumberAdjustment)
      density = res$density
      polygon.data = res$polygon.data
      
      
      #density = read.table(paste(full_outdir, samplename, "_DirichletProcessplotdensity.txt", sep=""), header=T)
      subclonal.fraction = dataset$mutation.copy.number / dataset$copyNumberAdjustment
      subclonal.fraction[is.nan(subclonal.fraction)] = 0
      clustering = oneDimensionalClustering(samplename, subclonal.fraction, GS.data, density, no.iters, no.iters.burn.in)
      
      # Disabled for now
      setwd(wd)
#       # Replot the data with cluster locations
#       plot1D(density=density, 
#              polygon.data=polygon.data, 
#              pngFile=paste(full_outdir, "/", samplename, "_DirichletProcessplot_with_cluster_locations.png", sep=""), 
#              density.from=0, 
#              x.max=1.5, 
#              mutationCopyNumber=dataset$mutation.copy.number, 
#              no.chrs.bearing.mut=dataset$copyNumberAdjustment,
#              samplename=samplename,
#              cluster.locations=clustering$cluster.locations,
#              mutation.assignments=clustering$best.node.assignments)
      
    } else if (mut.assignment.type == 2) {
      setwd(wd)
      clustering = mutation_assignment_em(mutCount=dataset$mutCount, WTCount=dataset$WTCount, node.assignments=GS.data$S.i, opts=opts)
      
    } else {
      print(paste("Unknown mutation assignment type", mut.assignment.type, sep=" "))
      q(save="no", status=1)
    }
#     # TODO: REMOVE - temp overwrite of no.iters.burn.in
#     no.iters.burn.in = no.iters.burn.in.supplied
    
    
  } else {
    print(paste("Unknown type of analysis",analysis_type))
    q(save="no", status=1)
  }

  if (all(!analysis_type %in% c('tree', 'replot_1d', 'replot_nd', 'reassign_muts_1d'))) {

	  # Check if mutation sampling has been done, if so, unpack and assign here
  	if (!is.na(most.similar.mut)) {
  	  res = unsample_mutations(dataset, clustering_result)
      dataset = res$dataset
      clustering_result = res$clustering_result
  	}

    # Write final output
    outfiles.prefix = paste(outdir, "/", samplename, "_", no.iters, "iters_", no.iters.burn.in, "burnin", sep="")
    output = cbind(dataset$chromosome[,1], dataset$position[,1]-1, dataset$position[,1], clustering$best.node.assignments, clustering$best.assignment.likelihoods)
    
    # Add the removed mutations back in
    for (i in dataset$removed_indices) {
      if (i==1) {
        output = rbind(c(dataset$chromosome.not.filtered[i], dataset$mut.position.not.filtered[i]-1, dataset$mut.position.not.filtered[i], NA, NA), output)
      } else if (i >= nrow(output)) {
        output = rbind(output, c(dataset$chromosome.not.filtered[i], dataset$mut.position.not.filtered[i]-1, dataset$mut.position.not.filtered[i], NA, NA))
      } else {
        output = rbind(output[1:(i-1),], c(dataset$chromosome.not.filtered[i], dataset$mut.position.not.filtered[i]-1, dataset$mut.position.not.filtered[i], NA, NA), output[i:nrow(output),])
      }
    }
    
    # Save the indices of the mutations that were not used during the analysis
    write.table(data.frame(mut.index=dataset$removed_indices), file=paste(outfiles.prefix,"_removedMutationsIndex.txt", sep=""), row.names=F, quote=F)

    # Save the consensus mutation assignments
    save(file=paste(outfiles.prefix, "_bestConsensusResults.RData", sep=""), output, clustering, samplename, outdir, no.iters, no.iters.burn.in)
    print(head(output))
    colnames(output) = c("chr", "start", "end", "cluster", "likelihood")
    write.table(output, file=paste(outfiles.prefix, "_bestConsensusAssignments.bed", sep=""), quote=F, row.names=F, sep="\t")

    # If tree based analysis, also save the tree
    if (analysis_type != 'nd_dp') {
      write.table(clustering$best.tree, file=paste(outfiles.prefix, "_bestConsensusTree.txt", sep=""), quote=F, row.names=F, sep="\t")
    }
  }
}
