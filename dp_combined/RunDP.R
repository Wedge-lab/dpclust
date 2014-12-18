source("DirichletProcessClustering.R")
source("PlotDensities.R")

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
                             shrinkage.threshold=shrinkage.threshold)
    
  } else if(analysis_type == "replot_1d") {
    ##############################
    # Replot 1D clustering
    ##############################
    density = read.table(paste(outdir, "/", samplename, "_DPoutput_1250iters_250burnin", "/", samplename, "_DirichletProcessplotdensity.txt", sep=""), header=T)
    polygon.data = read.table(paste(outdir, "/", samplename, "_DPoutput_1250iters_250burnin", "/", samplename, "_DirichletProcessplotpolygonData.txt", sep=""), header=T)
    pngFile = paste(outdir, "/", samplename, "_DPoutput_1250iters_250burnin", "/", samplename, "_DirichletProcessplot_replot.png", sep="")
    
    plot1D(density, 
           polygon.data, 
           pngFile=pngFile, 
           density.from=0, 
           #y.max=6, 
           x.max=2, 
           mutationCopyNumber=dataset$mutation.copy.number, 
           no.chrs.bearing.mut=dataset$copyNumberAdjustment,
           samplename=samplename)
    
  } else if(analysis_type == "replot_nd") {  
    ##############################
    # Replot nD clustering
    ##############################
    for(i in 1:(length(subsamples)-1)){
      for(j in (i+1):length(subsamples)){
        filename_prefix = paste(outdir,"/",samplename, "_DPoutput_25iters_2burnin/",samplename,subsamples[i],subsamples[j],"_iters",no.iters,"_concParam",conc_param,"_clusterWidth",1/cluster_conc,sep="")
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
    
  } else {
    print(paste("Unknown type of analysis",analysis_type))
    q(save="no", status=1)
  }

  if (analysis_type != 'tree' & analysis_type != 'replot_1d' & analysis_type != 'replot_nd') {
    # Write final output
    outfiles.prefix = paste(samplename, "_", no.iters, "iters_", no.iters.burn.in, "burnin", sep="")
    output = cbind(dataset$chromosome[,1], dataset$position[,1]-1, dataset$position[,1], clustering$best.node.assignments, clustering$best.assignment.likelihoods)
    
    # Add the removed mutations back in
    for (i in dataset$removed_indices) {
      if (i >= nrow(output)) {
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
