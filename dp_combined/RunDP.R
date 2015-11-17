# source("DirichletProcessClustering.R")
# source("PlotDensities.R")

RunDP <- function(analysis_type, dataset, samplename, subsamples, no.iters, no.iters.burn.in, outdir, conc_param, cluster_conc, resort.mutations, parallel, blockid, no.of.blocks, mut.assignment.type, annotation=vector(mode="character",length=nrow(dataset$mutCount)), init.alpha=0.01, shrinkage.threshold=0.1, remove.node.frequency=NA, remove.branch.frequency=NA, bin.size=NA, num_muts_sample=NA, cndata=NULL, add.conflicts=F, cna.conflicting.events.only=F, sample.snvs.only=F) {
  # Check if co-clustering of copy number data is in order
  if (!is.null(cndata)) {
    dataset = add.in.cn.as.snv.cluster(dataset, cndata, add.conflicts=add.conflicts, conflicting.events.only=cna.conflicting.events.only)
  } else if (!is.null(dataset$cndata)) {
    # In case of a rerun, pull out the cndata
    cndata = dataset$cndata
  }
  
  # Obtain the mutations that were not sampled, as these must be assigned to clusters separately
  if (!is.na(num_muts_sample) & num_muts_sample!="NA") {
    dataset = sample_mutations(dataset, num_muts_sample, sample.snvs.only=sample.snvs.only)
    most.similar.mut = dataset$most.similar.mut
  } else {
    most.similar.mut = NA
  }
  dataset$cndata = cndata
  save(file=paste(outdir, "/dataset.RData", sep=""), dataset)

  # Pick the analysis to run
  if (analysis_type == 'sample_muts') {
    # If only sampling then quit now. Use this when running parts of a method in parallel
    q(save="no")
  } else if (analysis_type == 'nd_dp') {
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
							                              most.similar.mut=most.similar.mut,
							                              mutationTypes=dataset$mutationType)
    
  } else if (analysis_type == "tree_dp" | analysis_type == 'tree' | analysis_type == 'cons') {

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
                             conflict.array=dataset$conflict.array)

  } else if (analysis_type == "replot_1d") {
    ##############################
    # Replot 1D clustering
    ##############################
    density = read.table(paste(outdir, "/", samplename, "_DirichletProcessplotdensity.txt", sep=""), header=T)
    polygon.data = read.table(paste(outdir, "/", samplename, "_DirichletProcessplotpolygonData.txt", sep=""), header=T)
    pngFile = paste(outdir, "/", samplename, "_DirichletProcessplot_replot.png", sep="")
    load(paste(outdir, "/", samplename, "_1250iters_250burnin_bestConsensusResults.RData", sep=""))
    cluster_locations = read.table(paste(outdir, "/", samplename, "_optimaInfo.txt", sep=""), header=T, stringsAsFactors=F)
    
    plot1D(density=density, 
           polygon.data=polygon.data[,1], 
           pngFile=pngFile, 
           density.from=0, 
           #y.max=6, 
           x.max=1.5, #2.7
           mutationCopyNumber=dataset$mutation.copy.number, 
           no.chrs.bearing.mut=dataset$copyNumberAdjustment,
           samplename=samplename,
           mutationTypes=dataset$mutationType,
           cluster.locations=cluster_locations,
           mutation.assignments=clustering$best.node.assignments)
    
  } else if (analysis_type == "replot_nd") {  
    ##############################
    # Replot nD clustering
    ##############################
    for(i in 1:(length(subsamples)-1)){
      for(j in (i+1):length(subsamples)){
        filename_prefix = paste(outdir,"/",samplename,subsamples[i],subsamples[j],"_iters",no.iters,"_concParam",conc_param,"_clusterWidth",1/cluster_conc,sep="")
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
    outfiles.prefix = paste(outdir, "/", samplename, "_", no.iters, "iters_", no.iters.burn.in, "burnin", sep="")
    # Dump a checkpoint
    save(file=paste(outfiles.prefix, "_bestConsensusResults.RData", sep=""), clustering, samplename, outdir, no.iters, no.iters.burn.in, dataset, cndata)
    
	  # Check if mutation sampling has been done, if so, unpack and assign here
  	if (!is.na(most.similar.mut)) {
  	 res = unsample_mutations(dataset, clustering)
      dataset = res$dataset
  	  clustering = res$clustering
  	}
    
    # If cndata was added, assign individual CNA events to clusters, instead of the pseudo-SNVs
    if (!is.null(cndata)) {
      cna_assignments = array(NA, c(nrow(cndata), 2))
      for (i in 1:nrow(cndata)) {
        is_overlapping_snv = dataset$chromosome==cndata[i,]$chr & dataset$position==cndata[i,]$startpos
        if (any(is_overlapping_snv)) {
          overlapping_snv_index = which(is_overlapping_snv)
          assignment_counts = table(clustering$best.node.assignments[overlapping_snv_index])
          cna_assignments[i,1] = as.numeric(names(assignment_counts)[which.max(assignment_counts)])
          cna_assignments[i,2] = mean(clustering$best.assignment.likelihoods[overlapping_snv_index])
        }
      }
      colnames(cna_assignments) = c("most.likely.cluster", "assignment.likelyhood")

      # Remove the pseudo SNVs from the clustering output
      pseudo_snv_index = which(dataset$mutationType=="CNA")
      clustering$best.node.assignments = clustering$best.node.assignments[-pseudo_snv_index]
      clustering$best.assignment.likelihoods = clustering$best.assignment.likelihoods[-pseudo_snv_index]

      # Replace the pseudo SNV clusters with the assigned CNA events
      dataset = remove_pseudo_snv_cna_clusters(dataset)
      cna_assignments = data.frame(cndata[,c("chr", "startpos", "endpos", "CNA")], cna_assignments)
      write.table(cna_assignments, file=paste(outfiles.prefix,"_bestConsensusCNAassignments.bed", sep=""), row.names=F, quote=F, sep="\t")
    }

    # Write final output of all mutation assignments
    output = cbind(dataset$chromosome[,1], dataset$position[,1]-1, dataset$position[,1], clustering$best.node.assignments, clustering$best.assignment.likelihoods, dataset$mutationType)

    # Add the removed mutations back in - Assuming here that only SNVs have been removed
    for (i in dataset$removed_indices) {
      if (i==1) {
        output = rbind(c(dataset$chromosome.not.filtered[i], dataset$mut.position.not.filtered[i]-1, dataset$mut.position.not.filtered[i], NA, NA, "SNV"), output)
      } else if (i >= nrow(output)) {
        output = rbind(output, c(dataset$chromosome.not.filtered[i], dataset$mut.position.not.filtered[i]-1, dataset$mut.position.not.filtered[i], NA, NA, "SNV"))
      } else {
        output = rbind(output[1:(i-1),], c(dataset$chromosome.not.filtered[i], dataset$mut.position.not.filtered[i]-1, dataset$mut.position.not.filtered[i], NA, NA, "SNV"), output[i:nrow(output),])
      }
    }
    
    # Save the indices of the mutations that were not used during the analysis
    write.table(data.frame(mut.index=dataset$removed_indices), file=paste(outfiles.prefix,"_removedMutationsIndex.txt", sep=""), row.names=F, quote=F)

    # Save the consensus mutation assignments
    print(head(output))
    colnames(output) = c("chr", "start", "end", "cluster", "likelihood", "mut_type")
    write.table(output, file=paste(outfiles.prefix, "_bestConsensusAssignments.bed", sep=""), quote=F, row.names=F, sep="\t")

    # If tree based analysis, also save the tree
    if (analysis_type != 'nd_dp') {
      write.table(clustering$best.tree, file=paste(outfiles.prefix, "_bestConsensusTree.txt", sep=""), quote=F, row.names=F, sep="\t")
    }
  }
}
