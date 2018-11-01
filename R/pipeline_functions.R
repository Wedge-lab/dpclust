#
# Convenience functions for the pipelie
#

replot_1D = function(outdir, outfiles.prefix, samplename, dataset, clustering, density, polygon.data) {
  # Old plot
  plot1D(density=density, 
         polygon.data=polygon.data[,1], 
         pngFile=paste(outdir, "/", samplename, "_DirichletProcessplot_with_cluster_locations_replot.png", sep=""), 
         density.from=0, 
         x.max=1.5, 
         mutationCopyNumber=dataset$mutation.copy.number, 
         no.chrs.bearing.mut=dataset$copyNumberAdjustment,
         samplename=samplename,
         cluster.locations=clustering$cluster.locations,
         mutation.assignments=clustering$best.node.assignments)
  
  # New plot
  plot1D_2(density=density, 
           polygon.data=polygon.data[,1],
           pngFile=paste(outdir, "/", samplename, "_DirichletProcessplot_with_cluster_locations_2_replot.png", sep=""), 
           density.from=0, 
           x.max=1.5, 
           mutationCopyNumber=dataset$mutation.copy.number, 
           no.chrs.bearing.mut=dataset$copyNumberAdjustment,
           samplename=samplename,
           cluster.locations=clustering$cluster.locations,
           mutation.assignments=clustering$best.node.assignments,
           mutationTypes=dataset$mutationType)
}

reassign_1D = function(outdir, samplename, no.iters, no.iters.burn.in, dataset, conc_param, cluster_conc, min.frac.snvs.cluster) {
  full_outdir = paste(outdir, "/", sep="")
  load(file=paste(full_outdir, samplename, "_gsdata.RData", sep=""))
  
  # Obtain the density
  wd = getwd()
  setwd(full_outdir)
  res = Gibbs.subclone.density.est.1d(GS.data, 
                                      paste(samplename,"_DirichletProcessplot.png", sep=''), 
                                      samplename=samplename,
                                      post.burn.in.start=no.iters.burn.in, 
                                      post.burn.in.stop=no.iters,
                                      y.max=15,
                                      x.max=NA, 
                                      mutationCopyNumber=dataset$mutation.copy.number, 
                                      no.chrs.bearing.mut=dataset$copyNumberAdjustment)
  clustering_density = res$density
  polygon.data = res$polygon.data
  opts = list(samplename=samplename, subsamplenames=subsamples, no.iters=no.iters, no.iters.burn.in=no.iters.burn.in, no.iters.post.burn.in=no.iters-no.iters.burn.in, outdir=outdir)
  
  # Use one of the three ways of assigning mutations
  if (mut.assignment.type == 1) {
    subclonal.fraction = dataset$mutation.copy.number / dataset$copyNumberAdjustment
    subclonal.fraction[is.nan(subclonal.fraction)] = 0
    clustering = oneDimensionalClustering(samplename, subclonal.fraction, GS.data, clustering_density, no.iters, no.iters.burn.in)
    setwd(wd)
    # Change the outfiles prefix to be able to identify the output files
    outfiles.prefix = paste(full_outdir, samplename, "_reassign_option_1", sep="")
    
  } else if (mut.assignment.type == 2) {
    setwd(wd)
    clustering = mutation_assignment_em(mutCount=dataset$mutCount, WTCount=dataset$WTCount, node.assignments=GS.data$S.i, opts=opts)
    # Change the outfiles prefix to be able to identify the output files
    outfiles.prefix = paste(full_outdir, samplename, "_reassign_option_2", sep="")
    
  } else if (mut.assignment.type == 3) {
    clustering = mutation_assignment_binom(clustering_density=clustering_density,
                                           mutCount=dataset$mutCount, 
                                           WTCount=dataset$WTCount, 
                                           copyNumberAdjustment=dataset$copyNumberAdjustment, 
                                           tumourCopyNumber=dataset$totalCopyNumber,
                                           normalCopyNumber=array(2, dim(dataset$mutCount)),
                                           cellularity=cellularity)
    setwd(wd) # Move back to work dir, steps after this can handle being in a different directory
    # Change the outfiles prefix to be able to identify the output files
    outfiles.prefix = paste(full_outdir, samplename, "_reassign_option_3", sep="")
    
  } else if (mut.assignment.type == 4) {
    clustering = mutation_assignment_mpear(GS.data=GS.data, 
                                               no.iters=no.iters, 
                                               no.iters.burn.in=no.iters.burn.in, 
                                               min.frac.snvs.cluster=min.frac.snvs.cluster, 
                                               dataset=dataset, 
                                               samplename=samplename, 
                                               outdir=full_outdir)
    setwd(wd) # Move back to work dir, steps after this can handle being in a different directory
    outfiles.prefix = paste(full_outdir, samplename, "_reassign_option_4", sep="")
    
  } else {
    print(paste("Unknown mutation assignment type", mut.assignment.type, sep=" "))
    q(save="no", status=1)
  }
  
  # Replot the data with cluster locations
  plot1D_2(density=clustering_density, 
           polygon.data=polygon.data, 
           pngFile=paste(outfiles.prefix, "_DirichletProcessplot_with_cluster_locations_2.png", sep=""), 
           density.from=0, 
           x.max=1.5, 
           mutationCopyNumber=dataset$mutation.copy.number, 
           no.chrs.bearing.mut=dataset$copyNumberAdjustment,
           mutationTypes=dataset$mutationType,
           samplename=samplename,
           cluster.locations=clustering$cluster.locations,
           mutation.assignments=clustering$best.node.assignments)
    
  return(list(clustering=clustering, outfiles.prefix=outfiles.prefix))
}


replot_nD = function(outdir, outfiles.prefix, samplename, subsamples, dataset, clustering, conc_param, cluster_conc) {
  for(i in 1:(length(subsamples)-1)){
    for(j in (i+1):length(subsamples)){
      filename_prefix = paste(outdir, "/", samplename,subsamples[i],subsamples[j],"_iters",no.iters,"_concParam",conc_param,"_clusterWidth",1/cluster_conc,sep="")
      pngFile = paste(filename_prefix, "_2D_binomial_with_cluster_locations_replot.png", sep="")
      xvals = read.table(paste(filename_prefix, "_2D_binomial_xvals.csv",sep=""), header=T, sep=",", row.names=1)
      yvals = read.table(paste(filename_prefix, "_2D_binomial_yvals.csv",sep=""), header=T, sep=",", row.names=1)
      zvals = read.table(paste(filename_prefix, "_2D_binomial_zvals.csv",sep=""), header=T, sep=",", row.names=1)
      
      plotnD(xvals=xvals, 
             yvals=yvals, 
             zvals=zvals, 
             subclonal.fraction_x=dataset$subclonal.fraction[,i], 
             subclonal.fraction_y=dataset$subclonal.fraction[,j], 
             pngFile=pngFile, 
             samplename_x=paste(samplename,subsamples[i], sep=""), 
             samplename_y=paste(samplename,subsamples[j], sep=""), 
             max.plotted.value=1.4,
             cluster.locations=clustering$cluster.locations[,c(i+1,j+1)]) # +1 because the first column contains the cluster numbers
    }
  }
}

reassign_nd = function(outdir, samplename, subsamples, no.iters, no.iters.burn.in, dataset, conc_param, cluster_conc) {
  # Load the output from the algorithm
  # GS.data = read_gsdata_object(outdir, no.iters=no.iters, conc_param=conc_param, cluster_conc=cluster_conc)
  load(file=file.path(outdir, paste(samplename, "_gsdata.RData", sep="")))
  
  # Assign mutations to clusters using one of the different assignment methods
  opts = list(samplename=samplename, subsamplenames=subsamples, no.iters=no.iters, no.iters.burn.in=no.iters.burn.in, no.iters.post.burn.in=no.iters-no.iters.burn.in, outdir=outdir)
  print("Assigning mutations to clusters")
  if (mut.assignment.type == 1) {
    clustering = multiDimensionalClustering(mutation.copy.number=dataset$mutation.copy.number, 
                                            copyNumberAdjustment=dataset$copyNumberAdjustment, 
                                            GS.data=GS.data, 
                                            density.smooth=0.01, 
                                            opts=opts)
    # Change the outfiles prefix to be able to identify the output files
    outfiles.prefix = file.path(outdir, paste(samplename, "_reassign_option_1", sep=""))
    
  } else if (mut.assignment.type == 2) {
    clustering = mutation_assignment_em(mutCount=dataset$mutCount, WTCount=dataset$WTCount, node.assignments=GS.data$S.i, opts=opts)
    # Change the outfiles prefix to be able to identify the output files
    outfiles.prefix = file.path(outdir, paste(samplename, "_reassign_option_2", sep=""))
    
  } else if (mut.assignment.type == 3) {  
    warning("binom mut assignment not implemented for multiple timepoints")
    q(save="no", status=1)
    
  } else {
    warning(paste("Unknown mutation assignment type", mut.assignment.type, sep=" "))
    q(save="no", status=1)
  }
  
  # Make the figure
  replot_nD(outdir, outfiles.prefix, samplename, subsamples, dataset, clustering, conc_param, cluster_conc)
  
  return(list(clustering=clustering, outfiles.prefix=outfiles.prefix))
}

#' Helper function to read in the DP algorithm files
#' @param indir The directory in which the files were stored
#' @return A GS.data object with fields S.i, V.h, pi.h and alpha
#' @author sd11
read_gsdata_object = function(indir, no.iters, conc_param, cluster_conc) {
  GS.data = list()
  filename_prefix = file.path(outdir, paste(samplename, "_2D_iters", no.iters, "_concParam", conc_param, "_clusterWidth", 1/cluster_conc, sep=""))
  GS.data$S.i = as.matrix(read.csv(paste(filename_prefix, "_states.csv", sep=""), row.names=1))
  GS.data$V.h = as.matrix(read.csv(paste(filename_prefix, "_stickbreaking_weights.csv", sep=""), row.names=1))
  GS.data$pi.h = as.matrix(read.csv(paste(filename_prefix, "_discreteMutationCopyNumbers.csv", sep=""), row.names=1))
  GS.data$alpha = read.csv(paste(filename_prefix, "_alpha.csv", sep=""), row.names=1)
  
  # Fix the dimensions of this array
  pi.h = array(NA, c(nrow(GS.data$pi.h), ncol(GS.data$pi.h), 1))
  pi.h[1:nrow(GS.data$pi.h), 1:ncol(GS.data$pi.h), 1] = GS.data$pi.h
  GS.data$pi.h = pi.h
  return(GS.data)
}