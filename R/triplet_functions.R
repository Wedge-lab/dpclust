
#########################################################################################################
# Pipeline functions
#########################################################################################################

#' Helper function to package triplet pipeline run parameters
#' @param no.iters The number of iterations that the MCMC chain should be run for (Default: 1000). Used by steps "core" and "density", not used by step "combine". Must be identical between steps "core" and "density"
#' @param no.iters.burn.in The number of iterations that should be discarded as burn in of the MCMC chain (Default: ceiling(no.iters/5)). Used by steps "core" and "density", not used by step "combine". Must be identical between steps "core" and "density"
#' @param density.smooth Smoothing parameter used during multi-dimensional density estimation (Default: 0.01) - internal tuning parameter, not intended to be changed by regular users. Must be identical between steps "density" and "combine"
#' @param keep_md_out Boolean, set to FALSE to remove the "md_out" working directory (the per-triplet intermediate files written by step "density") once step "combine" has copied its final result files out to outdir (Default: TRUE). Only used by step "combine"
#' @return A list containing these components
#' @author sd11
#' @export
make_triplet_params = function(no.iters=1000, no.iters.burn.in=NULL, density.smooth=0.01, keep_md_out=TRUE) {
  if (is.null(no.iters.burn.in)) { no.iters.burn.in = ceiling(no.iters/5) }
  return(list(no.iters=no.iters, no.iters.burn.in=no.iters.burn.in, density.smooth=density.smooth, keep_md_out=keep_md_out))
}

#' Main entry point for the triplet pipeline. This pipeline is used when a donor has more
#' samples than the native nD pipeline can practically cluster jointly. Work is split into
#' three steps, mirroring RunDP's analysis_type dispatch: "core" (step 1, run the DP core
#' algorithm on all subsamples), "density" (step 2, estimate density for one triplet of
#' subsamples - to be run once per possible triplet) and "combine" (step 3, combine the
#' density estimates of all triplets into a single consensus).
#' @param step One of "core", "density", "combine"
#' @param sample_params List with sample parameters (see make_sample_params). Required for all steps
#' @param triplet_params List with triplet run parameters (see make_triplet_params). Required for all steps
#' @param outdir Directory where the output will be saved
#' @param advanced_params List with advanced parameters (see make_advanced_params). Only used by step "core"
#' @param subsample.indices The indices (within sample_params$subsamples) of the subsamples that make up this triplet. Required for step "density", not used otherwise
#' @author sd11
#' @export
RunTripletDP = function(step, sample_params, triplet_params, outdir, advanced_params=NULL, subsample.indices=NULL) {

  supported_steps = c("core", "density", "combine")
  if (!(step %in% supported_steps)) {
    stop(paste("Unknown triplet pipeline step '", step, "', specify one of: ", paste(supported_steps, collapse=", "), sep=""))
  }

  if (!file.exists(outdir)) { dir.create(outdir, recursive=T) }

  if (step == "core") {
    return(run_dpc_core_triplet(sample_params=sample_params, triplet_params=triplet_params, advanced_params=advanced_params, outdir=outdir))
  } else if (step == "density") {
    if (is.null(subsample.indices)) { stop("subsample.indices is required for step 'density'") }
    return(run_parallel_density_est(sample_params=sample_params, triplet_params=triplet_params, outdir=outdir, subsample.indices=subsample.indices))
  } else if (step == "combine") {
    return(combine_triplet_runs(sample_params=sample_params, triplet_params=triplet_params, outdir=outdir))
  }
}

#' Step 1 of the triplet pipeline: run the DP core algorithm on all subsamples
#' @param sample_params List with sample parameters (see make_sample_params)
#' @param triplet_params List with triplet run parameters (see make_triplet_params)
#' @param advanced_params List with advanced parameters (see make_advanced_params)
#' @param outdir Directory where the output is to be stored
#' @author Naser Ansari-Pour (BDI, Oxford), sd11
#' @noRd
run_dpc_core_triplet = function(sample_params, triplet_params, advanced_params, outdir) {

  set.seed(advanced_params$seed)

  dat = load_triplet_data(sample_params$datafiles, sample_params$cellularity)

  DirichletProcessClustering(mutCount = dat$mutCount,
                             WTCount = dat$WTCount,
                             totalCopyNumber = dat$totalCopyNumber,
                             copyNumberAdjustment = dat$copyNumberAdjustment,
                             mutation.copy.number = dat$mutation.copy.number,
                             cellularity = sample_params$cellularity,
                             output_folder = outdir,
                             no.iters = triplet_params$no.iters,
                             no.iters.burn.in = triplet_params$no.iters.burn.in,
                             subsamplesrun = sample_params$subsamples,
                             samplename = sample_params$samplename,
                             conc_param = advanced_params$conc_param,
                             cluster_conc = advanced_params$cluster_conc,
                             mut.assignment.type = 1,
                             most.similar.mut = NA,
                             mutationTypes = NA,
                             max.considered.clusters = advanced_params$max.considered.clusters)
}

#' Step 2 of the triplet pipeline: estimate density for one triplet of subsamples
#' @param sample_params List with sample parameters (see make_sample_params)
#' @param triplet_params List with triplet run parameters (see make_triplet_params)
#' @param outdir Directory that contains the output of step "core", and where this step's output will be stored (in a "md_out" subdirectory)
#' @param subsample.indices The indices (within sample_params$subsamples) of the subsamples that make up this triplet
#' @author Naser Ansari-Pour (BDI, Oxford), sd11
#' @noRd
run_parallel_density_est = function(sample_params, triplet_params, outdir, subsample.indices) {

  gsdata_file = file.path(outdir, paste0(sample_params$samplename, "_gsdata.RData"))

  dat = load_triplet_data(sample_params$datafiles, sample_params$cellularity)
  mutation.copy.number = dat$mutation.copy.number[,subsample.indices,drop=F]
  copyNumberAdjustment = dat$copyNumberAdjustment[,subsample.indices,drop=F]

  load(gsdata_file)
  pi.h = GS.data$pi.h[,,subsample.indices,drop=F]
  GS.data = list(S.i=GS.data$S.i, V.h=GS.data$V.h, pi.h=pi.h)

  new_output_folder = file.path(outdir, "md_out")
  if (!file.exists(new_output_folder)) { dir.create(new_output_folder) }

  opts = list(samplename=sample_params$samplename,
              subsamplenames=sample_params$subsamples[subsample.indices],
              no.iters=triplet_params$no.iters,
              no.iters.burn.in=triplet_params$no.iters.burn.in,
              outdir=new_output_folder)
  res = multiDimensionalClustering(mutation.copy.number=mutation.copy.number,
                             copyNumberAdjustment=copyNumberAdjustment,
                             GS.data=GS.data,
                             density.smooth=triplet_params$density.smooth,
                             opts=opts,
                             outfilesuffix=paste0("_",paste(subsample.indices,collapse="_")))
}

#' Step 3 of the triplet pipeline: combine the density estimates of all triplets into a single consensus
#' @param sample_params List with sample parameters (see make_sample_params)
#' @param triplet_params List with triplet run parameters (see make_triplet_params)
#' @param outdir Directory that contains the output of step "density" (in its "md_out" subdirectory), and where this step's output will be stored
#' @author Naser Ansari-Pour (BDI, Oxford), sd11
#' @noRd
combine_triplet_runs = function(sample_params, triplet_params, outdir) {

  density.smooth = triplet_params$density.smooth
  new_output_folder = file.path(outdir, "md_out")

  loaded = load_triplet_data(sample_params$datafiles, sample_params$cellularity)
  dat = list(data=loaded$data, mutCount=loaded$mutCount, WTCount=loaded$WTCount,
             totalCopyNumber=loaded$totalCopyNumber, copyNumberAdjustment=loaded$copyNumberAdjustment,
             mutation.copy.number=loaded$mutation.copy.number, cellularity=sample_params$cellularity,
             samplename=sample_params$samplename, subsamples=sample_params$subsamples)
  attach(dat)

  no.subsamples = length(subsamples)
  choose.number = 3
  choose.from = no.subsamples
  no.perms = choose(no.subsamples,3)

  no.muts = nrow(data[[1]]) #naser: nrow(info) until runDP filter has changed from -exclude.indices to non.zero ??
  print(paste(no.muts))
  node.assignments=NULL
  for(choose.index in 1:no.perms){
    print(choose.index)
    subsample.indices = get.subsample.indices(choose.index,choose.number,choose.from)
    print(subsample.indices)
    perm.assignments = read.table(paste(new_output_folder,"/",samplename,"_",paste(subsample.indices,collapse="_"),"_DP_and_cluster_info_0.01.txt",sep=""),header=T,sep="\t")

    if(is.null(node.assignments)){
      node.assignments = array(data.matrix(perm.assignments[ncol(perm.assignments)]),c(no.muts,1))
    }else{
      node.assignments = cbind(node.assignments,data.matrix(perm.assignments[ncol(perm.assignments)]))
    }
  }

  #get different combinations of assignments from the different permutations
  unique.assignments = unique(node.assignments)
  write.table(unique.assignments,paste(new_output_folder,"/",samplename,"_matchedClustersInParallelRuns.txt",sep=""),sep="\t",row.names=F,quote=F)

  no.consensus.nodes = nrow(unique.assignments)
  print(paste("#consensus nodes = ",no.consensus.nodes,sep=""))
  consensus.assignments = vector(mode="numeric",length = no.consensus.nodes)
  for(n in 1:no.consensus.nodes){
    print(n)
    consensus.assignments[sapply(1:no.muts,function(a,u,i){all(u==a[i,])},a=node.assignments,u=unique.assignments[n,])]=n
  }

  #get probabilities of assignment to each set of assignments
  prob.consensus.assignments = array(1,c(no.muts,no.consensus.nodes))
  #prob.consensus.assignments = array(0,c(no.muts,no.consensus.nodes)) #get average probability, rather than product
  for(choose.index in 1:no.perms){
    print(choose.index)
    subsample.indices = get.subsample.indices(choose.index,choose.number,choose.from)
    perm.assignments = read.table(paste(new_output_folder,"/",samplename,"_",paste(subsample.indices,collapse="_"),"_DP_and_cluster_info_0.01.txt",sep=""),header=T,sep="\t")
    for(n in 1:no.consensus.nodes){
      #prob.consensus.assignments[,n] = prob.consensus.assignments[,n] * perm.assignments[,unique.assignments[n,choose.index]+2]
      #prob.consensus.assignments[,n] = prob.consensus.assignments[,n] + perm.assignments[,unique.assignments[n,choose.index]+2]

      prob.consensus.assignments[,n] = prob.consensus.assignments[,n] * perm.assignments[,unique.assignments[n,choose.index]]
    }
  }
  #prob.consensus.assignments = prob.consensus.assignments / no.perms


  out = data[[1]][,1:2]
  for(s in 1:no.subsamples){
    out = cbind(out,data[[s]]$subclonal.fraction)
  }
  out = cbind(out,prob.consensus.assignments)
  names(out) = c("chr","pos",paste(samplename,subsamples,"_subclonalFraction",sep=""),paste("prob.cluster",1:no.consensus.nodes,sep=""))
  write.table(out,paste(new_output_folder,"/",samplename,"_allClusterProbabilitiesFromParallelRuns_16Oct2014.txt",sep=""),sep="\t",row.names=F,quote=F)

  ##############################################################################################################################################
  all.CIs = array(NA,c(no.consensus.nodes,no.subsamples,no.perms,2))
  for(choose.index in 1:no.perms){
    subsample.indices = get.subsample.indices(choose.index,choose.number,choose.from)
    print(subsample.indices)
    CIs = data.matrix(read.table(paste(new_output_folder,"/",samplename,"_",paste(subsample.indices,collapse="_"),"_confInts_",density.smooth,".txt",sep=""),header=T,sep="\t",stringsAsFactors=F))
    for(n in 1:no.consensus.nodes){
      all.CIs[n,subsample.indices,choose.index,1] = CIs[unique.assignments[n,choose.index],2*(1:choose.number)-1]
      all.CIs[n,subsample.indices,choose.index,2] = CIs[unique.assignments[n,choose.index],2*(1:choose.number)]
    }
  }
  median.CIs = array(NA,c(no.consensus.nodes,no.subsamples,2))
  for(n in 1:no.consensus.nodes){
    for(s in 1:no.subsamples){
      for(c in 1:2){
        median.CIs[n,s,c] = median(all.CIs[n,s,,c],na.rm=T)
      }
    }
  }
  median.CIs.2D = t(sapply(1:no.consensus.nodes,function(m,i){as.vector(t(m[i,,]))},m=median.CIs))
  write.table(cbind(1:no.consensus.nodes,median.CIs.2D,table(consensus.assignments)),sep="\t",row.names=F,quote=F,col.names=c("cluster.no",paste(rep(paste(samplename,subsamples,sep=""),each=2),rep(c("lowerCI","upperCI"),times=no.subsamples),sep="_"),"no.muts"),paste(new_output_folder,"/",samplename,"_consensusClustersByParallelNodeAssignment_16Oct2014.txt",sep=""))
  out = data[[1]][,1:2]
  for(s in 1:no.subsamples){
    out = cbind(out,data[[s]]$subclonal.fraction)
  }
  out = cbind(out,consensus.assignments)
  names(out) = c("chr","pos",paste(samplename,subsamples,"_subclonalFraction",sep=""),"cluster.no")
  write.table(out,paste(new_output_folder,"/",samplename,"_allClusterassignmentsFromParallelRuns_23Nov2018.txt",sep=""),sep="\t",row.names=F,quote=F)

  pdf(paste(new_output_folder,"/","consensus_cluster_assignment_",samplename,"_combined_",density.smooth,"_23Nov2018.pdf",sep=""),height=4,width=4)
  #its hard to distinguish more than 8 different colours
  max.cols = 8
  cols = rainbow(min(max.cols,no.consensus.nodes))
  plot.data = mutation.copy.number/copyNumberAdjustment
  plot.data[is.na(plot.data)]=0

  for(i in 1:(no.subsamples-1)){
    for(j in (i+1):no.subsamples){
      plot(plot.data[,i],plot.data[,j],type = "n",xlab = paste(samplename,(subsamples[subsample.indices])[i]," subclonal fraction",sep=""), ylab = paste(samplename,subsamples[j]," subclonal fraction",sep=""),xlim = c(0,max(plot.data[,i])*1.25))
      for(n in 1:no.consensus.nodes){
        pch=20 + floor((n-1)/max.cols)
        #pch is not implmeneted above 25
        if(pch>25){
          pch=pch-20
        }
        points(plot.data[,i][consensus.assignments==n],plot.data[,j][consensus.assignments==n],col=cols[(n-1) %% max.cols + 1],pch=pch,cex=0.5)
      }
      pch=20 + floor((0:(no.consensus.nodes-1))/max.cols)
      pch[pch>25] = pch[pch>25]-20
      legend(max(plot.data[,i])*1.05,max(plot.data[,j]),legend = 1:no.consensus.nodes,col=cols[(0:(no.consensus.nodes-1)) %% max.cols + 1],pch=pch,cex=1)
    }
  }
  dev.off()

  plot_consensus_assignment_heatmap(plot.data=plot.data,
                                    consensus.assignments=consensus.assignments,
                                    subsamples=subsamples,
                                    pngFile=paste(new_output_folder,"/",samplename,"_consensusAssignmentHeatmap.png",sep=""))

  # Copy the final consensus result files out of the "md_out" working directory
  # (which also contains many per-triplet intermediate files) into outdir
  final_files = c(paste(samplename,"_allClusterProbabilitiesFromParallelRuns_16Oct2014.txt",sep=""),
                  paste(samplename,"_consensusClustersByParallelNodeAssignment_16Oct2014.txt",sep=""),
                  paste(samplename,"_allClusterassignmentsFromParallelRuns_23Nov2018.txt",sep=""),
                  paste("consensus_cluster_assignment_",samplename,"_combined_",density.smooth,"_23Nov2018.pdf",sep=""),
                  paste(samplename,"_consensusAssignmentHeatmap.png",sep=""))
  file.copy(file.path(new_output_folder, final_files), outdir, overwrite=TRUE)

  if (!triplet_params$keep_md_out) {
    unlink(new_output_folder, recursive=TRUE)
  }
}

#' Heatmap of per-mutation CCF/subclonal fraction across samples, with mutations
#' sorted and grouped by consensus cluster assignment (clusters ordered by
#' breadth of presence across samples, most-shared cluster at the top) and
#' samples ordered by hierarchical clustering on their per-cluster CCF profile,
#' so samples that share the most clusters sit adjacent. Uses lattice::levelplot,
#' consistent with plotnD() in R/PlotDensities.R, rather than introducing
#' pheatmap/RColorBrewer as new package dependencies.
#' @param plot.data Matrix (mutations x samples) of subclonal fraction values, NAs already replaced with 0
#' @param consensus.assignments Vector, length = nrow(plot.data), giving each mutation's consensus cluster number
#' @param subsamples Vector of sample identifiers, one per column of plot.data (already includes samplename)
#' @param pngFile Output file path
#' @author sd11
#' @noRd
plot_consensus_assignment_heatmap = function(plot.data, consensus.assignments, subsamples, pngFile) {

  colnames(plot.data) = subsamples

  # Order clusters by breadth of presence across samples (most-shared first),
  # ties broken by cluster size (more mutations first)
  presence_threshold = 0.1
  cluster_ids = sort(unique(consensus.assignments))
  breadth = sapply(cluster_ids, function(cl) {
    per_sample_mean = colMeans(plot.data[consensus.assignments==cl,,drop=F])
    sum(per_sample_mean > presence_threshold)
  })
  cluster_size = sapply(cluster_ids, function(cl) sum(consensus.assignments==cl))
  ordered_cluster_ids = cluster_ids[order(breadth, cluster_size)]
  rank_map = setNames(seq_along(ordered_cluster_ids), ordered_cluster_ids)
  sort_key = rank_map[as.character(consensus.assignments)]

  cluster_order = order(sort_key)
  sorted.data = plot.data[cluster_order,,drop=F]
  cluster_boundaries = which(diff(sort_key[cluster_order]) != 0) + 0.5

  # Order samples by hierarchical clustering on their per-cluster CCF profile,
  # so samples that share the most clusters end up adjacent
  cluster_profile = t(sapply(cluster_ids, function(cl) colMeans(plot.data[consensus.assignments==cl,,drop=F])))
  sample_hc = hclust(dist(t(cluster_profile)))
  sorted.data = sorted.data[, sample_hc$order, drop=F]

  colours = colorRampPalette(c("white","red"))
  column_boundaries = seq(1.5, ncol(sorted.data)-0.5, by=1)
  panel_function = function(...) {
    lattice::panel.levelplot(...)
    lattice::panel.abline(h=cluster_boundaries, col="black")
    lattice::panel.abline(v=column_boundaries, col="black")
  }

  png(filename=pngFile, width=800, height=1200)
  fig = lattice::levelplot(t(sorted.data),
                           aspect="fill",
                           xlab=list(label="Sample", cex=2.5),
                           ylab=list(label="Mutation (sorted by cluster)", cex=2.5),
                           ylab.right=list(label="CCF", cex=2.5),
                           scales=list(x=list(cex=1.8, rot=90), y=list(draw=FALSE)),
                           col.regions=colours,
                           colorkey=list(space="right", height=0.4, labels=list(cex=1.8)),
                           panel=panel_function)
  print(fig)
  dev.off()
}

#########################################################################################################
# Helper functions
#########################################################################################################

#' Determines the total number of triplets from the given number of samples
#'
#' @param num_samples Integer total number of samples
#' @return Integer total of triplets
#' @author Naser Ansari-Pour (BDI, Oxford)
#' @export
determine_num_triplets = function(num_samples) {
	return(choose(num_samples,3))
}

#' Determines the indices of the samples to be used for a certain triplet
#'
#' @param choose.index Of all possible combinations of samples, run this combination
#' @param choose.number The number within the total triplets to be picked
#' @param choose.from The total number of triplets
#' @return Vector of integers
#' @author Naser Ansari-Pour (BDI, Oxford)
#' @export
get.subsample.indices<-function(choose.index,choose.number,choose.from){
      subsample.indices = 1:choose.number
      temp.choose.index = choose.index
      for(i in 1:choose.number){
              last.subtotal = 0
              subtotal = choose(choose.from-subsample.indices[i],choose.number-i)
              while(temp.choose.index>subtotal){
                      subsample.indices[i] = subsample.indices[i] + 1
                      last.subtotal = subtotal
                      subtotal = subtotal + choose(choose.from-subsample.indices[i],choose.number-i)
              }
              if(i<choose.number){
                      subsample.indices[i+1] = subsample.indices[i]+1
              }
              temp.choose.index = temp.choose.index - last.subtotal
      }
      return(subsample.indices)
}

#########################################################################################################
# Data parser
#########################################################################################################

#' Shared data loader for the triplet pipeline. Reads each subsample's dp-input file and
#' combines them into per-datatype matrices, keeping only mutations that pass the triplet
#' pipeline's QC filter. This is the same reading/combining/filtering logic previously
#' duplicated across the three triplet steps - used identically by all of them.
#' @param datafiles Vector of paths to per-subsample dp-input files, one per subsample
#' @param cellularity Vector with the purity/cellularity for each subsample
#' @return A list with combined per-datatype matrices (mutCount, WTCount, totalCopyNumber,
#' copyNumberAdjustment, mutation.copy.number), cellularity, and the underlying per-subsample tables (data)
#' @author Naser Ansari-Pour (BDI, Oxford), sd11
#' @noRd
load_triplet_data = function(datafiles, cellularity) {
  no.subsamples = length(datafiles)

  data = list()
  for (s in 1:no.subsamples) {
    data[[s]] = read.table(datafiles[s], sep="\t", header=T)
  }

  combine_column = function(colname) {
    do.call(cbind, lapply(data, function(x) x[[colname]]))
  }

  WTCount = combine_column("WT.count")
  mutCount = combine_column("mut.count")
  totalCopyNumber = combine_column("subclonal.CN")
  copyNumberAdjustment = combine_column("no.chrs.bearing.mut")
  mutation.copy.number = combine_column("mutation.copy.number")

  #naser: exclude unwanted variants using the "non.zero" filter (as opposed to the original "-exclude.indices") as applied in the naser_multidimclustparallel_updated_final.R
  non.zero = which(rowSums(mutCount)>0 & !is.na(rowSums(totalCopyNumber)) & rowSums(copyNumberAdjustment==0)==0)
  mutCount = mutCount[non.zero,,drop=F]
  WTCount = WTCount[non.zero,,drop=F]
  totalCopyNumber = totalCopyNumber[non.zero,,drop=F]
  copyNumberAdjustment = copyNumberAdjustment[non.zero,,drop=F]
  mutation.copy.number = mutation.copy.number[non.zero,,drop=F]
  for(s in 1:no.subsamples){
    data[[s]] = data[[s]][non.zero,,drop=F]
  }

  return(list(data=data, mutCount=mutCount, WTCount=WTCount, totalCopyNumber=totalCopyNumber,
              copyNumberAdjustment=copyNumberAdjustment, mutation.copy.number=mutation.copy.number,
              cellularity=cellularity))
}
