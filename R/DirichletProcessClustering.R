#
# Functions to run DPClust in various modes
#

#' Helper function to package run parameters
#' @param no.iters The number of iterations that the MCMC chain should be run for
#' @param no.iters.burn.in The number of iterations that should be discarded as burn in of the MCMC chain
#' @param mut.assignment.type Mutation assignment type option
#' @param num_muts_sample The number of mutations from which to start downsampling
#' @param min_muts_cluster The minimum number of mutations required for a cluster to be kept in the final output (Default: NULL)
#' @param min_frac_muts_cluster The minimum fraction of mutations required for a cluster to be kept in the final output (Default: 0.01)
#' @param species Species (Default: Human)
#' @param is.male Boolean set to TRUE when the donor is male, female otherwise
#' @param assign_sampled_muts A boolean whether to assign mutations that have not been used for clustering due to downsampling (Default: TRUE)
#' @param supported_chroms Vector with chromosome names from which mutations can be used (Default: NULL)
#' @param keep_temp_files Set to TRUE to keep temporary files (Default: TRUE)
#' @param generate_cluster_ordering Set to TRUE to generate possible phylogenetic relationships between clusters (Default: FALSE)
#' @return A list containing these components
#' @author sd11
#' @export
make_run_params = function(no.iters, no.iters.burn.in, mut.assignment.type, num_muts_sample, is.male, min_muts_cluster=NULL, min_frac_muts_cluster=0.01, species="human", assign_sampled_muts=TRUE, supported_chroms=NULL, keep_temp_files=TRUE, generate_cluster_ordering=FALSE) {
  if (is.null(supported_chroms)) {
    if (species=="human" | species=="Human") {
      # Set the expected chromosomes based on the sex
      if (is.male) {
        supported_chroms = as.character(c(1:22, "X", "Y"))
      } else {
        supported_chroms = as.character(c(1:22, "X"))
      }
    } else {
      if (species=="mouse" | species=="Mouse") {
        # Set the expected chromosomes based on the sex
        if (is.male) {
          supported_chroms = as.character(c(1:19, "X", "Y"))
        } else {
          supported_chroms = as.character(c(1:19, "X"))
        }
      }
    }
  }
  
  return(list(no.iters=no.iters, no.iters.burn.in=no.iters.burn.in, mut.assignment.type=mut.assignment.type, 
              supported_chroms=supported_chroms, num_muts_sample=num_muts_sample, assign_sampled_muts=assign_sampled_muts, keep_temp_files=keep_temp_files,
              generate_cluster_ordering=generate_cluster_ordering, species=species, min_muts_cluster=min_muts_cluster, min_frac_muts_cluster=min_frac_muts_cluster))
}

#' Helper function to package sample parameters
#' @param datafiles Vector of data files to be read in
#' @param cellularity Vector with purity for all samples
#' @param is.male Boolean whether the sample is male
#' @param samplename Donorname, used in plots and to name output files
#' @param subsamples Samplenames, used for multi-dimensional clustering to denote individual samples (this information is thought of as a suffix to the samplename, so: samplename=PD4120 and subsamplename=c(a, c)
#' @return A list containing these components
#' @author sd11
#' @export
make_sample_params = function(datafiles, cellularity, is.male, samplename, subsamples, mutphasingfiles=NULL) {
  return(list(datafiles=datafiles, cellularity=cellularity, is.male=is.male, samplename=samplename, subsamples=subsamples, mutphasingfiles=mutphasingfiles))
}

#' Helper function to package advanced parameters - most of these will almost never need to be changed
#' @param seed The seed to use
#' @param conc_param Hyperparameter setting that affects the sampling of the alpha stick-breaking parameter
#' @param cluster_conc Legacy parameter, no longer used
#' @param max.considered.clusters The maximum number of clusters to be considered
#' @return A list containing these components
#' @author sd11
#' @export
make_advanced_params = function(seed, conc_param=0.01, cluster_conc=5, max.considered.clusters=20) {
  return(list(conc_param=conc_param, cluster_conc=cluster_conc, seed=seed, max.considered.clusters=max.considered.clusters))
}

#' Helper function to package CNA parameters - to be implemented
make_cna_params = function() {
  print("Not yet implemented")
}

#' Main DPClust function that handles the various pipelines
#' @param analysis_type Type of analysis to run: nd_dp (1d and nd clustering), replot_1d/replot_nd (recreate plots), reassign_muts_1d/reassign_muts_nd (reassign mutations)
#' @param run_params List with run parameters (see make_run_params)
#' @param sample_params List with sample parameters (see make_sample_params)
#' @param advanced_params List with advanced parameters (see make_advanced_params)
#' @param outdir Directory where the output will be saved
#' @param cna_params List with copy number parameters - currently unsupported (Default: NULL)
#' @param mutphasingfiles Mutation phasing files - currently unsupported (Default: NULL)
#' @author sd11
#' @export
RunDP <- function(analysis_type, run_params, sample_params, advanced_params, outdir, cna_params=NULL, mutphasingfiles=NULL) { 
  
  #####################################################################################
  # Unpack parameters
  #####################################################################################
  attach(run_params)
  attach(sample_params)
  attach(advanced_params)
  if (!is.null(cna_params)) {
    attach(cna_params)
  }
  
  #####################################################################################
  # Check input
  #####################################################################################
  # Check whether a supported analysis_type was supplied
  supported_commands = c('nd_dp', 'replot_1d', 'replot_nd', 'reassign_muts_1d', 'reassign_muts_nd')
  if (!(analysis_type %in% supported_commands)) {
    print(paste("Type of analysis", analysis_type, "unknown."))
    print(paste(c("Specify either ", supported_commands)), sep=" ")
    q(save="no", status=1)
  }
  # Check whether the mut.assignment.type is supported
  supported_mut.assignment.methods = c(1,2,3,4)
  if (!(mut.assignment.type %in% supported_mut.assignment.methods)) {
    print(paste("Type of mutation assignment method", mut.assignment.type, "unknown."))
    print(paste(c("Specify either ", supported_mut.assignment.methods)), sep=" ")
    q(save="no", status=1)
  }
  
  #####################################################################################
  # Setup
  #####################################################################################
  set.seed(seed)
  
  # Create the output directory
  if (!file.exists(outdir)) { dir.create(outdir) }
  
  # Path for output files
  outfiles.prefix = file.path(outdir, paste(samplename, "_", no.iters, "iters_", no.iters.burn.in, "burnin", sep=""))
  
  #####################################################################################
  # Loading data
  #####################################################################################
  print("Loading data...")
  # Load data from disk if there was already a dataset object, otherwise create a new one
  if (file.exists(paste(outdir, "/dataset.RData", sep=""))) {
    
    load(file.path(outdir, "dataset.RData"))
    cndata = dataset$cndata
    mutphasing = dataset$mutphasing
    
    if (!is.null(cndata)) {
      cndata_params = list()
      cndata_params$cndata = cndata
      cndata_params$add.conflicts = add.conflicts
      cndata_params$cna.conflicting.events.only = cna.conflicting.events.only
      cndata_params$num.clonal.events.to.add = num.clonal.events.to.add
      cndata_params$min.cna.size = min.cna.size
    } else {
      cndata_params = NULL
    }
    
    
  } else {
    list_of_datafiles = paste(datpath, datafiles, sep="/")
    cndatafiles = paste(datpath, cndatafiles, sep="")  
    # Note that the phase column is not used
    dataset = load.data(list_of_datafiles, 
                        cellularity=cellularity, 
                        Chromosome="chr", 
                        position="end",
                        WT.count="WT.count", 
                        mut.count="mut.count", 
                        subclonal.CN="subclonal.CN", 
                        no.chrs.bearing.mut="no.chrs.bearing.mut", 
                        mutation.copy.number="mutation.copy.number", 
                        subclonal.fraction="subclonal.fraction", 
                        phase=NULL, # disabled for now, as the data is not used
                        is.male=is.male,
                        is.vcf=F, # reading of VCF input files os disabled
                        ref.genome.version="hg19", # reading of VCF input files is disabled, this parameter is not used
                        supported_chroms=supported_chroms)
    
    if (co_cluster_cna & !is.na(cndatafiles)) {
      cndata = load.cn.data(cndatafiles)
      cndata_params = list()
      cndata_params$cndata = cndata
      cndata_params$add.conflicts = add.conflicts
      cndata_params$cna.conflicting.events.only = cna.conflicting.events.only
      cndata_params$num.clonal.events.to.add = num.clonal.events.to.add
      cndata_params$min.cna.size = min.cna.size
    } else {
      cndata = NULL
      cndata_params = NULL
    }
    
    if (!is.null(mutphasingfiles)) {
      print("Loading mutation phasing info")
      mutphasing = NULL
      for (infile in mutphasingfiles) {
        mutphasing = rbind(mutphasing, read.table(infile, header=T, stringsAsFactors=F))
      }
    } else {
      mutphasing = NULL
    }
    
    save(file=file.path(outdir, "dataset.RData"), dataset)
  }
  
  
  #####################################################################################
  # Loading CN data and deal with the aftermath of that
  #####################################################################################
  # Unpack the copy number inclusion parameters
  if (!is.null(cndata_params)) {
    cndata = cndata_params$cndata
    add.conflicts = cndata_params$add.conflicts
    cna.conflicting.events.only = cndata_params$cna.conflicting.events.only
    num.clonal.events.to.add = cndata_params$num.clonal.events.to.add
    min.cna.size = cndata_params$min.cna.size
  }
  
  # Check if co-clustering of copy number data is in order
  resave.dataset = F # A boolean that keeps track of whether the dataset should be saved again. Set this to TRUE if the dataset changes.
  if (!is.null(dataset$cndata)) {
    # In case of a rerun, pull out the cndata
    cndata = dataset$cndata
  } else if (!is.null(cndata)) {
    dataset = add.in.cn.as.snv.cluster(dataset, 
                                       cndata, 
                                       cellularity=cellularity,
                                       add.conflicts=add.conflicts, 
                                       conflicting.events.only=cna.conflicting.events.only, 
                                       num.clonal.events.to.add=num.clonal.events.to.add,
                                       min.cna.size=min.cna.size)
    resave.dataset = T
  }
  
  #####################################################################################
  # Add in mutation phasing
  #####################################################################################
  # Check for mutationphasing info
  if (!is.null(dataset$mutphasing)) {
    mutphasing = dataset$mutphasing
  } else if (!is.null(mutphasing)) {
    dataset = add.mutphasing(dataset, mutphasing, add.conflicts=add.conflicts)
    resave.dataset = T
  }
  
  # Perform sampling
  if (!is.na(num_muts_sample) & num_muts_sample!="NA") {
    if (is.na(dataset$full.data)) {
      dataset = sample_mutations(dataset, num_muts_sample, sample.snvs.only=sample.snvs.only, remove.snvs=remove.snvs)
      most.similar.mut = dataset$most.similar.mut
      resave.dataset = T
    }       
    most.similar.mut = dataset$most.similar.mut
  } else {
    most.similar.mut = NA
  }
  dataset$cndata = cndata
  # The dataset object was modified, so save it
  if (resave.dataset) { save(file=file.path(outdir, "dataset.RData"), dataset) }
  
  if (analysis_type == 'nd_dp') {
    print("Running DPClust...")
    ##############################
    # nD DP clustering
    ##############################
    clustering = DirichletProcessClustering(mutCount=dataset$mutCount, 
                                            WTCount=dataset$WTCount, 
                                            no.iters=no.iters, 
                                            no.iters.burn.in=no.iters.burn.in, 
                                            cellularity=cellularity, 
                                            totalCopyNumber=dataset$totalCopyNumber, 
                                            mutation.copy.number=dataset$mutation.copy.number,
                                            copyNumberAdjustment=dataset$copyNumberAdjustment, 
                                            mutationTypes=dataset$mutationType,
                                            samplename=samplename, 
                                            subsamplesrun=subsamples,
                                            output_folder=outdir, 
                                            conc_param=conc_param, 
                                            cluster_conc=cluster_conc,
                                            mut.assignment.type=mut.assignment.type,
                                            most.similar.mut=most.similar.mut,
                                            max.considered.clusters=max.considered.clusters)
    
  } else if (analysis_type == "replot_1d") {
    print("Running Remaking plots...")
    ##############################
    # Replot 1D clustering
    ##############################
    density = read.table(file.path(outdir, paste(samplename, "_DirichletProcessplotdensity.txt", sep="")), header=T)
    polygon.data = read.table(file.path(outdir, paste(samplename, "_DirichletProcessplotpolygonData.txt", sep="")), header=T)
    load(file=paste(outfiles.prefix, "_bestConsensusResults.RData", sep=""))
    replot_1D(outdir=outdir,
              outfiles.prefix=outfiles.prefix,
              samplename=samplename, 
              dataset=dataset, 
              clustering=clustering, 
              density=density, 
              polygon.data=polygon.data)
    
  } else if (analysis_type == "replot_nd") {  
    print("Remaking plots...")
    ##############################
    # Replot nD clustering
    ##############################
    load(file=file.path(outdir, paste(samplename, "_gsdata.RData", sep="")))
    replot_nD(outdir=outdir, 
              outfiles.prefix=outfiles.prefix, 
              samplename=samplename, 
              subsamples=subsamples, 
              dataset=dataset,
              no.iters=no.iters,
              clustering=clustering, 
              conc_param=conc_param, 
              cluster_conc=cluster_conc)
    
  } else if (analysis_type == "reassign_muts_1d") {
    print("Reassigning mutations...")
    ##############################
    # Reassign mutations to clusters using a previous 1D clustering run
    ##############################
    load(file=file.path(outdir, paste(samplename, "_gsdata.RData", sep="")))
    res = reassign_1D(outdir=outdir, 
                      samplename=samplename, 
                      no.iters=no.iters, 
                      no.iters.burn.in=no.iters.burn.in, 
                      dataset=dataset, 
                      cellularity=cellularity, 
                      GS.data=GS.data,
                      conc_param=conc_param, 
                      cluster_conc=cluster_conc,
                      mut.assignment.type=mut.assignment.type)
    clustering = res$clustering
    outfiles.prefix = res$outfiles.prefix
    
  } else if (analysis_type == "reassign_muts_nd") {
    print("Reassigning mutations...")
    ##############################
    # Reassign mutations to clusters using a previous nD clustering run
    ##############################
    load(file=paste(outfiles.prefix, "_bestConsensusResults.RData", sep=""))
    res = reassign_nd(outdir=outdir, 
                      samplename=samplename, 
                      subsamples=subsamples, 
                      no.iters=no.iters, 
                      no.iters.burn.in=no.iters.burn.in, 
                      dataset=dataset, 
                      GS.data=GS.data,
                      conc_param=conc_param, 
                      cluster_conc=cluster_conc,
                      mut.assignment.type=mut.assignment.type)
    
    clustering = res$clustering
    outfiles.prefix = res$outfiles.prefix
    
  } else {
    print(paste("Unknown type of analysis",analysis_type))
    q(save="no", status=1)
  }
  
  ####################################################################################################################
  # Write the final output
  ####################################################################################################################
  if (all(!analysis_type %in% c('replot_1d', 'replot_nd'))) {
    print("Writing out final output...")
    
    # Load the MCMC output as its needed to get cluster confidence intervals
    density_file = file.path(outdir, paste(samplename, "_DirichletProcessplotdensity.txt", sep=""))
    if (file.exists(density_file)) {
      density = read.table(density_file, header=T)
      colnames(density)[1] = "fraction.of.tumour.cells"
    } else {
      density = NA
    }
    
    polygon_file = file.path(outdir, paste(samplename, "_DirichletProcessplotpolygonData.txt", sep=""))
    if (file.exists(polygon_file)) {
      polygon.data = read.table(polygon_file, header=T)
    } else {
      polygon.data = NA
    }
    
    load(file.path(outdir, paste(samplename, "_gsdata.RData", sep="")))
    write_tree = analysis_type != 'nd_dp' & analysis_type != 'reassign_muts_1d' & analysis_type != 'reassign_muts_nd'
    writeStandardFinalOutput(clustering=clustering, 
                             dataset=dataset,
                             most.similar.mut=most.similar.mut,
                             outfiles.prefix=outfiles.prefix,
                             outdir=outdir,
                             samplename=samplename,
                             subsamplenames=subsamples,
                             GS.data=GS.data,
                             density=density,
                             polygon.data=polygon.data,
                             no.iters=no.iters, 
                             no.iters.burn.in=no.iters.burn.in,
                             assign_sampled_muts=assign_sampled_muts,
                             write_tree=write_tree,
                             generate_cluster_ordering=generate_cluster_ordering,
                             min_muts_cluster=min_muts_cluster,
                             min_frac_muts_cluster=min_frac_muts_cluster)
  }
  
  ####################################################################################################################
  # Remove intermediate files
  ####################################################################################################################
  if (!keep_temp_files) {
    .remove_file = function(filename) {
      if (file.exists(filename)) file.remove(filename)
    }
    
    # 1D method and general files
    .remove_file(paste(outfiles.prefix, "_removedMutationsIndex.txt", sep=""))
    .remove_file(file.path(outdir, paste(samplename, "_DP_and_cluster_info.txt", sep="")))
    .remove_file(file.path(outdir, paste(samplename, "_DirichletProcessplot.png", sep="")))
    .remove_file(file.path(outdir, paste(samplename, "_DirichletProcessplot_with_cluster_locations.png", sep="")))
    .remove_file(file.path(outdir, paste(samplename, "_DirichletProcessplotdensity.txt", sep="")))
    .remove_file(file.path(outdir, paste(samplename, "_DirichletProcessplotpolygonData.txt", sep="")))
    .remove_file(file.path(outdir, paste(samplename, "_localOptima.txt", sep="")))
    .remove_file(file.path(outdir, paste(samplename, "_optimaInfo.txt", sep="")))
    .remove_file(file.path(outdir, paste(samplename, "_gsdata.RData", sep="")))
    .remove_file(file.path(outdir, "dataset.RData"))
    
    # nD method files
    .remove_file(file.path(outdir, paste(samplename, "_DP_and cluster_info_0.01.txt", sep="")))
    .remove_file(file.path(outdir, paste(samplename, "_confInts_0.01.txt", sep="")))
    .remove_file(file.path(outdir, paste(samplename, "_localHighConfidenceMultidimensionalOptima_0.01.txt", sep="")))
    .remove_file(file.path(outdir, paste(samplename, "_localMultidimensionalOptima_0.01.txt", sep="")))
    .remove_file(file.path(outdir, paste(samplename, "_optimaInfo_0.01.txt", sep="")))
    
    for (i in 1:(length(subsamples)-1)) {
      for (j in (i+1):length(subsamples)) {
        .remove_file(file.path(outdir, paste(samplename, subsamples[i], subsamples[j], "_densityoutput.RData", sep="")))
        .remove_file(file.path(outdir, paste(samplename, subsamples[i], subsamples[j], "_densityoutput.csv", sep="")))
        .remove_file(file.path(outdir, pattern=glob2rx(paste(samplename, subsamples[i], subsamples[j], "*densityData1.csv", sep="")), full.names=T))
        density_csv_files = list.files(outdir, pattern=glob2rx(paste(samplename, subsamples[i], subsamples[j], "*vals.csv", sep="")), full.names=T)
        for (infile in density_csv_files) { .remove_file(infile) }
      }
    }
    
    nd_density_files = list.files(outdir, pattern="_2D_binomial_")
    if (length(nd_density_files) > 0) { file.remove(nd_density_files) }
  }
  print("Done.")
}

#' Function that stores the final output in a unified format on disk
#' @param clustering A clustering result
#' @param dataset The dataset that went into clustering
#' @param most.similar.mut Vector containing for each non-sampled mutation its most similar sampled mutation. The non-sampled mutation will be assigned to the same cluster
#' @param outfiles.prefix A prefix for the filenames
#' @param outdir Output directory where the replot of the 1D method is to be stored if clusters are removed due to being too small
#' @param samplename Overall samplename
#' @param subsamplenames Samplenames of the different timepoints
#' @param GS.data MCMC output
#' @param density Posterior density estimate across MCMC iterations
#' @param polygon.data 1D confidence interval, used for plotting the 1D method (set to NA for multi-D method)
#' @param no.iters Number of total iterations
#' @param no.iters.burn.in Number of iterations to be used as burn-in
#' @param min_muts_cluster The minimum number of mutations required for a cluster to be kept in the final output
#' @param min_frac_muts_cluster The minimum fraction of mutations required for a cluster to be kept in the final
#' @param assign_sampled_muts Boolean whether to assign the non-sampled mutations (Default: TRUE)
#' @param write_tree Boolean whether to write a tree to file. Not all clustering methods return a tree (Default: FALSE)
#' @param generate_cluster_ordering Boolean specifying whether a possible cluster ordering should be determined (Default: FALSE)
#' @param no.samples.cluster.order Number of mutations to sample (with replacement) to classify pairs of clusters into parent-offspring or siblings (Default: 1000)
#' @author sd11
writeStandardFinalOutput = function(clustering, dataset, most.similar.mut, outfiles.prefix, outdir, samplename, subsamplenames, GS.data, density, polygon.data, no.iters, no.iters.burn.in, min_muts_cluster, min_frac_muts_cluster, assign_sampled_muts=T, write_tree=F, generate_cluster_ordering=F, no.samples.cluster.order=1000) {
  num_samples = ncol(dataset$mutCount)
  
  if(num_samples > 1 & generate_cluster_ordering == TRUE){
    stop("If run dpclust for multisample, setting the parameter 'generate_cluster_orders' as TRUE will result in an error. Please set 'generate_cluster_orders = FALSE'.")
  }
  
  ########################################################################
  # Check for too small clusters
  ########################################################################
  if (nrow(clustering$cluster.locations) > 1 & (min_muts_cluster!=-1 | min_frac_muts_cluster!=-1)) {
    if (min_muts_cluster!=-1 & min_frac_muts_cluster!=-1) print("Found entries for both min_muts_cluster and min_frac_muts_cluster, used which yielded the largest number")
    
    # min_muts_cluster = ifelse(is.null(min_muts_cluster), -1, min_muts_cluster)
    min_frac_muts_cluster = ifelse(min_frac_muts_cluster==-1, -1, min_frac_muts_cluster*nrow(dataset$mutCount))
    
    # use min_muts_cluster variable further down, overwrite with min_frac_muts_cluster
    if (min_frac_muts_cluster > min_muts_cluster) {
      min_muts_cluster = min_frac_muts_cluster
    }
    
    clusters_to_remove = clustering$cluster.locations[,num_samples+2] < min_muts_cluster
    
    if (sum(clusters_to_remove) > 0 & sum(clusters_to_remove) < length(clusters_to_remove)) {
      # remove clusters that are too small
      clusterids_to_remove = clustering$cluster.locations[clusters_to_remove,1]
      new_cluster.locations = clustering$cluster.locations[!clustering$cluster.locations[,1] %in% clusterids_to_remove,,drop=F]
      new_all.assignment.likelihoods = clustering$all.assignment.likelihoods[,!clusters_to_remove, drop=F]
      new_best.assignment.likelihoods = clustering$best.assignment.likelihoods
      new_best.node.assignments = clustering$best.node.assignments
      # reset best likelihoods and hard assignments for mutations assigned to the removed cluster(s)
      new_best.assignment.likelihoods[new_best.node.assignments %in% clusterids_to_remove] = NA
      new_best.node.assignments[new_best.node.assignments %in% clusterids_to_remove] = NA
      
      clustering$cluster.locations = new_cluster.locations
      clustering$all.assignment.likelihoods = new_all.assignment.likelihoods
      clustering$best.node.assignments = new_best.node.assignments
      clustering$best.assignment.likelihoods = new_best.assignment.likelihoods
      
      # if 1D clustering, then replot without the removed cluster
      if (ncol(dataset$mutCount)==1) {
        # Old plot
        plot1D(density=density, 
               polygon.data=polygon.data[,1], 
               pngFile=paste(outdir, "/", samplename, "_DirichletProcessplot_with_cluster_locations.png", sep=""), 
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
                 pngFile=paste(outdir, "/", samplename, "_DirichletProcessplot_with_cluster_locations_2.png", sep=""), 
                 density.from=0, 
                 x.max=1.5, 
                 mutationCopyNumber=dataset$mutation.copy.number, 
                 no.chrs.bearing.mut=dataset$copyNumberAdjustment,
                 samplename=samplename,
                 cluster.locations=clustering$cluster.locations,
                 mutation.assignments=clustering$best.node.assignments,
                 mutationTypes=dataset$mutationType)
      }
    }
  }
    
  if (generate_cluster_ordering) {
    ########################################################################
    # Before doing anything else, calculate confidence intervals on the cluster locations using only the mutations used during clustering
    ########################################################################
    # Calc confidence intervals for the cluster locations
    conf = calc_cluster_conf_intervals(GS.data, 
                                       mut_assignments=clustering$best.node.assignments, 
                                       clusterids=clustering$cluster.locations[,1], 
                                       no.muts=nrow(dataset$mutCount), 
                                       no.timepoints=ncol(dataset$mutCount), 
                                       no.iters=no.iters, 
                                       no.iters.burn.in=no.iters.burn.in)
    conf = data.frame(conf)
    colnames(conf) = c("cluster.no", "timepoint", "loc_conf_0.025", "loc_conf_0.500", "loc_conf_0.975")
    write.table(conf, file=paste(outfiles.prefix, "_clusterConfidenceIntervals.txt", sep=""), quote=F, row.names=F, sep="\t")
    
    # Calc probs for cluster orders
    probs = calc_cluster_order_probs(GS.data=GS.data,
                                     density=density,
                                     mut_assignments=clustering$best.node.assignments,
                                     clusterids=clustering$cluster.locations[,1],
                                     cluster_ccfs=clustering$cluster.locations[,2],
                                     no.muts=nrow(dataset$mutCount),
                                     no.timepoints=ncol(dataset$mutCount),
                                     no.iters=no.iters,
                                     no.iters.burn.in=no.iters.burn.in,
                                     no.samples=no.samples.cluster.order)
    probs = flatten_3d_to_2d(probs$classification, c("timepoint", clustering$cluster.locations[,1]))
    write.table(probs, file=paste(outfiles.prefix, "_clusterOrderProbabilities.txt", sep=""), quote=F, row.names=F, sep="\t")
  }
  
  ########################################################################
  # Check if mutation sampling has been done, if so, unpack and assign here
  ########################################################################
  if (!is.na(most.similar.mut) && assign_sampled_muts) {
    res = unsample_mutations(dataset, clustering)
    dataset = res$dataset
    clustering = res$clustering
  }
  
  ########################################################################
  # Write out the final mutation-cluster probabilities with all mutations spiked in
  ########################################################################
  # if (!is.null(clustering$all.assignment.likelihoods) & !is.na(clustering$all.assignment.likelihoods)) {
  if ("all.assignment.likelihoods" %in% names(clustering)) {
    # Fetch and drop all columns that have just zeroes
    cols_all_zero = which(sapply(1:ncol(clustering$all.assignment.likelihoods), function(i) { max(clustering$all.assignment.likelihoods[,i])==0 }))
    if (length(cols_all_zero)!=0) {
      all_assignment_likelihoods = clustering$all.assignment.likelihoods[,-cols_all_zero, drop=F]
      cluster_colnames = (1:ncol(clustering$all.assignment.likelihoods))
      cluster_colnames = (1:ncol(clustering$all.assignment.likelihoods))[-cols_all_zero]
    } else {
      all_assignment_likelihoods = clustering$all.assignment.likelihoods
      cluster_colnames = 1:ncol(clustering$all.assignment.likelihoods)
    }
    
    all_assignment_likelihoods = data.frame(dataset$chromosome[,1], dataset$position[,1]-1, dataset$position[,1], all_assignment_likelihoods, clustering$best.node.assignments)
    colnames(all_assignment_likelihoods) = c("chr", "start", "end", paste("prob.cluster", cluster_colnames, sep="."), "most.likely.cluster")
    write.table(all_assignment_likelihoods[dataset$mutationType=="SNV",], file=paste(outfiles.prefix, "_mutationClusterLikelihoods.bed", sep=""), quote=F, row.names=F, sep="\t")
    
    if (any(dataset$mutationType=="CNA")) {
      write.table(all_assignment_likelihoods[dataset$mutationType=="CNA",], file=paste(outfiles.prefix, "_mutationClusterLikelihoodsPseudoSNV.bed", sep=""), quote=F, row.names=F, sep="\t")
    }
  }
  
  ########################################################################
  # Write final cluster locations
  ########################################################################
  # Add the confidence intervals to the final cluster locations
  if (ncol(clustering$cluster.locations) > 3) {
    # nD based clustering
    write.table(clustering$cluster.locations, paste(outfiles.prefix, "_bestClusterInfo.txt",sep=""), col.names=c("cluster.no", paste(samplename, subsamplenames, sep=""), "no.of.mutations"), sep="\t", quote=F, row.names=F)
  } else {
    # 1D based
    write.table(clustering$cluster.locations, paste(outfiles.prefix, "_bestClusterInfo.txt",sep=""), col.names=c("cluster.no", "location", "no.of.mutations"), row.names=F, sep="\t", quote=F)
  }
  
  ########################################################################
  # Add the removed mutations back in
  ########################################################################
  output = cbind(dataset$chromosome[,1], dataset$position[,1]-1, dataset$position[,1], clustering$best.node.assignments, clustering$best.assignment.likelihoods)
  output = add_removed_snvs(dataset, output)
  
  ########################################################################
  # Save the indices of the mutations that were not used during the analysis
  ########################################################################
  save(file=paste(outfiles.prefix, "_bestConsensusResults.RData", sep=""), output, clustering, GS.data, density, dataset, most.similar.mut, no.iters, no.iters.burn.in)
  write.table(data.frame(mut.index=dataset$removed_indices), file=paste(outfiles.prefix,"_removedMutationsIndex.txt", sep=""), row.names=F, quote=F)
  
  ########################################################################
  # Save the consensus mutation assignments
  ########################################################################
  colnames(output) = c("chr", "start", "end", "cluster", "likelihood")
  write.table(output[dataset$mutationType=="SNV",], file=paste(outfiles.prefix, "_bestConsensusAssignments.bed", sep=""), quote=F, row.names=F, sep="\t")
  
  ########################################################################
  # Save the CNA assignments separately
  ########################################################################
  if (any(dataset$mutationType=="CNA")) {
    write.table(output[dataset$mutationType=="CNA",], file=paste(outfiles.prefix, "_bestConsensusAssignmentsPseudoSNV.bed", sep=""), quote=F, row.names=F, sep="\t")
    # Assign the CNAs to clusters using their pseudoSNV representations
    cndata = assign_cnas_to_clusters(dataset$cndata, output)
    write.table(cndata, file=paste(outfiles.prefix, "_bestCNAassignments.txt", sep=""), quote=F, row.names=F, sep="\t")
    
    # if (!is.null(clustering$all.assignment.likelihoods) & !is.na(clustering$all.assignment.likelihoods)) {
    if ("all.assignment.likelihoods" %in% names(clustering)) {
      cna_assignment_likelihoods = get_cnas_cluster_probs(dataset$cndata, all_assignment_likelihoods[dataset$mutationType=="CNA",], colnames(all_assignment_likelihoods))
      write.table(cna_assignment_likelihoods, file=paste(outfiles.prefix, "_cnaClusterLikelihoods.bed", sep=""), quote=F, row.names=F, sep="\t")
    }
    
    # Create a new assignment table figure with the correct information
    # This removes pseudo SNVs as the assignmentTable will add an extra column for CNAs
    cluster_locations = clustering$cluster.locations
    cluster_locations[,3] = rep(0, nrow(cluster_locations))
    mut_assignments = table(output[dataset$mutationType=="SNV","cluster"])
    for (i in 1:nrow(cluster_locations)) {
      if (as.character(cluster_locations[i,1]) %in% names(mut_assignments)) {
        cluster_locations[i,3] = mut_assignments[as.character(cluster_locations[i,1])]
      }
    }
    plotAssignmentTable(cluster_locations, paste(outfiles.prefix, "_mutation_assignments.png", sep=""), cndata=cndata, num_samples=num_samples)
  } else {
    plotAssignmentTable(clustering$cluster.locations, paste(outfiles.prefix, "_mutation_assignments.png", sep=""), num_samples=num_samples)
  }
  
  ########################################################################
  # If tree based analysis, also save the tree
  ########################################################################
  if (write_tree) {
    write.table(clustering$best.tree, file=paste(outfiles.prefix, "_bestConsensusTree.txt", sep=""), quote=F, row.names=F, sep="\t")
  }
}

#' Add removed mutations back into the assignment table. SNVs will be assigned to the cluster of its most similar not-removed SNV
#' @param dataset A dataset object
#' @param snv_assignment_table Data frame with the mutation assignments
#' @return The snv_assignment_table with the removed mutations added into the position they were originally
#' @author sd11
add_removed_snvs = function(dataset, snv_assignment_table) {
  if (length(dataset$removed_indices) > 0) {
    for (i in dataset$removed_indices) {
      if (i==1) {
        snv_assignment_table = rbind(c(dataset$chromosome.not.filtered[i], dataset$mut.position.not.filtered[i]-1, dataset$mut.position.not.filtered[i], NA, NA), snv_assignment_table)
      } else if (i >= nrow(snv_assignment_table)) {
        snv_assignment_table = rbind(snv_assignment_table, c(dataset$chromosome.not.filtered[i], dataset$mut.position.not.filtered[i]-1, dataset$mut.position.not.filtered[i], NA, NA))
      } else {
        snv_assignment_table = rbind(snv_assignment_table[1:(i-1),], c(dataset$chromosome.not.filtered[i], dataset$mut.position.not.filtered[i]-1, dataset$mut.position.not.filtered[i], NA, NA), snv_assignment_table[i:nrow(snv_assignment_table),])
      }
    }
  }
  
  # Sort the output in the same order as the dataset
  chrpos_input = paste(dataset$chromosome, dataset$position, sep="_")
  chrpos_output = paste(snv_assignment_table[,1], snv_assignment_table[,3], sep="_")
  snv_assignment_table = snv_assignment_table[match(chrpos_input, chrpos_output),]
  return(snv_assignment_table)
}


#' Assign CNA events to clusters using their pseudoSNV representation
#' @param cndata Data frame with the CNA data
#' @param snv_assignment_table Data frame with the mutation assignments, with chromosome and position-start/end expected as first three columns
#' @return The cndata object with extra column cluster_assignment
#' @author sd11
assign_cnas_to_clusters = function(cndata, snv_assignment_table) {
  cndata$cluster_assignment = "NA"
  for (i in 1:nrow(cndata)) {
    # Fetch all pseudoSNVs that represent this CNA
    selection = snv_assignment_table[,1]==as.character(cndata$chr[i]) & snv_assignment_table[,3]==as.character(cndata$startpos[i])
    if (any(selection)) {
      # Work out which cluster has the most pseudoSNVs assigned. That will be the cluster to which the CNA event is assigned
      pseudo_snvs = snv_assignment_table[selection,,drop=F]
      assign_inventory = table(pseudo_snvs[,4])
      cndata$cluster_assignment[i] = names(assign_inventory)[which.max(assign_inventory)]
    }
  }
  return(cndata)
}

#' Use the Pseudo-SNV probabilities to obtain a probability of each CNA of each cluster
#' @param cndata The copy number data data.frame
#' @param snv_assignment_likelihoods Probabilities of the pseudo-SNVs
#' @param cluster_colnames The colnames of the output bedfile to be used to select the assignment columns
#' @return A data.frame with chr,start,end,probs_per_cluster
#' @author sd11
get_cnas_cluster_probs = function(cndata, snv_assignment_likelihoods, cluster_colnames) {
  cna_cluster_probs = matrix(NA, nrow=nrow(cndata), ncol=ncol(snv_assignment_likelihoods)-4)
  for (i in 1:nrow(cndata)) {
    # Fetch all pseudoSNVs that represent this CNA
    selection = snv_assignment_likelihoods[,1]==as.character(cndata$chr[i]) & snv_assignment_likelihoods[,3]==as.character(cndata$startpos[i])
    if (any(selection)) {
      # Work out which cluster has the most pseudoSNVs assigned. That will be the cluster to which the CNA event is assigned
      pseudo_snvs = snv_assignment_likelihoods[selection,,drop=F]
      
      cna_cluster_probs[i,] = sapply(which(grepl("prob", cluster_colnames)), function(j) {
        if (any(pseudo_snvs[,j]==0)) {
          0
        } else {
          # Combine p-values using fishers' method
          pchisq(-2 * sum(log(pseudo_snvs[,j])), df=length(pseudo_snvs[,j]), lower.tail=F)
        }
      })
    }
  }
  # Obtain most likely cluster
  # First get those CNAs for which we don't have any probabilities and set them to NA
  assignments = apply(cna_cluster_probs, 1, function(x) { all(is.na(x)) })
  assignments[assignments] = NA
  # Then assign those for which we have probabilities to the cluster with the highest probability
  assignments[!is.na(assignments)] = unlist(apply(cna_cluster_probs, 1, which.max))
  
  output = data.frame(cndata[, c("chr", "startpos", "endpos", "CNA")], cna_cluster_probs, assignments)
  colnames(output) = c("chr", "startpos", "endpos", "CNA", cluster_colnames[grepl("prob", cluster_colnames)], "most.likely.cluster")
  
  return(output)
}

#' Helper function that flattens a 3D array into a 2D one
#' @param data The data to be flattened
#' @param col_names The names of the columns in the output
#' @return A data.frame with the third column annotated as the first column
#' @author sd11
flatten_3d_to_2d = function(data, col_names) {
  no.timepoints = dim(data)[3]
  no.clusters = dim(data)[2]
  new_data = data.frame(array(NA, c(no.clusters*no.timepoints, no.clusters+1)))
  for (i in 1:no.timepoints) {
    row = ((i-1)*no.clusters) + 1
    new_data[row:(row+no.clusters-1), 2:(no.clusters+1)] = data[,,i]
    new_data[row:(row+no.clusters-1), 1] = i
  }
  colnames(new_data) = col_names
  return(new_data)
}

#' Main function to run subclonal reconstruction
#' 
#' Will perform clustering using the given data. The method
#' decides automatically whether the 1D or nD method is run based on the number of samples given at the input.
#' The number of samples is determined through the number of columns of the input.
#' @param mutCount Matrix with readcounts of the mutated allele
#' @param WTCount Matrix with readcounts of the wild-type allele
#' @param totalCopyNumber Matrix with total copynumber at each mutation locus
#' @param copyNumberAdjustment Matrix with multiplicity values
#' @param mutation.copy.number Matrix with mutation copy number values
#' @param cellularity Vector with sample purities
#' @param output_folder Directory where to write output
#' @param no.iters The number of iterations to run the MCMC chain for
#' @param no.iters.burn.in Number of iterations to discard as burn in
#' @param samplename Donor name, used in plots and to name output files
#' @param subsamplesrun Samplenames of individual samples for this donor
#' @param conc_param Hyperparameter setting that affects the sampling of the alpha stick-breaking parameter
#' @param cluster_conc Legacy parameter, no longer used
#' @param mut.assignment.type Type of mutation assignment to be used
#' @param most.similar.mut Vector with most similar mutation for mutations removed during sampling (if any)
#' @param mutationTypes Vector with mutation types, used for plotting
#' @param max.considered.clusters Maximum number of clusters to consider
#' @author sd11
DirichletProcessClustering <- function(mutCount, WTCount, totalCopyNumber, copyNumberAdjustment, mutation.copy.number, cellularity, output_folder, no.iters, no.iters.burn.in, subsamplesrun, samplename, conc_param, cluster_conc, mut.assignment.type, most.similar.mut, mutationTypes, max.considered.clusters) {

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
                                    cluster_conc=cluster_conc,
                                    C=max.considered.clusters)
  
  save(file=file.path(output_folder, paste(samplename, "_gsdata.RData", sep="")), GS.data)
  
  # nD dataset, plot sample versus sample
  if (ncol(mutCount) > 1) {
    ########################
    # Plot density and Assign mutations to clusters - nD
    ########################
    print("Estimating density between pairs of samples...")
    for (i in 1:(length(subsamplesrun)-1)) {
      for (j in (i+1):length(subsamplesrun)) {
        print(paste("Samples", subsamplesrun[i], "and", subsamplesrun[j], sep=" "))
        imageFile = file.path(output_folder, paste(samplename,subsamplesrun[i],subsamplesrun[j],"_iters",no.iters,"_concParam",conc_param,"_clusterWidth",1/cluster_conc,"_2D_binomial.png",sep=""))
        density = Gibbs.subclone.density.est(mutation.copy.number[,c(i,j)]/copyNumberAdjustment[,c(i,j)],
                                             GS.data,
                                             imageFile,
                                             post.burn.in.start = no.iters.burn.in,
                                             post.burn.in.stop = no.iters,
                                             samplenames = paste(samplename,subsamplesrun[c(i,j)],sep=""),
                                             indices=c(i,j))
        save(file=file.path(output_folder, paste(samplename, subsamplesrun[i], subsamplesrun[j], "_densityoutput.RData", sep="")), GS.data, density)
      }
    }
    
    # Assign mutations to clusters using one of the different assignment methods
    print("Assigning mutations to clusters...")
    opts = list(samplename=samplename, subsamplenames=subsamplesrun, no.iters=no.iters, no.iters.burn.in=no.iters.burn.in, no.iters.post.burn.in=no.iters-no.iters.burn.in, outdir=output_folder)
    if (mut.assignment.type == 1) {
      consClustering = multiDimensionalClustering(mutation.copy.number=mutation.copy.number, 
                                                  copyNumberAdjustment=copyNumberAdjustment, 
                                                  GS.data=GS.data, 
                                                  density.smooth=0.01, 
                                                  opts=opts)
    } else if (mut.assignment.type == 2) {
      consClustering = mutation_assignment_em(GS.data=GS.data,
                                              mutCount=mutCount, 
                                              WTCount=WTCount, 
                                              subclonal.fraction=mutation.copy.number/copyNumberAdjustment, 
                                              node.assignments=GS.data$S.i, 
                                              opts=opts)
      
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

    print("Assigning mutations to clusters...")
    if (mut.assignment.type == 1) {
      subclonal.fraction = mutation.copy.number / copyNumberAdjustment
      subclonal.fraction[is.nan(subclonal.fraction)] = 0
      consClustering = oneDimensionalClustering(samplename, subclonal.fraction, GS.data, density, no.iters, no.iters.burn.in)
      setwd(wd) 
    } else if (mut.assignment.type == 2) {
      setwd(wd) # set the wd back earlier. The oneD clustering and gibbs sampler do not play nice yet and need the switch, the em assignment doesnt
      consClustering = mutation_assignment_em(GS.data=GS.data,
                                              mutCount=mutCount, 
                                              WTCount=WTCount, 
                                              subclonal.fraction=mutation.copy.number/copyNumberAdjustment, 
                                              node.assignments=GS.data$S.i, 
                                              opts=opts)
      
    } else if (mut.assignment.type == 3) {  
      consClustering = mutation_assignment_binom(clustering_density=density,
                                                 mutCount=mutCount, 
                                                 WTCount=WTCount, 
                                                 copyNumberAdjustment=copyNumberAdjustment, 
                                                 tumourCopyNumber=totalCopyNumber,
                                                 normalCopyNumber=array(2, dim(mutCount)),
                                                 cellularity=cellularity,
                                                 samplename=samplename)
      
    } else {
      warning(paste("Unknown mutation assignment type", mut.assignment.type, sep=" "))
      q(save="no", status=1)
    }

    # Make a second set of figures with the mutation assignments showing
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
    
    plot1D_2(density=density, 
             polygon.data=polygon.data, 
             pngFile=paste(output_folder, "/", samplename, "_DirichletProcessplot_with_cluster_locations_2.png", sep=""), 
             density.from=0, 
             x.max=1.5, 
             mutationCopyNumber=mutation.copy.number, 
             no.chrs.bearing.mut=copyNumberAdjustment,
             samplename=samplename,
             cluster.locations=consClustering$cluster.locations,
             mutation.assignments=consClustering$best.node.assignments,
             mutationTypes=mutationTypes)
    
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
    if (i==1) {
      dat = cbind(rep(def.value, nrow(dat)), dat)
      dat = rbind(rep(def.value, ncol(dat)), dat)
    } else if (i >= ncol(dat)) {
      	dat = cbind(dat, rep(def.value, nrow(dat)))
      	dat = rbind(dat, rep(def.value, ncol(dat)))
    } else {
        dat = cbind(dat[,1:(i-1)], rep(def.value, nrow(dat)), dat[,i:ncol(dat)])
        dat = rbind(dat[1:(i-1),], rep(def.value, ncol(dat)), dat[i:nrow(dat),])
    }
  }
  return(dat)
}

