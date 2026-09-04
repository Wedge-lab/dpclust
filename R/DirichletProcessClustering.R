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
#' @param memory_limit_gb Optional hard memory guard rail in GB. DPClust runs if feasible under this limit, otherwise errors before starting.
#' @param num_threads Number of CPU threads for the Gibbs C++ core (OpenMP). If NA, use runtime default.
#' @param sample.snvs.only Boolean whether to only sample from SNVs (Default: TRUE)
#' @param remove.snvs Boolean whether to remove all SNVs (Default: FALSE)
#' @return A list containing these components
#' @author sd11
#' @export
make_run_params <- function(no.iters, no.iters.burn.in, mut.assignment.type, num_muts_sample, is.male, min_muts_cluster = NULL, min_frac_muts_cluster = 0.01, species = "human", assign_sampled_muts = TRUE, supported_chroms = NULL, keep_temp_files = TRUE, generate_cluster_ordering = FALSE, memory_limit_gb = NA_real_, num_threads = NA_integer_, sample.snvs.only = TRUE, remove.snvs = FALSE, prefix = NULL, conc_param = 0.01, density_smooth = NA_real_, x_max_cap = 3, hypercube_size = 5, cluster_conc = 5) {
  if (is.null(supported_chroms)) {
    if (species == "human" | species == "Human") {
      # Set the expected chromosomes based on the sex
      if (is.male) {
        supported_chroms <- as.character(c(1:22, "X", "Y"))
      } else {
        supported_chroms <- as.character(c(1:22, "X"))
      }
    } else {
      if (species == "mouse" | species == "Mouse") {
        # Set the expected chromosomes based on the sex
        if (is.male) {
          supported_chroms <- as.character(c(1:19, "X", "Y"))
        } else {
          supported_chroms <- as.character(c(1:19, "X"))
        }
      }
    }
  }

  return(list(
    no.iters = no.iters, no.iters.burn.in = no.iters.burn.in, mut.assignment.type = mut.assignment.type,
    is.male = is.male,
    supported_chroms = supported_chroms, num_muts_sample = num_muts_sample, assign_sampled_muts = assign_sampled_muts, keep_temp_files = keep_temp_files,
    generate_cluster_ordering = generate_cluster_ordering, species = species, min_muts_cluster = min_muts_cluster, min_frac_muts_cluster = min_frac_muts_cluster,
    memory_limit_gb = memory_limit_gb, num_threads = num_threads, sample.snvs.only = sample.snvs.only, remove.snvs = remove.snvs,
    prefix = prefix,
    conc_param = conc_param,
    density_smooth = density_smooth,
    x_max_cap = x_max_cap,
    hypercube_size = hypercube_size,
    cluster_conc = cluster_conc
  ))
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
make_sample_params <- function(datafiles, cellularity, is.male, samplename, subsamples, mutphasingfiles = NULL, datpath = "", cndatafiles = NA) {
  return(list(datafiles = datafiles, cellularity = cellularity, is.male = is.male, samplename = samplename, subsamples = subsamples, mutphasingfiles = mutphasingfiles, datpath = datpath, cndatafiles = cndatafiles))
}

#' Helper function to package advanced parameters - most of these will almost never need to be changed
#' @param seed The seed to use
#' @param conc_param Hyperparameter setting that affects the sampling of the alpha stick-breaking parameter
#' @param cluster_conc Legacy parameter, no longer used
#' @param max.considered.clusters The maximum number of clusters to be considered
#' @return A list containing these components
#' @author sd11
#' @export
make_advanced_params <- function(seed, conc_param = 0.01, cluster_conc = 5, max.considered.clusters = 20) {
  return(list(conc_param = conc_param, cluster_conc = cluster_conc, seed = seed, max.considered.clusters = max.considered.clusters))
}

#' Helper function to package CNA parameters - to be implemented
make_cna_params <- function() {
  log_info("Not yet implemented")
}

#' Internal helper that logs data characteristics relevant for cluster separability.
#' This does not alter clustering behavior; it only provides observability in logs.
.log_clusterability_diagnostics <- function(dataset, samplename = "") {
  if (is.null(dataset$mutation.copy.number) || is.null(dataset$copyNumberAdjustment)) {
    return(invisible(NULL))
  }
  ccf <- as.numeric(dataset$mutation.copy.number / dataset$copyNumberAdjustment)
  ccf <- ccf[is.finite(ccf)]
  if (length(ccf) == 0) {
    log_info("Clusterability diagnostics: no finite CCF values available.")
    return(invisible(NULL))
  }

  q <- stats::quantile(ccf, probs = c(0.01, 0.05, 0.50, 0.95, 0.99), na.rm = TRUE, names = FALSE)
  mad_ccf <- stats::mad(ccf, constant = 1, na.rm = TRUE)
  prop_gt2 <- mean(ccf > 2, na.rm = TRUE)
  prop_gt3 <- mean(ccf > 3, na.rm = TRUE)
  n <- length(ccf)
  tag <- if (nzchar(samplename)) paste0(" (", samplename, ")") else ""

  log_info(sprintf(
    "Clusterability diagnostics%s: n=%d, CCF q01=%.3f q05=%.3f q50=%.3f q95=%.3f q99=%.3f, MAD=%.3f, >2=%.2f%%, >3=%.2f%%",
    tag, n, q[1], q[2], q[3], q[4], q[5], mad_ccf, 100 * prop_gt2, 100 * prop_gt3
  ))

  if (mad_ccf < 0.08) {
    warning("CCF distribution is very narrow (low MAD). Single-cluster solutions are more likely.", call. = FALSE)
  }
  if (prop_gt3 > 0.01) {
    warning("More than 1% of CCF values are >3; check multiplicity/copy-number inputs for calibration issues.", call. = FALSE)
  }
  invisible(NULL)
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
RunDP <- function(analysis_type, run_params, sample_params, advanced_params, outdir, cna_params = NULL, mutphasingfiles = NULL) {
  #####################################################################################
  # Unpack parameters
  #####################################################################################
  # Parameters are explicitly extracted from the input lists to avoid scoping/attachment issues

  # Explicitly extract parameters to avoid scoping/attachment issues in a package context
  no.iters <- run_params$no.iters
  no.iters.burn.in <- run_params$no.iters.burn.in
  mut.assignment.type <- run_params$mut.assignment.type
  num_muts_sample <- run_params$num_muts_sample
  is.male <- if ("is.male" %in% names(run_params) && !is.null(run_params$is.male)) run_params$is.male else TRUE
  min_muts_cluster <- if ("min_muts_cluster" %in% names(run_params) && !is.null(run_params$min_muts_cluster)) run_params$min_muts_cluster else -1
  min_frac_muts_cluster <- if ("min_frac_muts_cluster" %in% names(run_params) && !is.null(run_params$min_frac_muts_cluster)) run_params$min_frac_muts_cluster else 0.01
  density_smooth <- if ("density_smooth" %in% names(run_params)) run_params$density_smooth else NA_real_
  x_max_cap <- if ("x_max_cap" %in% names(run_params)) run_params$x_max_cap else 3
  hypercube_size <- if ("hypercube_size" %in% names(run_params)) run_params$hypercube_size else 5
  cluster_conc <- if ("cluster_conc" %in% names(run_params)) run_params$cluster_conc else advanced_params$cluster_conc
  conc_param <- if ("conc_param" %in% names(run_params)) run_params$conc_param else advanced_params$conc_param
  max.considered.clusters <- if ("max.considered.clusters" %in% names(advanced_params)) advanced_params$max.considered.clusters else 20
  species <- if ("species" %in% names(run_params) && !is.null(run_params$species)) run_params$species else "human"
  assign_sampled_muts <- if ("assign_sampled_muts" %in% names(run_params) && !is.null(run_params$assign_sampled_muts)) run_params$assign_sampled_muts else TRUE
  supported_chroms <- if ("supported_chroms" %in% names(run_params)) run_params$supported_chroms else NULL
  keep_temp_files <- if ("keep_temp_files" %in% names(run_params) && !is.null(run_params$keep_temp_files)) run_params$keep_temp_files else TRUE
  generate_cluster_ordering <- if ("generate_cluster_ordering" %in% names(run_params) && !is.null(run_params$generate_cluster_ordering)) run_params$generate_cluster_ordering else FALSE
  memory_limit_gb <- if ("memory_limit_gb" %in% names(run_params)) run_params$memory_limit_gb else NA_real_
  num_threads <- if ("num_threads" %in% names(run_params)) run_params$num_threads else NA_integer_
  sample.snvs.only <- if ("sample.snvs.only" %in% names(run_params)) run_params$sample.snvs.only else TRUE
  remove.snvs <- if ("remove.snvs" %in% names(run_params)) run_params$remove.snvs else FALSE
  prefix <- if ("prefix" %in% names(run_params)) run_params$prefix else NULL

  samplename <- sample_params$samplename
  datafiles <- sample_params$datafiles
  subsamples <- sample_params$subsamples
  cellularity <- sample_params$cellularity
  datpath <- if ("datpath" %in% names(sample_params)) sample_params$datpath else ""
  cndatafiles <- if ("cndatafiles" %in% names(sample_params)) sample_params$cndatafiles else NA

  seed <- advanced_params$seed
  co_cluster_cna <- if (!is.null(cna_params) && ("co_cluster_cna" %in% names(cna_params))) cna_params$co_cluster_cna else FALSE
  add.conflicts <- if (!is.null(cna_params) && ("add.conflicts" %in% names(cna_params))) cna_params$add.conflicts else FALSE
  cna.conflicting.events.only <- if (!is.null(cna_params) && ("cna.conflicting.events.only" %in% names(cna_params))) cna_params$cna.conflicting.events.only else FALSE
  num.clonal.events.to.add <- if (!is.null(cna_params) && ("num.clonal.events.to.add" %in% names(cna_params))) cna_params$num.clonal.events.to.add else 1
  min.cna.size <- if (!is.null(cna_params) && ("min.cna.size" %in% names(cna_params))) cna_params$min.cna.size else 100

  #####################################################################################
  # Check input
  #####################################################################################
  # Check whether a supported analysis_type was supplied
  supported_commands <- c("nd_dp", "replot_1d", "replot_nd", "reassign_muts_1d", "reassign_muts_nd")
  if (!(analysis_type %in% supported_commands)) {
    log_info(paste("Type of analysis", analysis_type, "unknown."))
    log_info(paste(c("Specify either ", supported_commands), collapse = " "))
    stop(paste("Type of analysis", analysis_type, "unknown."))
  }
  # Check whether the mut.assignment.type is supported
  supported_mut.assignment.methods <- c(1, 2, 3, 4)
  if (!(mut.assignment.type %in% supported_mut.assignment.methods)) {
    log_info(paste("Type of mutation assignment method", mut.assignment.type, "unknown."))
    log_info(paste(c("Specify either ", supported_mut.assignment.methods), collapse = " "))
    stop(paste("Type of mutation assignment method", mut.assignment.type, "unknown."))
  }

  #####################################################################################
  # Setup
  #####################################################################################
  set.seed(seed)

  # Create the output directory (use absolute, recursive to work in containers)
  outdir <- normalizePath(outdir, mustWork = FALSE)
  if (!dir.exists(outdir)) {
    log_info(paste("Creating output directory:", outdir))
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  }

  # Path for output files
  prefix_str <- if (!is.null(prefix) && !is.na(prefix) && nchar(prefix) > 0) paste0("_", prefix) else ""
  outfiles.prefix <- file.path(outdir, paste0(samplename, prefix_str, "_", no.iters, "iters_", no.iters.burn.in, "burnin"))

  #####################################################################################
  # Loading data
  #####################################################################################
  log_info("Loading data...")
  # Load data from disk if there was already a dataset object, otherwise create a new one
  # Use a stable filename based on prefix (if any) to ensure we can overwrite/reuse instead of creating timestamped clutter.
  prefix_tag <- if (!is.null(prefix) && !is.na(prefix) && nchar(prefix) > 0) paste0(prefix, "_") else ""
  rdata_file_name <- paste0(prefix_tag, "Seed-", seed, "_dataset.RData")
  if (file.exists(paste(outdir, "/", rdata_file_name, sep = ""))) {
    load(file.path(outdir, rdata_file_name))
    cndata <- dataset$cndata
    mutphasing <- dataset$mutphasing

    if (!is.null(cndata)) {
      cndata_params <- list()
      cndata_params$cndata <- cndata
      cndata_params$add.conflicts <- add.conflicts
      cndata_params$cna.conflicting.events.only <- cna.conflicting.events.only
      cndata_params$num.clonal.events.to.add <- num.clonal.events.to.add
      cndata_params$min.cna.size <- min.cna.size
    } else {
      cndata_params <- NULL
    }
  } else {
    # Build data file paths - if datpath is non-empty, prepend it, otherwise use datafiles as-is
    # This ensures both absolute paths (e.g. /data/file.txt) and relative paths work
    list_of_datafiles <- if (nchar(datpath) > 0) file.path(datpath, datafiles) else datafiles
    cndatafiles <- if (nchar(datpath) > 0 && !all(is.na(cndatafiles))) file.path(datpath, cndatafiles) else cndatafiles
    # Note that the phase column is not used
    dataset <- load.data(list_of_datafiles,
      cellularity = cellularity,
      Chromosome = "chr",
      position = "end",
      WT.count = "WT.count",
      mut.count = "mut.count",
      subclonal.CN = "subclonal.CN",
      no.chrs.bearing.mut = "no.chrs.bearing.mut",
      mutation.copy.number = "mutation.copy.number",
      subclonal.fraction = "subclonal.fraction",
      phase = NULL, # disabled for now, as the data is not used
      is.male = is.male,
      is.vcf = FALSE, # reading of VCF input files os disabled
      ref.genome.version = "hg19", # reading of VCF input files is disabled, this parameter is not used
      supported_chroms = supported_chroms
    )

    if (co_cluster_cna & !is.na(cndatafiles)) {
      cndata <- load.cn.data(cndatafiles)
      cndata_params <- list()
      cndata_params$cndata <- cndata
      cndata_params$add.conflicts <- add.conflicts
      cndata_params$cna.conflicting.events.only <- cna.conflicting.events.only
      cndata_params$num.clonal.events.to.add <- num.clonal.events.to.add
      cndata_params$min.cna.size <- min.cna.size
    } else {
      cndata <- NULL
      cndata_params <- NULL
    }

    if (!is.null(mutphasingfiles)) {
      log_info("Loading mutation phasing info")
      mutphasing <- NULL
      for (infile in mutphasingfiles) {
        mutphasing <- rbind(mutphasing, read.table(infile, header = TRUE, stringsAsFactors = FALSE))
      }
    } else {
      mutphasing <- NULL
    }

    save(file = file.path(outdir, rdata_file_name), dataset)
  }


  #####################################################################################
  # Loading CN data and deal with the aftermath of that
  #####################################################################################
  # Unpack the copy number inclusion parameters
  if (!is.null(cndata_params)) {
    cndata <- cndata_params$cndata
    add.conflicts <- cndata_params$add.conflicts
    cna.conflicting.events.only <- cndata_params$cna.conflicting.events.only
    num.clonal.events.to.add <- cndata_params$num.clonal.events.to.add
    min.cna.size <- cndata_params$min.cna.size
  }

  # Check if co-clustering of copy number data is in order
  resave.dataset <- FALSE # A boolean that keeps track of whether the dataset should be saved again. Set this to TRUE if the dataset changes.
  if (!is.null(dataset$cndata)) {
    # In case of a rerun, pull out the cndata
    cndata <- dataset$cndata
  } else if (!is.null(cndata)) {
    dataset <- add.in.cn.as.snv.cluster(dataset,
      cndata,
      cellularity = cellularity,
      add.conflicts = add.conflicts,
      conflicting.events.only = cna.conflicting.events.only,
      num.clonal.events.to.add = num.clonal.events.to.add,
      min.cna.size = min.cna.size
    )
    resave.dataset <- TRUE
  }

  #####################################################################################
  # Add in mutation phasing
  #####################################################################################
  # Check for mutationphasing info
  if (!is.null(dataset$mutphasing)) {
    mutphasing <- dataset$mutphasing
  } else if (!is.null(mutphasing)) {
    dataset <- add.mutphasing(dataset, mutphasing, add.conflicts = add.conflicts)
    resave.dataset <- TRUE
  }

  # Perform sampling
  if (!is.na(num_muts_sample) && !identical(num_muts_sample, "NA")) {
    if (.is_na_sentinel(dataset$full.data)) {
      dataset <- sample_mutations(dataset, num_muts_sample, sample.snvs.only = sample.snvs.only, remove.snvs = remove.snvs)
      most.similar.mut <- dataset$most.similar.mut
      resave.dataset <- TRUE
    }
    most.similar.mut <- dataset$most.similar.mut
  } else {
    most.similar.mut <- NA
  }
  dataset$cndata <- cndata
  .log_clusterability_diagnostics(dataset, samplename = samplename)
  # The dataset object was modified, so save it
  if (resave.dataset) {
    save(file = file.path(outdir, rdata_file_name), dataset)
  }

  if (analysis_type == "nd_dp") {
    log_info("Running DPClust...")
    thin_s_i <- mut.assignment.type %in% c(1, 4) && !generate_cluster_ordering
    keep_aux_fields <- FALSE
    memory_guard <- .memory_guard_plan(
      no.muts = nrow(dataset$mutCount),
      no.samples = ncol(dataset$mutCount),
      no.iters = no.iters,
      no.iters.burn.in = no.iters.burn.in,
      max.considered.clusters = max.considered.clusters,
      thin_s_i = thin_s_i,
      keep_aux_fields = keep_aux_fields,
      memory_limit_gb = memory_limit_gb,
      verbose = TRUE
    )
    thin_s_i <- memory_guard$thin_s_i
    keep_aux_fields <- memory_guard$keep_aux_fields
    ##############################
    # nD DP clustering
    ##############################
    clustering <- DirichletProcessClustering(
      mutCount = dataset$mutCount,
      WTCount = dataset$WTCount,
      no.iters = no.iters,
      no.iters.burn.in = no.iters.burn.in,
      cellularity = cellularity,
      totalCopyNumber = dataset$totalCopyNumber,
      mutation.copy.number = dataset$mutation.copy.number,
      copyNumberAdjustment = dataset$copyNumberAdjustment,
      mutationTypes = dataset$mutationType,
      samplename = samplename,
      subsamplesrun = subsamples,
      output_folder = outdir,
      conc_param = conc_param,
      cluster_conc = cluster_conc,
      mut.assignment.type = mut.assignment.type,
      most.similar.mut = most.similar.mut,
      max.considered.clusters = max.considered.clusters,
      thin_s_i = thin_s_i,
      keep_aux_fields = keep_aux_fields,
      num_threads = num_threads,
      conflict.array = dataset$conflict.array,
      keep_temp_files = keep_temp_files,
      density_smooth = density_smooth,
      x_max_cap = x_max_cap,
      hypercube_size = hypercube_size
    )
    GS.data <- clustering$GS.data
  } else if (analysis_type == "replot_1d") {
    log_info("Running Remaking plots...")
    ##############################
    # Replot 1D clustering
    ##############################
    density <- read.table(file.path(outdir, paste(samplename, "_DirichletProcessplotdensity.txt", sep = "")), header = TRUE)
    polygon.data <- read.table(file.path(outdir, paste(samplename, "_DirichletProcessplotpolygonData.txt", sep = "")), header = TRUE)
    load(file = paste(outfiles.prefix, "_bestConsensusResults.RData", sep = ""))
    replot_1D(
      outdir = outdir,
      outfiles.prefix = outfiles.prefix,
      samplename = samplename,
      dataset = dataset,
      clustering = clustering,
      density = density,
      polygon.data = polygon.data,
      x_max_cap = x_max_cap
    )
  } else if (analysis_type == "replot_nd") {
    log_info("Remaking plots...")
    ##############################
    # Replot nD clustering
    ##############################
    load(file = file.path(outdir, paste(samplename, "_gsdata.RData", sep = "")))
    replot_nD(
      outdir = outdir,
      outfiles.prefix = outfiles.prefix,
      samplename = samplename,
      subsamples = subsamples,
      dataset = dataset,
      no.iters = no.iters,
      clustering = clustering,
      conc_param = conc_param,
      cluster_conc = cluster_conc
    )
  } else if (analysis_type == "reassign_muts_1d") {
    log_info("Reassigning mutations...")
    ##############################
    # Reassign mutations to clusters using a previous 1D clustering run
    ##############################
    load(file = file.path(outdir, paste(samplename, "_gsdata.RData", sep = "")))
    res <- reassign_1D(
      outdir = outdir,
      samplename = samplename,
      no.iters = no.iters,
      no.iters.burn.in = no.iters.burn.in,
      dataset = dataset,
      cellularity = cellularity,
      GS.data = GS.data,
      conc_param = conc_param,
      cluster_conc = cluster_conc,
      mut.assignment.type = mut.assignment.type,
      x_max_cap = x_max_cap
    )
    clustering <- res$clustering
    outfiles.prefix <- res$outfiles.prefix
  } else if (analysis_type == "reassign_muts_nd") {
    log_info("Reassigning mutations...")
    ##############################
    # Reassign mutations to clusters using a previous nD clustering run
    ##############################
    load(file = paste(outfiles.prefix, "_bestConsensusResults.RData", sep = ""))
    if (!exists("GS.data")) {
      gsdata_file <- file.path(outdir, paste(samplename, "_gsdata.RData", sep = ""))
      if (file.exists(gsdata_file)) {
        load(file = gsdata_file)
      }
    }
    if (!exists("GS.data")) {
      stop("Cannot reassign mutations: GS.data not found in either bestConsensusResults or *_gsdata.RData.")
    }
    res <- reassign_nd(
      outdir = outdir,
      samplename = samplename,
      subsamples = subsamples,
      no.iters = no.iters,
      no.iters.burn.in = no.iters.burn.in,
      dataset = dataset,
      GS.data = GS.data,
      conc_param = conc_param,
      cluster_conc = cluster_conc,
      mut.assignment.type = mut.assignment.type
    )

    clustering <- res$clustering
    outfiles.prefix <- res$outfiles.prefix
  } else {
    log_info(paste("Unknown type of analysis", analysis_type))
    stop(paste("Unknown type of analysis", analysis_type))
  }

  ####################################################################################################################
  # Write the final output
  ####################################################################################################################
  if (all(!analysis_type %in% c("replot_1d", "replot_nd"))) {
    log_info("Writing out final output...")

    # Load the MCMC output as its needed to get cluster confidence intervals
    density_file <- file.path(outdir, paste(samplename, "_DirichletProcessplotdensity.txt", sep = ""))
    if (file.exists(density_file)) {
      density <- read.table(density_file, header = TRUE)
      colnames(density)[1] <- "fraction.of.tumour.cells"
    } else {
      density <- NA
    }

    polygon_file <- file.path(outdir, paste(samplename, "_DirichletProcessplotpolygonData.txt", sep = ""))
    if (file.exists(polygon_file)) {
      polygon.data <- read.table(polygon_file, header = TRUE)
    } else {
      polygon.data <- NA
    }

    if (!exists("GS.data") || is.null(GS.data)) {
      GS.data <- NULL
    }
    if (generate_cluster_ordering || !keep_temp_files) {
      if (is.null(GS.data)) {
        load(file.path(outdir, paste(samplename, "_gsdata.RData", sep = "")))
      }
    }
    write_tree <- analysis_type != "nd_dp" & analysis_type != "reassign_muts_1d" & analysis_type != "reassign_muts_nd"
    writeStandardFinalOutput(
      clustering = clustering,
      dataset = dataset,
      most.similar.mut = most.similar.mut,
      outfiles.prefix = outfiles.prefix,
      outdir = outdir,
      samplename = samplename,
      subsamplenames = subsamples,
      GS.data = GS.data,
      density = density,
      polygon.data = polygon.data,
      no.iters = no.iters,
      no.iters.burn.in = no.iters.burn.in,
      assign_sampled_muts = assign_sampled_muts,
      write_tree = write_tree,
      generate_cluster_ordering = generate_cluster_ordering,
      min_muts_cluster = min_muts_cluster,
      min_frac_muts_cluster = min_frac_muts_cluster,
      num_threads = num_threads,
      x_max_cap = x_max_cap
    )
  }

  ####################################################################################################################
  # Remove intermediate files
  ####################################################################################################################
  if (!keep_temp_files) {
    .remove_file <- function(filename) {
      if (file.exists(filename)) file.remove(filename)
    }

    # 1D method and general files
    .remove_file(paste(outfiles.prefix, "_removedMutationsIndex.txt", sep = ""))
    .remove_file(file.path(outdir, paste(samplename, "_DP_and_cluster_info.txt", sep = "")))
    .remove_file(file.path(outdir, paste(samplename, "_DirichletProcessplot.png", sep = "")))
    .remove_file(file.path(outdir, paste(samplename, "_DirichletProcessplot_with_cluster_locations.png", sep = "")))
    .remove_file(file.path(outdir, paste(samplename, "_DirichletProcessplotdensity.txt", sep = "")))
    .remove_file(file.path(outdir, paste(samplename, "_DirichletProcessplotpolygonData.txt", sep = "")))
    .remove_file(file.path(outdir, paste(samplename, "_localOptima.txt", sep = "")))
    .remove_file(file.path(outdir, paste(samplename, "_optimaInfo.txt", sep = "")))
    .remove_file(file.path(outdir, paste(samplename, "_gsdata.RData", sep = "")))
    .remove_file(file.path(outdir, rdata_file_name))

    # nD method files
    .remove_file(file.path(outdir, paste(samplename, "_DP_and cluster_info_0.01.txt", sep = "")))
    .remove_file(file.path(outdir, paste(samplename, "_confInts_0.01.txt", sep = "")))
    .remove_file(file.path(outdir, paste(samplename, "_localHighConfidenceMultidimensionalOptima_0.01.txt", sep = "")))
    .remove_file(file.path(outdir, paste(samplename, "_localMultidimensionalOptima_0.01.txt", sep = "")))
    .remove_file(file.path(outdir, paste(samplename, "_optimaInfo_0.01.txt", sep = "")))

    if (length(subsamples) > 1) {
      for (i in 1:(length(subsamples) - 1)) {
        for (j in (i + 1):length(subsamples)) {
          .remove_file(file.path(outdir, paste(samplename, subsamples[i], subsamples[j], "_densityoutput.RData", sep = "")))
          .remove_file(file.path(outdir, paste(samplename, subsamples[i], subsamples[j], "_densityoutput.csv", sep = "")))
          density_data1_files <- list.files(outdir, pattern = glob2rx(paste(samplename, subsamples[i], subsamples[j], "*densityData1.csv", sep = "")), full.names = TRUE)
          for (infile in density_data1_files) {
            .remove_file(infile)
          }
          density_csv_files <- list.files(outdir, pattern = glob2rx(paste(samplename, subsamples[i], subsamples[j], "*vals.csv", sep = "")), full.names = TRUE)
          for (infile in density_csv_files) {
            .remove_file(infile)
          }
        }
      }
    }

    nd_density_files <- list.files(outdir, pattern = "_2D_binomial_")
    if (length(nd_density_files) > 0) {
      file.remove(nd_density_files)
    }
  }
  log_info("Done.")
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
writeStandardFinalOutput <- function(clustering, dataset, most.similar.mut, outfiles.prefix, outdir, samplename, subsamplenames, GS.data, density, polygon.data, no.iters, no.iters.burn.in, min_muts_cluster, min_frac_muts_cluster, assign_sampled_muts = TRUE, write_tree = FALSE, generate_cluster_ordering = FALSE, no.samples.cluster.order = 1000, num_threads = NA_integer_, x_max_cap = 3) {
  num_samples <- ncol(dataset$mutCount)

  if (num_samples > 1 & generate_cluster_ordering == TRUE) {
    stop("If run dpclust for multisample, setting the parameter 'generate_cluster_orders' as TRUE will result in an error. Please set 'generate_cluster_orders = FALSE'.")
  }

  ########################################################################
  # Check for too small clusters
  ########################################################################
  if (nrow(clustering$cluster.locations) > 1 & (min_muts_cluster != -1 | min_frac_muts_cluster != -1)) {
    if (min_muts_cluster != -1 & min_frac_muts_cluster != -1) log_info("Found entries for both min_muts_cluster and min_frac_muts_cluster, used which yielded the largest number")

    # min_muts_cluster = ifelse(is.null(min_muts_cluster), -1, min_muts_cluster)
    min_frac_muts_cluster <- ifelse(min_frac_muts_cluster == -1, -1, min_frac_muts_cluster * nrow(dataset$mutCount))

    # use min_muts_cluster variable further down, overwrite with min_frac_muts_cluster
    if (min_frac_muts_cluster > min_muts_cluster) {
      min_muts_cluster <- min_frac_muts_cluster
    }

    clusters_to_remove <- clustering$cluster.locations[, num_samples + 2] < min_muts_cluster

    if (sum(clusters_to_remove) > 0 & sum(clusters_to_remove) < length(clusters_to_remove)) {
      # remove clusters that are too small
      clusterids_to_remove <- clustering$cluster.locations[clusters_to_remove, 1]
      new_cluster.locations <- clustering$cluster.locations[!clustering$cluster.locations[, 1] %in% clusterids_to_remove, , drop = FALSE]
      kept_clusterids <- new_cluster.locations[, 1]
      new_cluster.locations[, 1] <- seq_len(nrow(new_cluster.locations))
      if (.has_assignment_likelihoods(clustering)) {
        new_all.assignment.likelihoods <- clustering$all.assignment.likelihoods[, !clusters_to_remove, drop = FALSE]
      } else {
        new_all.assignment.likelihoods <- clustering$all.assignment.likelihoods
      }
      new_best.assignment.likelihoods <- clustering$best.assignment.likelihoods
      new_best.node.assignments <- clustering$best.node.assignments
      removed_assignment_mask <- !is.na(new_best.node.assignments) & (new_best.node.assignments %in% clusterids_to_remove)
      kept_mask <- !is.na(new_best.node.assignments) & !removed_assignment_mask
      new_best.node.assignments[kept_mask] <- match(new_best.node.assignments[kept_mask], kept_clusterids)
      # reset best likelihoods and hard assignments for mutations assigned to the removed cluster(s)
      new_best.assignment.likelihoods[removed_assignment_mask] <- NA
      new_best.node.assignments[removed_assignment_mask] <- NA

      clustering$cluster.locations <- new_cluster.locations
      clustering$all.assignment.likelihoods <- new_all.assignment.likelihoods
      clustering$best.node.assignments <- new_best.node.assignments
      clustering$best.assignment.likelihoods <- new_best.assignment.likelihoods

      # if 1D clustering, then replot without the removed cluster
      if (ncol(dataset$mutCount) == 1) {
        # Old plot
        plot1D(
          density = density,
          polygon.data = polygon.data[, 1],
          pngFile = paste(outdir, "/", samplename, "_DirichletProcessplot_with_cluster_locations.png", sep = ""),
          density.from = 0,
          x.max = NA,
          x.max.cap = x_max_cap,
          mutationCopyNumber = dataset$mutation.copy.number,
          no.chrs.bearing.mut = dataset$copyNumberAdjustment,
          samplename = samplename,
          cluster.locations = clustering$cluster.locations,
          mutation.assignments = clustering$best.node.assignments
        )

        # New plot
        plot1D_2(
          density = density,
          polygon.data = polygon.data[, 1],
          pngFile = paste(outdir, "/", samplename, "_DirichletProcessplot_with_cluster_locations_2.png", sep = ""),
          density.from = 0,
          x.max = NA,
          x.max.cap = x_max_cap,
          mutationCopyNumber = dataset$mutation.copy.number,
          no.chrs.bearing.mut = dataset$copyNumberAdjustment,
          samplename = samplename,
          cluster.locations = clustering$cluster.locations,
          mutation.assignments = clustering$best.node.assignments,
          mutationTypes = dataset$mutationType
        )
      }
    }
  }

  if (ncol(dataset$mutCount) == 1 && any(is.na(clustering$best.node.assignments))) {
    if (!is.null(GS.data) && !.is_na_sentinel(GS.data) && !is.null(density) && !.is_na_sentinel(density)) {
      log_info(sprintf(
        "Reassigning %d mutation(s) with NA cluster after 1D small-cluster filtering.",
        sum(is.na(clustering$best.node.assignments))
      ))
      clustering <- reassign_1d_na_mutations(
        clustering = clustering,
        GS.data = GS.data,
        density = density,
        no.iters = no.iters,
        no.iters.burn.in = no.iters.burn.in
      )
    } else {
      log_info("Skipping 1D NA reassignment because GS.data or density is unavailable.")
    }
  }

  if (generate_cluster_ordering) {
    if (is.null(GS.data) || .is_na_sentinel(GS.data)) {
      stop("GS.data is required when generate_cluster_ordering is TRUE.")
    }
    ########################################################################
    # Before doing anything else, calculate confidence intervals on the cluster locations using only the mutations used during clustering
    ########################################################################
    # Calc confidence intervals for the cluster locations
    conf <- calc_cluster_conf_intervals(GS.data,
      mut_assignments = clustering$best.node.assignments,
      clusterids = clustering$cluster.locations[, 1],
      no.muts = nrow(dataset$mutCount),
      no.timepoints = ncol(dataset$mutCount),
      no.iters = no.iters,
      no.iters.burn.in = no.iters.burn.in,
      num_threads = num_threads
    )
    conf <- data.frame(conf)
    colnames(conf) <- c("cluster.no", "timepoint", "loc_conf_0.025", "loc_conf_0.500", "loc_conf_0.975")
    write.table(conf, file = paste(outfiles.prefix, "_clusterConfidenceIntervals.txt", sep = ""), quote = FALSE, row.names = FALSE, sep = "\t")

    # Calc probs for cluster orders
    probs <- calc_cluster_order_probs(
      GS.data = GS.data,
      density = density,
      mut_assignments = clustering$best.node.assignments,
      clusterids = clustering$cluster.locations[, 1],
      cluster_ccfs = clustering$cluster.locations[, 2],
      no.muts = nrow(dataset$mutCount),
      no.timepoints = ncol(dataset$mutCount),
      no.iters = no.iters,
      no.iters.burn.in = no.iters.burn.in,
      no.samples = no.samples.cluster.order,
      num_threads = num_threads
    )
    probs <- flatten_3d_to_2d(probs$classification, c("timepoint", clustering$cluster.locations[, 1]))
    fwrite(probs, file = paste(outfiles.prefix, "_clusterOrderProbabilities.txt", sep = ""), quote = FALSE, row.names = FALSE, sep = "\t", na = "NA")
  }

  ########################################################################
  # Check if mutation sampling has been done, if so, unpack and assign here
  ########################################################################
  if (.has_value(most.similar.mut) && assign_sampled_muts) {
    res <- unsample_mutations(dataset, clustering)
    dataset <- res$dataset
    clustering <- res$clustering
  }

  ########################################################################
  # Write out the final mutation-cluster probabilities with all mutations spiked in
  ########################################################################
  cna.assignment.likelihoods <- NULL
  cluster_prob_colnames <- NULL
  if (.has_assignment_likelihoods(clustering)) {
    # Fetch and drop all columns that have just zeroes
    cols_all_zero <- which(apply(clustering$all.assignment.likelihoods, 2, max) == 0)
    if (length(cols_all_zero) != 0) {
      filtered.assignment.likelihoods <- clustering$all.assignment.likelihoods[, -cols_all_zero, drop = FALSE]
      cluster_colnames <- (1:ncol(clustering$all.assignment.likelihoods))[-cols_all_zero]
    } else {
      filtered.assignment.likelihoods <- clustering$all.assignment.likelihoods
      cluster_colnames <- 1:ncol(clustering$all.assignment.likelihoods)
    }

    cluster_prob_colnames <- paste("prob.cluster", cluster_colnames, sep = ".")
    snv_index <- dataset$mutationType == "SNV"
    dt_filtered <- as.data.table(filtered.assignment.likelihoods)
    snv_assignment_likelihoods <- cbind(
      as.data.table(dataset$chromosome[snv_index, 1]),
      dataset$position[snv_index, 1] - 1,
      dataset$position[snv_index, 1],
      dt_filtered[snv_index, ],
      clustering$best.node.assignments[snv_index]
    )
    setDT(snv_assignment_likelihoods)
    colnames(snv_assignment_likelihoods) <- c("chr", "start", "end", cluster_prob_colnames, "most.likely.cluster")
    fwrite(snv_assignment_likelihoods, file = paste(outfiles.prefix, "_mutationClusterLikelihoods.bed", sep = ""), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

    if (any(dataset$mutationType == "CNA")) {
      cna_index <- dataset$mutationType == "CNA"
      cna.assignment.likelihoods <- cbind(
        as.data.table(dataset$chromosome[cna_index, 1]),
        dataset$position[cna_index, 1] - 1,
        dataset$position[cna_index, 1],
        dt_filtered[cna_index, ],
        clustering$best.node.assignments[cna_index]
      )
      setDT(cna.assignment.likelihoods)
      colnames(cna.assignment.likelihoods) <- c("chr", "start", "end", cluster_prob_colnames, "most.likely.cluster")
      fwrite(cna.assignment.likelihoods, file = paste(outfiles.prefix, "_mutationClusterLikelihoodsPseudoSNV.bed", sep = ""), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
    }
  }

  ########################################################################
  # Write final cluster locations
  ########################################################################
  # Add the confidence intervals to the final cluster locations
  if (ncol(clustering$cluster.locations) > 3) {
    # nD based clustering
    write.table(clustering$cluster.locations, paste(outfiles.prefix, "_bestClusterInfo.txt", sep = ""), col.names = c("cluster.no", paste(samplename, subsamplenames, sep = ""), "no.of.mutations"), sep = "\t", quote = FALSE, row.names = FALSE)
  } else {
    # 1D based
    write.table(clustering$cluster.locations, paste(outfiles.prefix, "_bestClusterInfo.txt", sep = ""), col.names = c("cluster.no", "location", "no.of.mutations"), row.names = FALSE, sep = "\t", quote = FALSE)
  }

  ########################################################################
  # Add the removed mutations back in
  ########################################################################
  output <- cbind(dataset$chromosome[, 1], dataset$position[, 1] - 1, dataset$position[, 1], clustering$best.node.assignments, clustering$best.assignment.likelihoods)
  output <- add_removed_snvs(dataset, output)

  ########################################################################
  # Save the indices of the mutations that were not used during the analysis
  ########################################################################
  objects_to_save <- c("output", "clustering", "density", "dataset", "most.similar.mut", "no.iters", "no.iters.burn.in")
  if (!is.null(GS.data) && !.is_na_sentinel(GS.data)) {
    objects_to_save <- c(objects_to_save, "GS.data")
  }
  .parallel_save(objects_to_save, paste(outfiles.prefix, "_bestConsensusResults.RData", sep = ""), num_threads = num_threads)
  write.table(data.frame(mut.index = dataset$removed_indices), file = paste(outfiles.prefix, "_removedMutationsIndex.txt", sep = ""), row.names = FALSE, quote = FALSE)

  output <- as.data.table(output)
  colnames(output) <- c("chr", "start", "end", "cluster", "likelihood")

  # Construct a mutation type vector that matches the expanded output length
  num_orig_snvs <- length(dataset$chromosome.not.filtered)
  num_pseudo_snvs <- sum(dataset$mutationType != "SNV")
  full_mutation_type <- c(rep("SNV", num_orig_snvs), as.character(dataset$mutationType[dataset$mutationType != "SNV"]))
  
  fwrite(output[full_mutation_type == "SNV", ], file = paste(outfiles.prefix, "_bestConsensusAssignments.bed", sep = ""), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

  ########################################################################
  # Save the CNA assignments separately
  ########################################################################
  if (num_pseudo_snvs > 0) {
    fwrite(output[full_mutation_type != "SNV", ], file = paste(outfiles.prefix, "_bestConsensusAssignmentsPseudoSNV.bed", sep = ""), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
    # Assign the CNAs to clusters using their pseudoSNV representations
    cndata <- assign_cnas_to_clusters(dataset$cndata, output)
    fwrite(cndata, file = paste(outfiles.prefix, "_bestCNAassignments.txt", sep = ""), quote = FALSE, row.names = FALSE, sep = "\t", na = "NA")

    if (!is.null(cna.assignment.likelihoods)) {
      cna_assignment_likelihoods <- get_cnas_cluster_probs(as.data.table(dataset$cndata), as.data.table(cna.assignment.likelihoods), c("chr", "start", "end", cluster_prob_colnames, "most.likely.cluster"))
      fwrite(cna_assignment_likelihoods, file = paste(outfiles.prefix, "_cnaClusterLikelihoods.bed", sep = ""), quote = FALSE, row.names = FALSE, sep = "\t", na = "NA")
    }

    # Create a new assignment table figure with the correct information
    # This removes pseudo SNVs as the assignmentTable will add an extra column for CNAs
    cluster_locations <- clustering$cluster.locations
    cluster_locations[, 3] <- rep(0, nrow(cluster_locations))
    mut_assignments <- table(output[full_mutation_type == "SNV", cluster]) # cluster is a column name in data.table output
    for (i in seq_len(nrow(cluster_locations))) {
      if (as.character(cluster_locations[i, 1]) %in% names(mut_assignments)) {
        cluster_locations[i, 3] <- mut_assignments[as.character(cluster_locations[i, 1])]
      }
    }
    plotAssignmentTable(cluster_locations, paste(outfiles.prefix, "_mutation_assignments.png", sep = ""), cndata = cndata, num_samples = num_samples)
  } else {
    plotAssignmentTable(clustering$cluster.locations, paste(outfiles.prefix, "_mutation_assignments.png", sep = ""), num_samples = num_samples)
  }

  ########################################################################
  # If tree based analysis, also save the tree
  ########################################################################
  if (write_tree) {
    write.table(clustering$best.tree, file = paste(outfiles.prefix, "_bestConsensusTree.txt", sep = ""), quote = FALSE, row.names = FALSE, sep = "\t")
  }
}

#' Add removed mutations back into the assignment table. SNVs will be assigned to the cluster of its most similar not-removed SNV
#' @param dataset A dataset object
#' @param snv_assignment_table Data frame with the mutation assignments
#' @return The snv_assignment_table with the removed mutations added into the position they were originally
#' @author sd11
add_removed_snvs <- function(dataset, snv_assignment_table) {
  num_original_snvs <- length(dataset$chromosome.not.filtered)
  num_removed_snvs <- length(dataset$removed_indices)
  
  # Identify SNVs vs Pseudo-SNVs in the current table using logical indexing
  is_snv <- dataset$mutationType == "SNV"
  
  snvs_in_table <- snv_assignment_table[is_snv, , drop = FALSE]
  pseudos_in_table <- snv_assignment_table[!is_snv, , drop = FALSE]
  
  # Pre-allocate full SNV set
  full_snvs <- as.data.frame(matrix(NA, nrow = num_original_snvs, ncol = ncol(snv_assignment_table)))
  
  removed_idx <- dataset$removed_indices
  non_removed_idx <- setdiff(seq_len(num_original_snvs), removed_idx)
  
  # Fill non-removed SNVs
  full_snvs[non_removed_idx, ] <- snvs_in_table
  
  # Fill removed SNVs basic info
  if (num_removed_snvs > 0) {
    full_snvs[removed_idx, 1] <- dataset$chromosome.not.filtered[removed_idx]
    full_snvs[removed_idx, 2] <- dataset$mut.position.not.filtered[removed_idx] - 1
    full_snvs[removed_idx, 3] <- dataset$mut.position.not.filtered[removed_idx]
  }
  
  # Re-append pseudo-SNVs
  if (nrow(pseudos_in_table) > 0) {
    output <- rbind(full_snvs, pseudos_in_table)
  } else {
    output <- full_snvs
  }
  
  return(output)
}


#' Assign CNA events to clusters using their pseudoSNV representation
#' @param cndata Data frame with the CNA data
#' @param snv_assignment_table Data frame with the mutation assignments, with chromosome and position-start/end expected as first three columns
#' @return The cndata object with extra column cluster_assignment
#' @author sd11
assign_cnas_to_clusters <- function(cndata, snv_assignment_table) {
  dt_cna <- as.data.table(cndata)
  dt_snv <- as.data.table(snv_assignment_table)
  colnames(dt_snv) <- c("chr", "start", "end", "cluster", "likelihood")
  
  # Vote for majority cluster by chr and mutation position (snv end)
  # In pseudo-SNV conversion (load_data/__init_pseudo_snvs), the SNV end position matches cndata startpos
  res_dt <- dt_snv[, .(cluster_assignment = names(which.max(table(cluster)))), by = .(chr, end)]
  dt_cna <- merge(dt_cna, res_dt, by.x = c("chr", "startpos"), by.y = c("chr", "end"), all.x = TRUE)
  dt_cna[is.na(cluster_assignment), cluster_assignment := "NA"]
  
  return(as.data.frame(dt_cna))
}

#' Use the Pseudo-SNV probabilities to obtain a probability of each CNA of each cluster
#' @param cndata The copy number data data.frame
#' @param snv_assignment_likelihoods Probabilities of the pseudo-SNVs
#' @param cluster_colnames The colnames of the output bedfile to be used to select the assignment columns
#' @return A data.frame with chr,start,end,probs_per_cluster
#' @author sd11
get_cnas_cluster_probs <- function(cndata, snv_assignment_likelihoods, cluster_colnames) {
  dt_cna <- as.data.table(cndata)
  dt_snv <- as.data.table(snv_assignment_likelihoods)
  
  prob_cols <- setdiff(colnames(dt_snv), c("chr", "start", "end", "most.likely.cluster"))
  
  # Group pseudo-SNVs by their representative CNA (chr/startpos) and average probabilities
  # In pseudo-SNV conversion, the end position of the pseudo-SNV matches the startpos of the CNA segment.
  res_dt <- dt_snv[, lapply(.SD, mean), by = .(chr, end), .SDcols = prob_cols]
  dt_cna <- merge(dt_cna, res_dt, by.x = c("chr", "startpos"), by.y = c("chr", "end"), all.x = TRUE)
  
  # Identify most likely cluster from averaged percentages
  # (Since result might have many columns, using which.max across the probability columns)
  if (nrow(dt_cna) > 0) {
    # Extract just the prob columns to find max
    prob_subset <- as.matrix(dt_cna[, ..prob_cols])
    dt_cna$most.likely.cluster <- apply(prob_subset, 1, function(x) {
        if (all(is.na(x))) return(NA)
        return(which.max(x))
    })
  }

  return(dt_cna)
}

#' Helper function that flattens a 3D array into a 2D one
#' @param data The data to be flattened
#' @param col_names The names of the columns in the output
#' @return A data.frame with the third column annotated as the first column
#' @author sd11
flatten_3d_to_2d <- function(data, col_names) {
  no.timepoints <- dim(data)[3]
  no.clusters <- dim(data)[2]
  new_data <- data.frame(array(NA, c(no.clusters * no.timepoints, no.clusters + 1)))
  for (i in 1:no.timepoints) {
    row <- ((i - 1) * no.clusters) + 1
    new_data[row:(row + no.clusters - 1), 2:(no.clusters + 1)] <- data[, , i]
    new_data[row:(row + no.clusters - 1), 1] <- i
  }
  colnames(new_data) <- col_names
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
DirichletProcessClustering <- function(mutCount, WTCount, totalCopyNumber, copyNumberAdjustment, mutation.copy.number, cellularity, output_folder, no.iters, no.iters.burn.in, subsamplesrun, samplename, conc_param, cluster_conc, mut.assignment.type, most.similar.mut, mutationTypes, max.considered.clusters, thin_s_i = FALSE, keep_aux_fields = FALSE, num_threads = NA_integer_, conflict.array = .init_conflicts(), keep_temp_files = TRUE, density_smooth = NA_real_, x_max_cap = 3, hypercube_size = 5) {
  output_folder <- normalizePath(output_folder, mustWork = FALSE)
  if (!dir.exists(output_folder)) {
    dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)
  }
  density_smooth_nd <- if (is.na(density_smooth)) 0.01 else density_smooth
  density_smooth_1d <- if (is.na(density_smooth)) 0.1 else density_smooth
  stored_iters <- integer(0)
  if (thin_s_i) {
    stored_iters <- (no.iters.burn.in + 1):no.iters
    stored_iters <- stored_iters[stored_iters != 1]
    if (length(stored_iters) > 1000) {
      stored_iters <- floor(no.iters.burn.in + (1:1000) * (no.iters - no.iters.burn.in) / 1000)
    }
    stored_iters <- as.integer(sort(unique(stored_iters)))
    log_info(paste("Memory optimization: storing", length(stored_iters), "state iterations instead of", no.iters))
  }
  log_info("Entering Gibbs sampler...")
  GS.data <- subclone.dirichlet.gibbs(
    mutCount = mutCount,
    WTCount = WTCount,
    totalCopyNumber = totalCopyNumber,
    copyNumberAdjustment = copyNumberAdjustment,
    cellularity = cellularity,
    iter = no.iters,
    conc_param = conc_param,
    cluster_conc = cluster_conc,
    C = max.considered.clusters,
    keep_aux_fields = keep_aux_fields,
    num_threads = num_threads,
    stored_iters = stored_iters,
    conflict.array = conflict.array
  )

  if (keep_temp_files) {
    .parallel_save("GS.data", file.path(output_folder, paste(samplename, "_gsdata.RData", sep = "")), num_threads = num_threads)
  }

  # nD dataset, plot sample versus sample
  if (ncol(mutCount) > 1) {
    ########################
    # Plot density and Assign mutations to clusters - nD
    ########################
    log_info("Estimating density between pairs of samples...")
    for (i in 1:(length(subsamplesrun) - 1)) {
      for (j in (i + 1):length(subsamplesrun)) {
        log_info(paste("Samples", subsamplesrun[i], "and", subsamplesrun[j], sep = " "))
        imageFile <- file.path(output_folder, paste(samplename, subsamplesrun[i], subsamplesrun[j], "_iters", no.iters, "_concParam", conc_param, "_clusterWidth", 1 / cluster_conc, "_2D_binomial.png", sep = ""))
        density <- Gibbs.subclone.density.est(mutation.copy.number[, c(i, j)] / copyNumberAdjustment[, c(i, j)],
          GS.data,
          imageFile,
          post.burn.in.start = no.iters.burn.in,
          post.burn.in.stop = no.iters,
          samplenames = paste(samplename, subsamplesrun[c(i, j)], sep = ""),
          indices = c(i, j),
          density.smooth = density_smooth_nd
        )
        save(file = file.path(output_folder, paste(samplename, subsamplesrun[i], subsamplesrun[j], "_densityoutput.RData", sep = "")), density)
      }
    }

    # Assign mutations to clusters using one of the different assignment methods
    log_info("Assigning mutations to clusters...")
    opts <- list(samplename = samplename, subsamplenames = subsamplesrun, no.iters = no.iters, no.iters.burn.in = no.iters.burn.in, no.iters.post.burn.in = no.iters - no.iters.burn.in, outdir = output_folder)
    if (mut.assignment.type == 1) {
      consClustering <- multiDimensionalClustering(
        mutation.copy.number = mutation.copy.number,
        copyNumberAdjustment = copyNumberAdjustment,
        GS.data = GS.data,
        density.smooth = density_smooth_nd,
        opts = opts,
        num_threads = num_threads
      )
    } else if (mut.assignment.type == 2) {
      consClustering <- mutation_assignment_em(
        GS.data = GS.data,
        mutCount = mutCount,
        WTCount = WTCount,
        subclonal.fraction = mutation.copy.number / copyNumberAdjustment,
        node.assignments = GS.data$S.i,
        opts = opts
      )
    } else if (mut.assignment.type == 3) {
      warning("binom mut assignment not implemented for multiple timepoints")
      stop("binom mut assignment not implemented for multiple timepoints")
    } else {
      warning(paste("Unknown mutation assignment type", mut.assignment.type, sep = " "))
      stop(paste("Unknown mutation assignment type", mut.assignment.type, sep = " "))
    }
    consClustering$GS.data <- GS.data
    return(consClustering)

    # 1D dataset, plot just the single density
  } else {
    ########################
    # Plot density and Assign mutations to clusters - 1D
    ########################
    # 1D dataset, plot just the single density
    # Use absolute path - avoids setwd() which breaks inside Singularity containers
    res <- Gibbs.subclone.density.est.1d(GS.data,
      file.path(output_folder, paste0(samplename, "_DirichletProcessplot.png")),
      samplename = samplename,
      post.burn.in.start = no.iters.burn.in,
      post.burn.in.stop = no.iters,
      y.max = 15,
      x.max = NA,
      x.max.cap = NA,
      mutationCopyNumber = mutation.copy.number,
      no.chrs.bearing.mut = copyNumberAdjustment,
      density.smooth = density_smooth_1d
    )
    density <- res$density
    polygon.data <- res$polygon.data

    # Assign mutations to clusters using one of the different assignment methods
    opts <- list(samplename = samplename, subsamplenames = subsamplesrun, no.iters = no.iters, no.iters.burn.in = no.iters.burn.in, no.iters.post.burn.in = no.iters - no.iters.burn.in, outdir = output_folder)

    log_info("Assigning mutations to clusters...")
    if (mut.assignment.type == 1) {
      subclonal.fraction <- mutation.copy.number / copyNumberAdjustment
      subclonal.fraction[is.nan(subclonal.fraction)] <- 0
      consClustering <- oneDimensionalClustering(samplename, subclonal.fraction, GS.data, density, no.iters, no.iters.burn.in, outdir = output_folder, num_threads = num_threads, hypercube.size = hypercube_size)
    } else if (mut.assignment.type == 2) {
      consClustering <- mutation_assignment_em(
        GS.data = GS.data,
        mutCount = mutCount,
        WTCount = WTCount,
        subclonal.fraction = mutation.copy.number / copyNumberAdjustment,
        node.assignments = GS.data$S.i,
        opts = opts
      )
    } else if (mut.assignment.type == 3) {
      consClustering <- mutation_assignment_binom(
        clustering_density = density,
        mutCount = mutCount,
        WTCount = WTCount,
        copyNumberAdjustment = copyNumberAdjustment,
        tumourCopyNumber = totalCopyNumber,
        normalCopyNumber = array(2, dim(mutCount)),
        cellularity = cellularity,
        samplename = samplename,
        outdir = output_folder
      )
    } else {
      warning(paste("Unknown mutation assignment type", mut.assignment.type, sep = " "))
      stop(paste("Unknown mutation assignment type", mut.assignment.type, sep = " "))
    }

    # Make a second set of figures with the mutation assignments showing
    # Replot the data with cluster locations
    plot1D(
      density = density,
      polygon.data = polygon.data,
      pngFile = file.path(output_folder, paste(samplename, "_DirichletProcessplot_with_cluster_locations.png", sep = "")),
      density.from = 0,
      x.max = NA,
      x.max.cap = x_max_cap,
      mutationCopyNumber = mutation.copy.number,
      no.chrs.bearing.mut = copyNumberAdjustment,
      samplename = samplename,
      cluster.locations = consClustering$cluster.locations,
      mutation.assignments = consClustering$best.node.assignments
    )

    plot1D_2(
      density = density,
      polygon.data = polygon.data,
      pngFile = file.path(output_folder, paste(samplename, "_DirichletProcessplot_with_cluster_locations_2.png", sep = "")),
      density.from = 0,
      x.max = NA,
      x.max.cap = x_max_cap,
      mutationCopyNumber = mutation.copy.number,
      no.chrs.bearing.mut = copyNumberAdjustment,
      samplename = samplename,
      cluster.locations = consClustering$cluster.locations,
      mutation.assignments = consClustering$best.node.assignments,
      mutationTypes = mutationTypes
    )

    consClustering$GS.data <- GS.data
    return(consClustering)
  }
}

write.strengths.table <- function(dat, removed_indices, filename) {
  #
  # Adds in an empty column/row for each of the removed_indices and writes
  # the subsequent matrix to disk
  #
  write.table(add.muts.back.in(dat, removed_indices), filename, sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)
}

add.muts.back.in <- function(dat, removed_indices, def.value = 0) {
  #
  # Adds in empty columns and rows for mutations that were removed.
  # Mutations are added sequentially, so this method expects indices
  # of removed mutations in the original (full) matrix. The empty
  # mutations will be added in the place where they were removed,
  # keeping the order in tact.
  #
  for (i in removed_indices) {
    if (i == 1) {
      dat <- cbind(rep(def.value, nrow(dat)), dat)
      dat <- rbind(rep(def.value, ncol(dat)), dat)
    } else if (i >= ncol(dat)) {
      dat <- cbind(dat, rep(def.value, nrow(dat)))
      dat <- rbind(dat, rep(def.value, ncol(dat)))
    } else {
      dat <- cbind(dat[, 1:(i - 1)], rep(def.value, nrow(dat)), dat[, i:ncol(dat)])
      dat <- rbind(dat[1:(i - 1), ], rep(def.value, ncol(dat)), dat[i:nrow(dat), ])
    }
  }
  return(dat)
}
