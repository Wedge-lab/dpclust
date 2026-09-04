#
# Convenience functions for the pipelie
#

replot_1D <- function(outdir, outfiles.prefix, samplename, dataset, clustering, density, polygon.data, x_max_cap = 3) {
  # Old plot
  plot1D(
    density = density,
    polygon.data = polygon.data[, 1],
    pngFile = paste(outdir, "/", samplename, "_DirichletProcessplot_with_cluster_locations_replot.png", sep = ""),
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
    pngFile = paste(outdir, "/", samplename, "_DirichletProcessplot_with_cluster_locations_2_replot.png", sep = ""),
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

reassign_1D <- function(outdir, samplename, no.iters, no.iters.burn.in, dataset, cellularity, GS.data, conc_param, cluster_conc, mut.assignment.type, x_max_cap = 3) {
  # Use absolute path so this works both locally and inside containers
  outdir <- normalizePath(outdir, mustWork = FALSE)

  # Obtain the density using the absolute path for the PNG
  res <- Gibbs.subclone.density.est.1d(GS.data,
    file.path(outdir, paste0(samplename, "_DirichletProcessplot.png")),
    samplename = samplename,
    post.burn.in.start = no.iters.burn.in,
    post.burn.in.stop = no.iters,
    y.max = 15,
    x.max = NA,
    x.max.cap = NA,
    mutationCopyNumber = dataset$mutation.copy.number,
    no.chrs.bearing.mut = dataset$copyNumberAdjustment
  )
  clustering_density <- res$density
  polygon.data <- res$polygon.data
  opts <- list(samplename = samplename, subsamplenames = NA, no.iters = no.iters, no.iters.burn.in = no.iters.burn.in, no.iters.post.burn.in = no.iters - no.iters.burn.in, outdir = outdir)

  # Use one of the three ways of assigning mutations
  if (mut.assignment.type == 1) {
    subclonal.fraction <- dataset$mutation.copy.number / dataset$copyNumberAdjustment
    subclonal.fraction[is.nan(subclonal.fraction)] <- 0
    clustering <- oneDimensionalClustering(samplename, subclonal.fraction, GS.data, clustering_density, no.iters, no.iters.burn.in, outdir = outdir)
    # Change the outfiles prefix to be able to identify the output files
    outfiles.prefix <- file.path(outdir, paste0(samplename, "_reassign_option_1"))
  } else if (mut.assignment.type == 2) {
    clustering <- mutation_assignment_em(mutCount = dataset$mutCount, WTCount = dataset$WTCount, node.assignments = GS.data$S.i, opts = opts)
    # Change the outfiles prefix to be able to identify the output files
    outfiles.prefix <- file.path(outdir, paste0(samplename, "_reassign_option_2"))
  } else if (mut.assignment.type == 3) {
    clustering <- mutation_assignment_binom(
      clustering_density = clustering_density,
      mutCount = dataset$mutCount,
      WTCount = dataset$WTCount,
      copyNumberAdjustment = dataset$copyNumberAdjustment,
      tumourCopyNumber = dataset$totalCopyNumber,
      normalCopyNumber = array(2, dim(dataset$mutCount)),
      cellularity = cellularity,
      samplename = samplename
    )
    # Change the outfiles prefix to be able to identify the output files
    outfiles.prefix <- file.path(outdir, paste0(samplename, "_reassign_option_3"))
  } else {
    log_info(paste("Unknown mutation assignment type", mut.assignment.type, sep = " "))
    stop(paste("Unknown mutation assignment type", mut.assignment.type, sep = " "))
  }

  # Replot the data with cluster locations
  plot1D_2(
    density = clustering_density,
    polygon.data = polygon.data,
    pngFile = paste0(outfiles.prefix, "_DirichletProcessplot_with_cluster_locations_2.png"),
    density.from = 0,
    x.max = NA,
    x.max.cap = x_max_cap,
    mutationCopyNumber = dataset$mutation.copy.number,
    no.chrs.bearing.mut = dataset$copyNumberAdjustment,
    mutationTypes = dataset$mutationType,
    samplename = samplename,
    cluster.locations = clustering$cluster.locations,
    mutation.assignments = clustering$best.node.assignments
  )

  return(list(clustering = clustering, outfiles.prefix = outfiles.prefix))
}


replot_nD <- function(outdir, outfiles.prefix, samplename, subsamples, dataset, no.iters, clustering, conc_param, cluster_conc) {
  for (i in 1:(length(subsamples) - 1)) {
    for (j in (i + 1):length(subsamples)) {
      filename_prefix <- file.path(outdir, paste0(samplename, subsamples[i], subsamples[j], "_iters", no.iters, "_concParam", conc_param, "_clusterWidth", 1 / cluster_conc))
      pngFile <- paste0(filename_prefix, "_2D_binomial_with_cluster_locations_replot.png")
      xvals <- read.table(paste(filename_prefix, "_2D_binomial_xvals.csv", sep = ""), header = TRUE, sep = ",", row.names = 1)
      yvals <- read.table(paste(filename_prefix, "_2D_binomial_yvals.csv", sep = ""), header = TRUE, sep = ",", row.names = 1)
      zvals <- read.table(paste(filename_prefix, "_2D_binomial_zvals.csv", sep = ""), header = TRUE, sep = ",", row.names = 1)

      plotnD(
        xvals = xvals,
        yvals = yvals,
        zvals = zvals,
        subclonal.fraction_x = dataset$subclonal.fraction[, i],
        subclonal.fraction_y = dataset$subclonal.fraction[, j],
        pngFile = pngFile,
        samplename_x = paste(samplename, subsamples[i], sep = ""),
        samplename_y = paste(samplename, subsamples[j], sep = ""),
        max.plotted.value = 1.4,
        cluster.locations = clustering$cluster.locations[, c(i + 1, j + 1)]
      ) # +1 because the first column contains the cluster numbers
    }
  }
}

reassign_nd <- function(outdir, samplename, subsamples, no.iters, no.iters.burn.in, dataset, GS.data, conc_param, cluster_conc, mut.assignment.type) {
  # Load the output from the algorithm
  # GS.data = read_gsdata_object(outdir, no.iters=no.iters, conc_param=conc_param, cluster_conc=cluster_conc)
  # load(file=file.path(outdir, paste(samplename, "_gsdata.RData", sep="")))

  # Assign mutations to clusters using one of the different assignment methods
  opts <- list(samplename = samplename, subsamplenames = subsamples, no.iters = no.iters, no.iters.burn.in = no.iters.burn.in, no.iters.post.burn.in = no.iters - no.iters.burn.in, outdir = outdir)
  log_info("Assigning mutations to clusters")
  if (mut.assignment.type == 1) {
    clustering <- multiDimensionalClustering(
      mutation.copy.number = dataset$mutation.copy.number,
      copyNumberAdjustment = dataset$copyNumberAdjustment,
      GS.data = GS.data,
      density.smooth = 0.01,
      opts = opts
    )
    # Change the outfiles prefix to be able to identify the output files
    outfiles.prefix <- file.path(outdir, paste(samplename, "_reassign_option_1", sep = ""))
  } else if (mut.assignment.type == 2) {
    clustering <- mutation_assignment_em(mutCount = dataset$mutCount, WTCount = dataset$WTCount, node.assignments = GS.data$S.i, opts = opts)
    # Change the outfiles prefix to be able to identify the output files
    outfiles.prefix <- file.path(outdir, paste(samplename, "_reassign_option_2", sep = ""))
  } else if (mut.assignment.type == 3) {
    warning("binom mut assignment not implemented for multiple timepoints")
    stop("binom mut assignment not implemented for multiple timepoints")
  } else {
    log_info(paste("Unknown mutation assignment type", mut.assignment.type, sep = " "))
    stop(paste("Unknown mutation assignment type", mut.assignment.type, sep = " "))
  }

  # Make the figure
  replot_nD(outdir, outfiles.prefix, samplename, subsamples, dataset, clustering, conc_param, cluster_conc)

  return(list(clustering = clustering, outfiles.prefix = outfiles.prefix))
}
