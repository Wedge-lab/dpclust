#
# DPClust core algorithm
#

subclone.dirichlet.gibbs <- function(mutCount, WTCount, totalCopyNumber = array(1, dim(mutCount)), normalCopyNumber = array(2, dim(mutCount)), copyNumberAdjustment = array(1, dim(mutCount)), C = 30, cellularity = rep(1, ncol(mutCount)), iter = 1000, conc_param = 1, cluster_conc = 10, keep_aux_fields = FALSE, num_threads = NA_integer_, stored_iters = integer(0), conflict.array = .init_conflicts()) {
  if (is.null(copyNumberAdjustment)) {
    copyNumberAdjustment <- array(1, dim(mutCount))
  }

  log_info("Converting data to matrices for Rcpp...")
  ensure_numeric_matrix <- function(x) {
    if (!is.matrix(x)) {
      x <- as.matrix(x)
    }
    if (!is.double(x)) {
      storage.mode(x) <- "double"
    }
    x
  }
  # Ensure inputs are in the correct format for Rcpp while avoiding unnecessary copies.
  mutCount <- ensure_numeric_matrix(mutCount)
  WTCount <- ensure_numeric_matrix(WTCount)
  totalCopyNumber <- ensure_numeric_matrix(totalCopyNumber)
  normalCopyNumber <- ensure_numeric_matrix(normalCopyNumber)
  copyNumberAdjustment <- ensure_numeric_matrix(copyNumberAdjustment)
  if (!is.double(cellularity)) {
    storage.mode(cellularity) <- "double"
  }

  log_info(paste("Calling Rcpp Gibbs sampler for", nrow(mutCount), "mutations and", iter, "iterations..."))
  # Call C++ implementation
  cpp_threads <- if (is.na(num_threads)) as.integer(-1) else as.integer(num_threads)

  # Report the actual thread situation before launching so it's visible in the log
  omp_available <- isTRUE(tryCatch(DPClust:::omp_thread_count_cpp() > 0, error = function(e) FALSE))
  if (omp_available) {
    actual_threads <- DPClust:::omp_thread_count_cpp(cpp_threads)
    log_info(paste("OpenMP: ENABLED — Gibbs sampler will use", actual_threads, "thread(s)"))
  } else {
    log_info("OpenMP: NOT AVAILABLE in this build — running single-threaded (recompile with OpenMP support for parallelism)")
  }

  res <- subclone_dirichlet_gibbs_cpp(mutCount, WTCount, totalCopyNumber, normalCopyNumber, copyNumberAdjustment, C, cellularity, iter, conc_param, cluster_conc, keep_aux_fields, cpp_threads, as.integer(stored_iters), as.integer(conflict.array$i), as.integer(conflict.array$j), as.numeric(conflict.array$w), log_info)
  if (!keep_aux_fields) {
    # mutBurdens is large and not used by the current downstream pipeline.
    res$mutBurdens <- NULL
  }
  log_info("Gibbs sampler completed.")

  return(res)
}
