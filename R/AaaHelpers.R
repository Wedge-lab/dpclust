# Small internal helpers to detect placeholder NA fields in dataset/clustering objects.
.is_na_sentinel <- function(x) {
  if (is.null(x) || length(x) == 0) {
    return(TRUE)
  }
  if (is.atomic(x) && length(x) == 1 && is.na(x)) {
    return(TRUE)
  }
  FALSE
}

.has_value <- function(x) {
  !.is_na_sentinel(x)
}

.init_conflicts <- function() {
  list(i = integer(0), j = integer(0), w = numeric(0))
}

.subset_conflicts <- function(conflicts, selection) {
  if (!.has_value(conflicts) || length(conflicts$i) == 0) {
    return(conflicts)
  }
  new_i <- match(conflicts$i, selection)
  new_j <- match(conflicts$j, selection)
  keep <- !is.na(new_i) & !is.na(new_j)
  list(
    i = as.integer(new_i[keep]),
    j = as.integer(new_j[keep]),
    w = as.numeric(conflicts$w[keep])
  )
}

.has_assignment_likelihoods <- function(clustering) {
  if (!("all.assignment.likelihoods" %in% names(clustering))) {
    return(FALSE)
  }
  all.assignment.likelihoods <- clustering$all.assignment.likelihoods
  if (!.has_value(all.assignment.likelihoods)) {
    return(FALSE)
  }
  if (is.null(dim(all.assignment.likelihoods))) {
    return(FALSE)
  }
  nrow(all.assignment.likelihoods) > 0 && ncol(all.assignment.likelihoods) > 0
}

.dpclust_num_stored_iters <- function(no.iters, no.iters.burn.in, thin_s_i) {
  if (!thin_s_i) {
    return(as.integer(no.iters))
  }
  stored_iters <- (no.iters.burn.in + 1):no.iters
  stored_iters <- stored_iters[stored_iters != 1]
  if (length(stored_iters) > 1000) {
    stored_iters <- floor(no.iters.burn.in + (1:1000) * (no.iters - no.iters.burn.in) / 1000)
  }
  as.integer(length(unique(stored_iters)))
}

.estimate_dpclust_memory_gb <- function(no.muts, no.samples, no.iters, no.iters.burn.in, max.considered.clusters, thin_s_i, keep_aux_fields) {
  bytes_per_double <- 8
  bytes_per_int <- 4
  stored_iters <- .dpclust_num_stored_iters(no.iters, no.iters.burn.in, thin_s_i)

  # Base R overhead + DPClust and dependency loading (approx 250 MB)
  bytes_base <- 250 * 1024^2

  bytes_pi_h <- as.numeric(no.iters) * as.numeric(max.considered.clusters) * as.numeric(no.samples) * bytes_per_double
  bytes_v_h <- as.numeric(no.iters) * as.numeric(max.considered.clusters) * bytes_per_double
  bytes_alpha <- as.numeric(no.iters) * bytes_per_double
  bytes_s_i <- as.numeric(stored_iters) * as.numeric(no.muts) * bytes_per_int

  # The dataset object holds ~10 matrices of doubles [no.muts x no.samples]
  # (WTCount, mutCount, totalCopyNumber, copyNumberAdjustment, kappa, chromosome, etc)
  # Plus R data.frame overhead.
  bytes_dataset <- 12 * as.numeric(no.muts) * as.numeric(no.samples) * bytes_per_double

  # Mutation preferences matrix created during assignment (assume ~30 potential clusters)
  bytes_preferences <- as.numeric(no.muts) * 30 * bytes_per_double

  # Buffer for temporary copies and GC margin (20% of the major objects)
  bytes_buffer <- (bytes_s_i + bytes_dataset) * 0.2

  bytes_aux <- 0
  if (keep_aux_fields) {
    bytes_aux <- 2 * as.numeric(no.muts) * as.numeric(no.samples) * bytes_per_double
  }

  bytes_total <- bytes_base + bytes_pi_h + bytes_v_h + bytes_alpha + bytes_s_i + bytes_dataset + bytes_preferences + bytes_buffer + bytes_aux
  gb_raw <- bytes_total / (1024^3)

  list(
    total_gb = gb_raw,
    raw_gb = gb_raw,
    stored_iters = stored_iters,
    components_gb = c(
      base = bytes_base / (1024^3),
      pi_h = bytes_pi_h / (1024^3),
      S_i = bytes_s_i / (1024^3),
      dataset = bytes_dataset / (1024^3),
      preferences = bytes_preferences / (1024^3),
      buffer = bytes_buffer / (1024^3),
      aux = bytes_aux / (1024^3)
    )
  )
}

.memory_guard_plan <- function(no.muts, no.samples, no.iters, no.iters.burn.in, max.considered.clusters, thin_s_i, keep_aux_fields, memory_limit_gb, verbose = TRUE) {
  estimate <- .estimate_dpclust_memory_gb(
    no.muts = no.muts,
    no.samples = no.samples,
    no.iters = no.iters,
    no.iters.burn.in = no.iters.burn.in,
    max.considered.clusters = max.considered.clusters,
    thin_s_i = thin_s_i,
    keep_aux_fields = keep_aux_fields
  )

  if (is.na(memory_limit_gb)) {
    return(list(thin_s_i = thin_s_i, keep_aux_fields = keep_aux_fields, estimate = estimate, limited = FALSE))
  }
  if (!is.finite(memory_limit_gb) || memory_limit_gb <= 0) {
    stop("memory_limit_gb must be a positive finite number when provided.")
  }

  if (verbose) {
    log_info(sprintf("Memory guard estimate: %.2f GB (limit %.2f GB)", estimate$total_gb, memory_limit_gb))
  }

  if (estimate$total_gb <= memory_limit_gb) {
    return(list(thin_s_i = thin_s_i, keep_aux_fields = keep_aux_fields, estimate = estimate, limited = FALSE))
  }

  detail_msg <- sprintf(
    paste0(
      "Estimated peak memory %.2f GB exceeds memory_limit_gb %.2f GB. ",
      "Estimated components (GB): base=%.2f, pi.h=%.2f, S.i=%.2f, dataset=%.2f, preferences=%.2f, buffer=%.2f. ",
      "Try reducing no.iters, max.considered.clusters, or num_muts_sample."
    ),
    estimate$total_gb, memory_limit_gb,
    estimate$components_gb["base"], estimate$components_gb["pi_h"], estimate$components_gb["S_i"],
    estimate$components_gb["dataset"], estimate$components_gb["preferences"], estimate$components_gb["buffer"]
  )
  stop(detail_msg)
}
# High-performance parallel save that preserves standard .RData compatibility.
# Uses a pipe to pigz (Parallel Implementation of GZip) if available.
.parallel_save <- function(list_to_save, file_name, num_threads = NA_integer_) {
  if (is.na(num_threads)) {
    num_threads <- 4
  }

  # Priority list of parallel compression tools
  parallel_tools <- list(
    list(cmd = "pigz", args = "-6 -p %d"),
    list(cmd = "bgzip", args = "-@ %d -c")
  )

  selected_tool <- NULL
  for (tool in parallel_tools) {
    # Sys.which is the standard R way to find an executable in the PATH.
    found_path <- Sys.which(tool$cmd)
    if (found_path != "" && file.access(found_path, 1) == 0) {
      selected_tool <- tool
      selected_tool$full_path <- as.character(found_path)
      break
    }
  }

  # Capture the environment of the caller
  caller_env <- parent.frame()

  if (!is.null(selected_tool)) {
    # Parallel hack: R outputs binary to a pipe, tool compresses it using many cores.
    log_info(sprintf("Parallelizing .RData compression using %s (%d threads)...", selected_tool$cmd, num_threads))
    log_debug(sprintf("Using parallel tool at: %s", selected_tool$full_path))

    # Standard balance between speed and size (Gzip level 6)
    # We use the full path to avoid any ambiguity during execution.
    con_cmd <- sprintf("%s %s > %s", selected_tool$full_path, sprintf(selected_tool$args, num_threads), file_name)
    con <- pipe(con_cmd, "wb")

    on.exit(close(con))
    save(list = list_to_save, file = con, compress = FALSE, envir = caller_env)
  } else {
    # Fallback: standard synchronous compression.
    log_info("Parallel compression tool (pigz/bgzip) not found in PATH.")
    log_info(sprintf("Current PATH: %s", Sys.getenv("PATH")))
    log_info("Falling back to standard synchronous (slow) compression...")
    save(list = list_to_save, file = file_name, compress = TRUE, envir = caller_env)
  }
}
