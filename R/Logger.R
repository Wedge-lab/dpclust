#' DPClust Logging Helpers
#' @keywords internal
log_setup <- function(path, verbose) {
  options(dpclust.verbose = verbose)
  if (!is.null(path) && path != "") {
    # If it's a directory, create a default log file name
    if (dir.exists(path) || !grepl("\\.log$", path)) {
      if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
      path <- file.path(path, "dpclust_pipeline.log")
    } else {
      # If it's a file path, ensure parent dir exists
      log_dir <- dirname(path)
      if (!dir.exists(log_dir)) dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
    }
    options(dpclust.log_file = path)
    # Clear existing log or start fresh
    cat(paste0("--- DPClust Log Started: ", Sys.time(), " ---\n"), file = path, append = FALSE)
  }
}

#' @keywords internal
#' @export
log_info <- function(msg) {
  formatted_msg <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%OS6"), " | ", msg)
  message(formatted_msg)
  log_file <- getOption("dpclust.log_file")
  if (!is.null(log_file)) {
    cat(formatted_msg, "\n", file = log_file, append = TRUE)
  }
}

#' @keywords internal
#' @export
log_failure <- function(msg) {
  formatted_msg <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%OS6"), " | Fatal Error: ", msg)
  message(formatted_msg)
  log_file <- getOption("dpclust.log_file")
  if (!is.null(log_file)) {
    cat(formatted_msg, "\n", file = log_file, append = TRUE)
  }
}

#' @keywords internal
#' @export
log_debug <- function(msg) {
  if (isTRUE(getOption("dpclust.verbose"))) {
    formatted_msg <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%OS6"), " | [DEBUG] ", msg)
    message(formatted_msg)
    log_file <- getOption("dpclust.log_file")
    if (!is.null(log_file)) {
      cat(formatted_msg, "\n", file = log_file, append = TRUE)
    }
  }
}
