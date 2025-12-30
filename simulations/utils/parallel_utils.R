# ============================================================================
# Parallel Processing Utilities
# ============================================================================
# Functions for:
#   - Setting up parallel processing with future package
#   - Generating reproducible seed sequences
#   - Splitting work across workers
# ============================================================================

#' Setup Parallel Processing
#'
#' Configures future plan for parallel execution
#'
#' @param strategy Character. future plan strategy (default: "multisession")
#'   Options: "sequential", "multisession", "multicore" (Unix only)
#' @param workers Integer. Number of parallel workers (default: detectCores() - 1)
#' @param verbose Logical. Print setup information
#' @return Integer. Number of workers configured
#' @details
#'   Uses future package for parallel processing.
#'   multisession: Works on all OS, uses separate R sessions
#'   multicore: Unix only, uses forking (more efficient but limited)
#'   sequential: No parallelization (for debugging)
#' @examples
#'   n_workers <- setup_parallel("multisession", workers = 4)
#'   # Now use future.apply functions for parallel execution
setup_parallel <- function(strategy = "multisession", workers = NULL, verbose = TRUE) {

  if (!requireNamespace("future", quietly = TRUE)) {
    stop("Package 'future' required for parallel processing")
  }

  # Determine number of workers
  if (is.null(workers)) {
    workers <- max(1, future::availableCores() - 1)
  }

  # Check strategy validity
  valid_strategies <- c("sequential", "multisession", "multicore")
  if (!strategy %in% valid_strategies) {
    stop(paste("Strategy must be one of:", paste(valid_strategies, collapse = ", ")))
  }

  # Warn about multicore on Windows
  if (strategy == "multicore" && .Platform$OS.type == "windows") {
    warning("multicore not available on Windows; switching to multisession")
    strategy <- "multisession"
  }

  # Set future plan
  if (strategy == "sequential") {
    future::plan(future::sequential)
    actual_workers <- 1
  } else if (strategy == "multisession") {
    future::plan(future::multisession, workers = workers)
    actual_workers <- workers
  } else if (strategy == "multicore") {
    future::plan(future::multicore, workers = workers)
    actual_workers <- workers
  }

  # Print information
  if (verbose) {
    cat("Parallel processing configured:\n")
    cat(sprintf("  Strategy: %s\n", strategy))
    cat(sprintf("  Workers: %d\n", actual_workers))
    cat(sprintf("  Available cores: %d\n", future::availableCores()))
  }

  return(invisible(actual_workers))
}


#' Generate Reproducible Seed Sequence
#'
#' Creates a sequence of random seeds for parallel execution
#'
#' @param n_seeds Integer. Number of seeds to generate
#' @param base_seed Integer. Base seed for reproducibility (optional)
#' @return Integer vector of length n_seeds
#' @details
#'   Ensures reproducibility even when tasks are run in parallel.
#'   Each parallel task gets a unique seed from this sequence.
#'   If base_seed is NULL, uses current RNG state.
#' @examples
#'   seeds <- generate_seeds(1000, base_seed = 12345)
#'   # Use seeds[i] for i-th replication
generate_seeds <- function(n_seeds, base_seed = NULL) {

  if (n_seeds < 1) stop("n_seeds must be at least 1")

  # Set base seed if provided
  if (!is.null(base_seed)) {
    set.seed(base_seed)
  }

  # Generate sequence of seeds
  # Use large range to avoid collisions
  seeds <- sample.int(.Machine$integer.max, size = n_seeds, replace = FALSE)

  return(seeds)
}


#' Split Work Across Workers
#'
#' Distributes task indices across parallel workers
#'
#' @param total_tasks Integer. Total number of tasks to complete
#' @param n_workers Integer. Number of parallel workers
#' @return List of length n_workers, each element is a vector of task indices
#' @details
#'   Balances workload across workers as evenly as possible.
#'   Useful for manual parallelization with lapply/map patterns.
#' @examples
#'   task_lists <- split_work(total_tasks = 1000, n_workers = 4)
#'   # Worker 1 gets task_lists[[1]], Worker 2 gets task_lists[[2]], etc.
split_work <- function(total_tasks, n_workers) {

  if (total_tasks < 1) stop("total_tasks must be at least 1")
  if (n_workers < 1) stop("n_workers must be at least 1")

  # If more workers than tasks, reduce n_workers
  if (n_workers > total_tasks) {
    warning(sprintf("More workers (%d) than tasks (%d); using %d workers",
                    n_workers, total_tasks, total_tasks))
    n_workers <- total_tasks
  }

  # Distribute tasks evenly
  task_indices <- 1:total_tasks
  task_lists <- split(task_indices, cut(task_indices, breaks = n_workers, labels = FALSE))

  return(task_lists)
}


#' Check Future Package Availability
#'
#' Checks if required parallel processing packages are installed
#'
#' @param also_check Character vector. Additional packages to check (optional)
#' @return Logical. TRUE if all packages available
#' @details
#'   Checks for 'future' and optionally 'future.apply'.
#'   Provides informative error message if packages missing.
#' @examples
#'   if (check_parallel_packages()) {
#'     # Proceed with parallel setup
#'   }
check_parallel_packages <- function(also_check = c("future.apply")) {

  required_pkgs <- c("future", also_check)

  missing_pkgs <- character(0)
  for (pkg in required_pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      missing_pkgs <- c(missing_pkgs, pkg)
    }
  }

  if (length(missing_pkgs) > 0) {
    msg <- sprintf(
      "Required package(s) not installed: %s\nInstall with: install.packages(c(%s))",
      paste(missing_pkgs, collapse = ", "),
      paste(sprintf("\"%s\"", missing_pkgs), collapse = ", ")
    )
    stop(msg)
  }

  return(TRUE)
}


#' Reset to Sequential Processing
#'
#' Convenience function to return to sequential (non-parallel) execution
#'
#' @param verbose Logical. Print message
#' @return NULL (invisibly)
#' @examples
#'   reset_sequential()
reset_sequential <- function(verbose = TRUE) {

  if (!requireNamespace("future", quietly = TRUE)) {
    stop("Package 'future' required")
  }

  future::plan(future::sequential)

  if (verbose) {
    cat("Switched to sequential execution (no parallelization)\n")
  }

  return(invisible(NULL))
}


#' Get Current Parallel Configuration
#'
#' Returns information about current future plan
#'
#' @return List with strategy name and number of workers
#' @examples
#'   config <- get_parallel_config()
#'   cat(sprintf("Current strategy: %s with %d workers\n",
#'               config$strategy, config$workers))
get_parallel_config <- function() {

  if (!requireNamespace("future", quietly = TRUE)) {
    stop("Package 'future' required")
  }

  plan_obj <- future::plan()
  strategy_name <- class(plan_obj)[1]

  # Try to get number of workers
  workers <- tryCatch(
    future::nbrOfWorkers(),
    error = function(e) NA
  )

  return(list(
    strategy = strategy_name,
    workers = workers
  ))
}


#' Parallel Apply with Progress (Optional)
#'
#' Wrapper for future_lapply with optional progress bar
#'
#' @param X Vector or list to iterate over
#' @param FUN Function to apply
#' @param seeds Integer vector of seeds (length must equal length(X))
#' @param show_progress Logical. Show progress bar (requires progressr package)
#' @param ... Additional arguments passed to FUN
#' @return List of results
#' @details
#'   Uses future.apply::future_lapply for parallel execution.
#'   Each iteration gets a unique seed for reproducibility.
#' @examples
#'   seeds <- generate_seeds(100, base_seed = 123)
#'   results <- parallel_apply_with_seeds(
#'     X = 1:100,
#'     FUN = function(i, seed) { set.seed(seed); rnorm(10) },
#'     seeds = seeds
#'   )
parallel_apply_with_seeds <- function(X, FUN, seeds, show_progress = FALSE, ...) {

  if (!requireNamespace("future.apply", quietly = TRUE)) {
    stop("Package 'future.apply' required")
  }

  if (length(seeds) != length(X)) {
    stop("Length of seeds must equal length of X")
  }

  # Wrap FUN to accept seed
  FUN_with_seed <- function(i) {
    seed_i <- seeds[[i]]
    FUN(X[[i]], seed = seed_i, ...)
  }

  # Optional progress
  if (show_progress && requireNamespace("progressr", quietly = TRUE)) {
    progressr::with_progress({
      p <- progressr::progressor(steps = length(X))
      results <- future.apply::future_lapply(seq_along(X), function(i) {
        result <- FUN_with_seed(i)
        p()
        return(result)
      }, future.seed = TRUE)
    })
  } else {
    results <- future.apply::future_lapply(seq_along(X), FUN_with_seed, future.seed = TRUE)
  }

  return(results)
}
