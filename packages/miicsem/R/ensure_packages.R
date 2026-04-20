#' Ensure Required Runtime Packages Are Installed
#'
#' Checks whether each of the runtime dependencies (\code{lavaan},
#' \code{mice}, \code{MASS}, \code{SeedMaker}, \code{pbapply}) is
#' available in the current R library path. Installs any that are
#' missing via \code{install.packages()}.
#'
#' Intended for HPC use: after \code{devtools::install_github(...)} of
#' miicsem itself, a single call to \code{run_simulation()} will ensure
#' all downstream packages are present before the main loop starts.
#'
#' @param install Logical. If TRUE (default), missing packages are
#'   installed via \code{install.packages()}. If FALSE, the function
#'   errors out listing what is missing.
#' @param repos Character. CRAN mirror used by \code{install.packages()}
#'   when a repo isn't already set by the R session. Default
#'   "https://cloud.r-project.org".
#' @param verbose Logical. Emit progress messages. Default TRUE.
#' @return Invisible character vector of packages that were installed
#'   (empty if none needed).
#' @export
ensure_required_packages <- function(install = TRUE,
                                     repos  = "https://cloud.r-project.org",
                                     verbose = TRUE) {
  required <- c("lavaan", "mice", "MASS", "SeedMaker", "pbapply")

  missing_pkgs <- required[!vapply(required, requireNamespace,
                                   logical(1), quietly = TRUE)]

  if (length(missing_pkgs) == 0) {
    if (verbose) message("All required packages are available.")
    return(invisible(character(0)))
  }

  if (!install) {
    stop("Required packages not installed: ",
         paste(missing_pkgs, collapse = ", "),
         ". Call ensure_required_packages(install = TRUE) or run ",
         "install.packages() manually.")
  }

  # Use session repos if already set; otherwise fall back to the argument
  current_repos <- getOption("repos")
  repos_to_use <- if (is.null(current_repos) ||
                      any(current_repos == "@CRAN@") ||
                      any(current_repos == "")) {
    repos
  } else {
    current_repos
  }

  if (verbose) {
    message(sprintf("Installing missing packages (%s) from %s ...",
                    paste(missing_pkgs, collapse = ", "),
                    if (is.character(repos_to_use)) repos_to_use[1] else "session repos"))
  }

  utils::install.packages(missing_pkgs, repos = repos_to_use)

  still_missing <- missing_pkgs[!vapply(missing_pkgs, requireNamespace,
                                        logical(1), quietly = TRUE)]
  if (length(still_missing) > 0) {
    stop("Failed to install: ", paste(still_missing, collapse = ", "),
         ". Install manually and retry.")
  }

  if (verbose) {
    message(sprintf("Successfully installed: %s",
                    paste(missing_pkgs, collapse = ", ")))
  }

  invisible(missing_pkgs)
}
