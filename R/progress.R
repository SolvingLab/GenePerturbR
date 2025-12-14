#' Progress Bar System for GenePerturbR
#'
#' Unified progress bar for batch operations
#'
#' @keywords internal

#' Create Progress Bar
#'
#' @param n Integer. Total number of iterations
#' @param title Character. Progress bar title
#' @param show_progress Logical. Whether to show progress
#'
#' @return Progress bar object or NULL
#' @keywords internal
.gpdb_progress_bar <- function(n, title = "Processing", show_progress = TRUE) {
  if (!show_progress || n <= 1) {
    return(NULL)
  }

  # Use simple text progress bar (built-in, no dependencies)
  utils::txtProgressBar(
    min = 0,
    max = n,
    style = 3,
    width = 50
  )
}

#' Update Progress Bar
#'
#' @param pb Progress bar object
#' @param value Current value
#'
#' @keywords internal
.gpdb_update_progress <- function(pb, value) {
  if (!is.null(pb)) {
    utils::setTxtProgressBar(pb, value)
  }
}

#' Close Progress Bar
#'
#' @param pb Progress bar object
#'
#' @keywords internal
.gpdb_close_progress <- function(pb) {
  if (!is.null(pb)) {
    close(pb)
  }
}
