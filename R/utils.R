# Internal utilities for categoricalGWAS

#' Validate the n_cores argument
#'
#' @param n_cores Value supplied by the caller.
#' @return A length-1 integer >= 1L.
#' @keywords internal
.validate_n_cores <- function(n_cores) {
  if (length(n_cores) != 1L) {
    stop("'n_cores' must be a single positive integer, not a vector of length ",
         length(n_cores), ".", call. = FALSE)
  }
  n_cores <- suppressWarnings(as.integer(n_cores))
  if (is.na(n_cores) || n_cores < 1L) {
    stop("'n_cores' must be a positive integer >= 1.", call. = FALSE)
  }
  n_cores
}
