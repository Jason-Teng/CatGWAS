# ============================================================
# Ordinal link function
# ============================================================
# Cumulative ordinal models differ only in the latent error
# distribution used to form category probabilities and the
# derivative matrix D:
#
#   cprobit -> standard normal   (pnorm / dnorm)
#   clogit  -> standard logistic (plogis / dlogis)
#
# Everything else in the IBLUP / pseudo-response loop is unchanged.
# ============================================================

#' Ordinal cumulative link function
#'
#' Returns the CDF and PDF of the latent error distribution used by
#' the cumulative ordinal model.
#'
#' @param link One of \code{"cprobit"} (cumulative probit, the default)
#'   or \code{"clogit"} (cumulative logit / proportional odds).
#'   \code{"probit"} is accepted as an alias of \code{"cprobit"}.
#'   \code{"logit"} is accepted as an alias of \code{"clogit"}.
#' @return A list with components \code{name} (character), \code{p}
#'   (CDF), and \code{d} (PDF).
#' @keywords internal
ordinal_link <- function(link = "cprobit") {
  if (length(link) != 1L) {
    stop("'link' must be a single character string.", call. = FALSE)
  }
  if (!is.character(link) || is.na(link)) {
    stop("'link' must be a character string: \"cprobit\" or \"clogit\".",
         call. = FALSE)
  }

  link <- tolower(link)
  if (identical(link, "probit")) {
    link <- "cprobit"
  }
  if (identical(link, "logit")) {
    link <- "clogit"
  }

  link <- match.arg(link, choices = c("cprobit", "clogit"))

  if (identical(link, "cprobit")) {
    list(name = "cprobit", p = stats::pnorm, d = stats::dnorm)
  } else {
    list(name = "clogit", p = stats::plogis, d = stats::dlogis)
  }
}
