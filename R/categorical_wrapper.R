#' Run categorical GWAS
#'
#' @param y Phenotype object.
#' @param zz Marker genotype matrix, markers in rows and individuals in columns.
#' @param kk Kinship matrix.
#' @param x0 Optional fixed-effect design matrix.
#' @param trait_type Either "nominal" or "ordinal".
#' @param method Character vector of methods.
#' @param null_method Null-model fitting method.
#' @param null_fit Precomputed null model object (optional).
#' @param vc Precomputed variance components (optional).
#' @param ids Individual identifiers (optional).
#' @param outdir Directory for optional output files.
#' @param maxiter Maximum iterations for iterative fitting.
#' @param minerr Convergence tolerance.
#' @param n_cores Number of cores for parallel score-test scanning.
#'   Must be a single positive integer. Use \code{1L} (the default) for
#'   serial computation. Parallelisation uses base-R PSOCK clusters and
#'   therefore works on Windows, macOS, and Linux. The null model is always
#'   fitted serially; only the per-marker score calculations are parallelised.
#' @param link Ordinal cumulative link. Used only when
#'   \code{trait_type = "ordinal"}. \code{"cprobit"} (the default) is the
#'   cumulative probit model; \code{"clogit"} is the cumulative logit /
#'   proportional-odds model. \code{"probit"} is accepted as an alias of
#'   \code{"cprobit"}. \code{"logit"} is accepted as an alias of
#'   \code{"clogit"}. Ignored for nominal traits.
#' @return A list of GWAS results.
#' @export
categorical_gwas <- function(y,
                             zz,
                             kk = NULL,
                             x0 = NULL,
                             trait_type = c("ordinal", "nominal"),
                             method = NULL,
                             null_method = NULL,
                             null_fit = NULL,
                             vc = NULL,
                             ids = NULL,
                             outdir = NULL,
                             maxiter = 100,
                             minerr = 1e-8,
                             n_cores = 1L,
                             link = "cprobit") {

  trait_type <- match.arg(trait_type)
  n_cores <- .validate_n_cores(n_cores)

  if (trait_type == "nominal" && !identical(link, "cprobit")) {
    warning("'link' is only used for ordinal traits and is ignored for nominal traits.",
            call. = FALSE)
  }

  if (trait_type == "ordinal") {
    if (is.null(method)) {
      method <- c("score", "p3d", "psr", "psrsd", "exact", "glm")
    }
    if (is.null(null_method)) {
      null_method <- "pseudo"
    }

    return(ordinal_gwas(
      y = y,
      zz = zz,
      kk = kk,
      x0 = x0,
      method = method,
      null_method = null_method,
      null_fit = null_fit,
      vc = vc,
      ids = ids,
      outdir = outdir,
      maxiter = maxiter,
      minerr = minerr,
      n_cores = n_cores,
      link = link
    ))
  }

  if (trait_type == "nominal") {
    if (is.null(method)) {
      method <- c("score", "p3d", "psr", "exact", "glm")
    }
    if (is.null(null_method)) {
      null_method <- "pseudo"
    }

    # The current nominal_gwas() accepts only one method at a time.
    # This dispatcher allows multiple nominal methods by looping.
    allowed_nominal <- c("score", "p3d", "psr", "exact", "glm")
    method <- match.arg(method, choices = allowed_nominal, several.ok = TRUE)

    nominal_results <- list()
    nominal_times <- list()
    nominal_error_snps <- list()
    shared_null_fit <- null_fit

    for (one_method in method) {
      fit <- nominal_gwas(
        y = y,
        zz = zz,
        kk = kk,
        x0 = x0,
        method = one_method,
        null_method = null_method,
        null_fit = shared_null_fit,
        vc = vc,
        ids = ids,
        outdir = outdir,
        maxiter = maxiter,
        minerr = minerr,
        n_cores = n_cores
      )

      nominal_results[[one_method]] <- fit$result
      nominal_times[[one_method]] <- fit$gwas_scanning_time
      nominal_error_snps[[one_method]] <- fit$error_snps

      if (is.null(shared_null_fit) && !is.null(fit$null_fit)) {
        shared_null_fit <- fit$null_fit
      }
    }

    if (length(method) == 1) {
      return(list(
        trait_type = "nominal",
        method = method,
        result = nominal_results[[method]],
        error_snps = nominal_error_snps[[method]],
        null_fit = if (method %in% c("score", "p3d", "psr")) shared_null_fit else NULL,
        gwas_scanning_time = nominal_times[[method]]
      ))
    }

    return(list(
      trait_type = "nominal",
      method = method,
      results = nominal_results,
      error_snps = nominal_error_snps,
      null_fit = if (any(method %in% c("score", "p3d", "psr"))) shared_null_fit else NULL,
      gwas_scanning_time = nominal_times
    ))
  }
}
