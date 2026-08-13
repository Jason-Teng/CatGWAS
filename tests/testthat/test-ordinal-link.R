# tests/testthat/test-ordinal-link.R
#
# Tests for the ordinal cumulative link function.
# Default remains cumulative probit (pnorm / dnorm).
# link = "clogit" switches to cumulative logit (plogis / dlogis).

make_synthetic <- function(seed = 42L, n = 20L, m = 4L, c = 3L) {
  set.seed(seed)
  kk <- crossprod(matrix(rnorm(n * n), n, n))
  kk <- kk / mean(diag(kk))

  zz <- matrix(sample(0:2, m * n, replace = TRUE), nrow = m, ncol = n)
  rownames(zz) <- paste0("SNP", seq_len(m))

  y_int <- sample(seq_len(c), n, replace = TRUE)
  y     <- factor(y_int, levels = seq_len(c))

  list(kk = kk, zz = zz, y = y, n = n, m = m, c = c)
}

# =======================================================================
# 1. ordinal_link() itself
# =======================================================================
test_that("ordinal_link returns CDF/PDF for probit and clogit", {
  probit <- CatGWAS:::ordinal_link("probit")
  clogit <- CatGWAS:::ordinal_link("clogit")
  logit  <- CatGWAS:::ordinal_link("logit")

  expect_equal(probit$name, "probit")
  expect_equal(clogit$name, "clogit")
  expect_equal(logit$name,  "clogit")

  x <- c(-2, -0.5, 0, 0.5, 2)
  expect_equal(probit$p(x), pnorm(x))
  expect_equal(probit$d(x), dnorm(x))
  expect_equal(clogit$p(x), plogis(x))
  expect_equal(clogit$d(x), dlogis(x))
})

test_that("ordinal_link rejects invalid values", {
  expect_error(CatGWAS:::ordinal_link("identity"), regexp = "probit")
  expect_error(CatGWAS:::ordinal_link(c("probit", "clogit")), regexp = "single")
  expect_error(CatGWAS:::ordinal_link(1), regexp = "character")
})

# =======================================================================
# 2. Default remains probit; clogit is a distinct model
# =======================================================================
test_that("default ordinal glm matches explicit probit", {
  d <- make_synthetic()

  res_default <- categorical_gwas(
    y = d$y, zz = d$zz, kk = d$kk,
    trait_type = "ordinal", method = "glm"
  )
  res_probit <- categorical_gwas(
    y = d$y, zz = d$zz, kk = d$kk,
    trait_type = "ordinal", method = "glm",
    link = "probit"
  )

  expect_equal(res_default$link, "probit")
  expect_equal(res_probit$link,  "probit")
  expect_equal(res_default$result$Wald, res_probit$result$Wald, tolerance = 1e-10)
  expect_equal(res_default$result$p,    res_probit$result$p,    tolerance = 1e-10)
})

test_that("clogit glm runs and differs from probit", {
  d <- make_synthetic()

  res_probit <- categorical_gwas(
    y = d$y, zz = d$zz, kk = d$kk,
    trait_type = "ordinal", method = "glm",
    link = "probit"
  )
  res_clogit <- categorical_gwas(
    y = d$y, zz = d$zz, kk = d$kk,
    trait_type = "ordinal", method = "glm",
    link = "clogit"
  )

  expect_equal(res_clogit$link, "clogit")
  expect_equal(nrow(res_clogit$result), d$m)
  expect_true(all(res_clogit$result$p >= 0 & res_clogit$result$p <= 1))
  expect_false(isTRUE(all.equal(res_probit$result$p, res_clogit$result$p)))
})

test_that("clogit score test runs and stores the link on the null fit", {
  d <- make_synthetic()

  res <- categorical_gwas(
    y = d$y, zz = d$zz, kk = d$kk,
    trait_type = "ordinal", method = "score",
    link = "clogit"
  )

  expect_equal(res$link, "clogit")
  expect_equal(res$null_fit$link, "clogit")
  expect_equal(nrow(res$result), d$m)
  expect_true(all(res$result$p >= 0 & res$result$p <= 1))
})

test_that("logit is accepted as an alias of clogit", {
  d <- make_synthetic(m = 2L)

  res <- categorical_gwas(
    y = d$y, zz = d$zz, kk = d$kk,
    trait_type = "ordinal", method = "glm",
    link = "logit"
  )

  expect_equal(res$link, "clogit")
})

test_that("invalid ordinal link is rejected", {
  d <- make_synthetic(m = 2L)
  expect_error(
    categorical_gwas(
      y = d$y, zz = d$zz, kk = d$kk,
      trait_type = "ordinal", method = "glm",
      link = "identity"
    ),
    regexp = "probit"
  )
})
