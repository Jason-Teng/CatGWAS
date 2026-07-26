# tests/testthat/test-parallel-score.R
#
# Tests for the n_cores parallel score-test implementation.
# Unit tests use a small deterministic synthetic dataset.
# Integration tests use the IMF2 example data shipped with the package.

# -----------------------------------------------------------------------
# Helper: build a tiny but deterministic dataset
# n = 20 individuals, m = 6 markers, c = 3 categories
# -----------------------------------------------------------------------
make_synthetic <- function(seed = 42L, n = 20L, m = 6L, c = 3L) {
  set.seed(seed)
  kk <- crossprod(matrix(rnorm(n * n), n, n))
  kk <- kk / mean(diag(kk))

  zz <- matrix(sample(0:2, m * n, replace = TRUE), nrow = m, ncol = n)
  rownames(zz) <- paste0("SNP", seq_len(m))

  y_int <- sample(seq_len(c), n, replace = TRUE)
  y     <- factor(y_int, levels = seq_len(c))

  list(kk = kk, zz = zz, y = y, n = n, m = m, c = c)
}

# -----------------------------------------------------------------------
# Helper: load IMF2 example data from extdata CSV files
# -----------------------------------------------------------------------
load_imf2 <- function() {
  ext <- system.file("extdata", package = "categoricalGWAS")
  if (!nzchar(ext)) return(NULL)

  gen_path <- file.path(ext, "IMF2-Genotypes.csv")
  phe_path <- file.path(ext, "IMF2-Phenotypes.csv")

  if (!file.exists(gen_path) || !file.exists(phe_path)) return(NULL)

  gen <- read.csv(gen_path)
  phe <- read.csv(phe_path)

  zz_full <- as.matrix(gen[, -c(1:5)])
  # Compute kinship from genotypes (Van Raden style, centred + scaled)
  Z <- scale(t(zz_full), center = TRUE, scale = FALSE)
  kk <- tcrossprod(Z) / ncol(Z)

  nominal <- as.factor(phe$Nominal)
  ordinal <- as.factor(phe$Ordinal)

  list(zz = zz_full, kk = kk, nominal = nominal, ordinal = ordinal)
}

# =======================================================================
# 1. Validation tests
# =======================================================================
test_that(".validate_n_cores rejects invalid values", {
  expect_error(categoricalGWAS:::.validate_n_cores(0L),   regexp = "positive integer")
  expect_error(categoricalGWAS:::.validate_n_cores(-1L),  regexp = "positive integer")
  expect_error(categoricalGWAS:::.validate_n_cores(NA),   regexp = "positive integer")
  expect_error(categoricalGWAS:::.validate_n_cores("two"), regexp = "positive integer")
  expect_error(categoricalGWAS:::.validate_n_cores(c(1L, 2L)), regexp = "single positive integer")
})

test_that(".validate_n_cores accepts valid values", {
  expect_equal(categoricalGWAS:::.validate_n_cores(1L),  1L)
  expect_equal(categoricalGWAS:::.validate_n_cores(2L),  2L)
  expect_equal(categoricalGWAS:::.validate_n_cores(1.0), 1L)  # coerced
  expect_equal(categoricalGWAS:::.validate_n_cores(4),   4L)
})

test_that("categorical_gwas errors on invalid n_cores", {
  d <- make_synthetic()
  expect_error(
    categorical_gwas(y = d$y, zz = d$zz, kk = d$kk,
                     trait_type = "ordinal", method = "score",
                     n_cores = 0L),
    regexp = "positive integer"
  )
})

# =======================================================================
# 2. Ordinal score test: serial == parallel (synthetic data)
# =======================================================================
test_that("ordinal score: serial and parallel give identical results", {
  d <- make_synthetic()

  res1 <- categorical_gwas(y = d$y, zz = d$zz, kk = d$kk,
                           trait_type = "ordinal", method = "score",
                           n_cores = 1L)
  res2 <- categorical_gwas(y = d$y, zz = d$zz, kk = d$kk,
                           trait_type = "ordinal", method = "score",
                           n_cores = 2L)

  tab1 <- res1$result
  tab2 <- res2$result

  # Same dimensions and column names
  expect_equal(dim(tab1), dim(tab2))
  expect_equal(colnames(tab1), colnames(tab2))

  # Same marker order
  expect_equal(tab1$SNP, tab2$SNP)
  expect_equal(tab1$SNP, seq_len(d$m))

  # No duplicates
  expect_equal(length(unique(tab1$SNP)), d$m)

  # Numerically equal statistics and p-values
  expect_equal(tab1$Score,  tab2$Score,  tolerance = 1e-10)
  expect_equal(tab1$p,      tab2$p,      tolerance = 1e-10)
  expect_equal(tab1$Effect, tab2$Effect, tolerance = 1e-10)
  expect_equal(tab1$StdErr, tab2$StdErr, tolerance = 1e-10)
})

# =======================================================================
# 3. Nominal score test: serial == parallel (synthetic data)
# =======================================================================
test_that("nominal score: serial and parallel give identical results", {
  d <- make_synthetic()

  res1 <- categorical_gwas(y = d$y, zz = d$zz, kk = d$kk,
                           trait_type = "nominal", method = "score",
                           n_cores = 1L)
  res2 <- categorical_gwas(y = d$y, zz = d$zz, kk = d$kk,
                           trait_type = "nominal", method = "score",
                           n_cores = 2L)

  tab1 <- res1$result
  tab2 <- res2$result

  expect_equal(dim(tab1), dim(tab2))
  expect_equal(colnames(tab1), colnames(tab2))
  expect_equal(tab1$SNP, tab2$SNP)
  expect_equal(tab1$SNP, seq_len(d$m))
  expect_equal(length(unique(tab1$SNP)), d$m)
  expect_equal(tab1$Score, tab2$Score, tolerance = 1e-10)
  expect_equal(tab1$p,     tab2$p,     tolerance = 1e-10)
})

# =======================================================================
# 4. Non-score methods are not affected by n_cores (ordinal)
# =======================================================================
test_that("ordinal glm gives same result regardless of n_cores", {
  d <- make_synthetic(m = 3L)

  res1 <- categorical_gwas(y = d$y, zz = d$zz, kk = d$kk,
                           trait_type = "ordinal", method = "glm",
                           n_cores = 1L)
  res2 <- categorical_gwas(y = d$y, zz = d$zz, kk = d$kk,
                           trait_type = "ordinal", method = "glm",
                           n_cores = 2L)

  expect_equal(res1$result$Wald, res2$result$Wald, tolerance = 1e-10)
  expect_equal(res1$result$p,    res2$result$p,    tolerance = 1e-10)
})

# =======================================================================
# 5. Integration test: IMF2 example data (full scan, ordinal + nominal)
# =======================================================================
test_that("integration: ordinal score serial == parallel on IMF2 data", {
  d <- load_imf2()
  skip_if(is.null(d), "IMF2 example data not available")

  # Subset to first 50 markers for reasonable test speed
  m_sub <- 50L
  zz_sub <- d$zz[seq_len(m_sub), , drop = FALSE]

  t1 <- system.time({
    res1 <- categorical_gwas(y = d$ordinal, zz = zz_sub, kk = d$kk,
                             trait_type = "ordinal", method = "score",
                             n_cores = 1L)
  })

  t2 <- system.time({
    res2 <- categorical_gwas(y = d$ordinal, zz = zz_sub, kk = d$kk,
                             trait_type = "ordinal", method = "score",
                             null_fit  = res1$null_fit,   # reuse null model
                             n_cores = 2L)
  })

  tab1 <- res1$result
  tab2 <- res2$result

  max_diff_score <- max(abs(tab1$Score - tab2$Score))
  max_diff_p     <- max(abs(tab1$p     - tab2$p    ))

  message(sprintf(
    "[IMF2 ordinal] serial=%.2fs  parallel=%.2fs  max|score diff|=%.2e  max|p diff|=%.2e",
    t1["elapsed"], t2["elapsed"], max_diff_score, max_diff_p
  ))

  expect_equal(tab1$SNP,   tab2$SNP)
  expect_equal(nrow(tab1), m_sub)
  expect_equal(length(unique(tab2$SNP)), m_sub)
  expect_equal(tab1$Score, tab2$Score, tolerance = 1e-10)
  expect_equal(tab1$p,     tab2$p,     tolerance = 1e-10)
})

test_that("integration: nominal score serial == parallel on IMF2 data", {
  d <- load_imf2()
  skip_if(is.null(d), "IMF2 example data not available")

  m_sub <- 50L
  zz_sub <- d$zz[seq_len(m_sub), , drop = FALSE]

  t1 <- system.time({
    res1 <- categorical_gwas(y = d$nominal, zz = zz_sub, kk = d$kk,
                             trait_type = "nominal", method = "score",
                             n_cores = 1L)
  })

  t2 <- system.time({
    res2 <- categorical_gwas(y = d$nominal, zz = zz_sub, kk = d$kk,
                             trait_type = "nominal", method = "score",
                             null_fit  = res1$null_fit,
                             n_cores = 2L)
  })

  tab1 <- res1$result
  tab2 <- res2$result

  max_diff_score <- max(abs(tab1$Score - tab2$Score))
  max_diff_p     <- max(abs(tab1$p     - tab2$p    ))

  message(sprintf(
    "[IMF2 nominal] serial=%.2fs  parallel=%.2fs  max|score diff|=%.2e  max|p diff|=%.2e",
    t1["elapsed"], t2["elapsed"], max_diff_score, max_diff_p
  ))

  expect_equal(tab1$SNP,   tab2$SNP)
  expect_equal(nrow(tab1), m_sub)
  expect_equal(length(unique(tab2$SNP)), m_sub)
  expect_equal(tab1$Score, tab2$Score, tolerance = 1e-10)
  expect_equal(tab1$p,     tab2$p,     tolerance = 1e-10)
})

# =======================================================================
# 6. Full-data integration test (all markers) reporting runtime
# =======================================================================
test_that("integration: full IMF2 ordinal scan runtime and correctness report", {
  d <- load_imf2()
  skip_if(is.null(d), "IMF2 example data not available")

  m_full <- nrow(d$zz)

  t1 <- system.time({
    res1 <- categorical_gwas(y = d$ordinal, zz = d$zz, kk = d$kk,
                             trait_type = "ordinal", method = "score",
                             n_cores = 1L)
  })

  t2 <- system.time({
    res2 <- categorical_gwas(y = d$ordinal, zz = d$zz, kk = d$kk,
                             trait_type = "ordinal", method = "score",
                             null_fit = res1$null_fit,
                             n_cores = 2L)
  })

  tab1 <- res1$result
  tab2 <- res2$result

  n_failed <- sum(is.na(tab1$Score) | is.na(tab2$Score))
  max_diff  <- max(abs(tab1$Score - tab2$Score), na.rm = TRUE)

  message(sprintf(
    "[FULL IMF2 ordinal] markers=%d  serial=%.1fs  parallel=%.1fs  max|diff|=%.2e  failed=%d",
    m_full, t1["elapsed"], t2["elapsed"], max_diff, n_failed
  ))

  expect_equal(tab1$SNP,   tab2$SNP)
  expect_equal(nrow(tab2), m_full)
  expect_equal(tab1$Score, tab2$Score, tolerance = 1e-10)
  expect_equal(n_failed,   0L)
})

test_that("integration: full IMF2 nominal scan runtime and correctness report", {
  d <- load_imf2()
  skip_if(is.null(d), "IMF2 example data not available")

  m_full <- nrow(d$zz)

  t1 <- system.time({
    res1 <- categorical_gwas(y = d$nominal, zz = d$zz, kk = d$kk,
                             trait_type = "nominal", method = "score",
                             n_cores = 1L)
  })

  t2 <- system.time({
    res2 <- categorical_gwas(y = d$nominal, zz = d$zz, kk = d$kk,
                             trait_type = "nominal", method = "score",
                             null_fit = res1$null_fit,
                             n_cores = 2L)
  })

  tab1 <- res1$result
  tab2 <- res2$result

  n_failed <- sum(is.na(tab1$Score) | is.na(tab2$Score))
  max_diff  <- max(abs(tab1$Score - tab2$Score), na.rm = TRUE)

  message(sprintf(
    "[FULL IMF2 nominal] markers=%d  serial=%.1fs  parallel=%.1fs  max|diff|=%.2e  failed=%d",
    m_full, t1["elapsed"], t2["elapsed"], max_diff, n_failed
  ))

  expect_equal(tab1$SNP,   tab2$SNP)
  expect_equal(nrow(tab2), m_full)
  expect_equal(tab1$Score, tab2$Score, tolerance = 1e-10)
  expect_equal(n_failed,   0L)
})
