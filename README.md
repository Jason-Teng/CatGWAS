# CatGWAS

`CatGWAS` is an R package for genome-wide association studies (GWAS) with **categorical phenotypes**, including both **nominal** and **ordinal** traits.

The package provides a unified interface for categorical mixed-model GWAS and supports multiple marker-scanning methods, including score tests, pseudo-response methods, P3D, GLM, and exact marker-specific fitting.

---

## Installation

```r
install.packages("remotes")
remotes::install_github("Jason-Teng/CatGWAS")
library(CatGWAS)
```

---

## Quick start

Load the example data:

```r
data(z)
data(kk)
data(nominal)
data(ordinal)

zz <- z[1:5, ]   # use a small number of markers for a fast example
```

### Ordinal GWAS

For ordinal traits, we recommend starting with `score` and `psrsd`.

```r
res_ord <- categorical_gwas(
  y = ordinal,
  zz = zz,
  kk = kk,
  trait_type = "ordinal",
  method = c("score", "psrsd")
)

res_ord$results$score
res_ord$results$psrsd
```

### Nominal GWAS

For nominal traits, we recommend starting with the `score` test.

```r
res_nom <- categorical_gwas(
  y = nominal,
  zz = zz,
  kk = kk,
  trait_type = "nominal",
  method = "score"
)

res_nom$result
```
---

## Main function

All analyses are run through `categorical_gwas()`.

```r
categorical_gwas(
  y,
  zz,
  kk = NULL,
  trait_type = c("ordinal", "nominal"),
  method,
  null_method = "pseudo",
  null_fit = NULL,
  vc = NULL,
  link = "cprobit",
  n_cores = 1L
)
```

| Argument | Description |
|---|---|
| `y` | Categorical phenotype |
| `zz` | Genotype matrix, **markers × individuals** |
| `kk` | Kinship matrix, **individuals × individuals** |
| `trait_type` | `"ordinal"` or `"nominal"` |
| `method` | One or more scanning methods |
| `null_method` | Null-model method for nominal; default `"pseudo"` |
| `null_fit` | Precomputed null model (optional) |
| `vc` | Precomputed variance component (optional) |
| `link` | Ordinal cumulative link (ignored for nominal); default `"cprobit"` |
| `n_cores` | Number of cores for parallel score-test scanning (default `1L`) |

## `method`

| Value | Description | Nominal | Ordinal |
|---|---|---:|---:|
| `score` | Fast score test using the fitted null model | Yes | Yes |
| `p3d` | Population parameters previously determined | Yes | Yes |
| `psr` | Pseudo-response marker scan | Yes | Yes |
| `psrsd` | Pseudo-response with ordinal-specific structure | No | Yes |
| `glm` | GLM without kinship correction | Yes | Yes |
| `exact` | Marker-specific model fitting | Yes | Yes |

`score`, `p3d`, `psr`, and `psrsd` need a null model (`kk` required). `glm` and `exact` do not. If `null_fit` / `vc` are not supplied, the null is fitted automatically.

## `null_method` (nominal only)

Default is `"pseudo"`. Ordinal traits always use pseudo and ignore this argument.

| Value | Description |
|---|---|
| `"pseudo"` | Pseudo-response variance-component estimation (default) |
| `"laplace"` | Laplace approximation |

The pseudo-response approach is faster. The Laplace approximation is usually less biased, but takes more time. In practice, this extra cost is often acceptable because the null model only needs to be fitted once. The fitted null model or estimated variance component can then be reused for genome-wide marker scanning.

A nominal GWAS using the Laplace null model can be run as:

```r
res_nom_laplace <- categorical_gwas(
  y = nominal,
  zz = zz,
  kk = kk,
  trait_type = "nominal",
  method = "score",
  null_method = "laplace"
)

res_nom_laplace$result
```
## `null_fit`
If a fitted null model is already available, it can be reused:

```r
res_nom_reuse <- categorical_gwas(
  y = nominal,
  zz = zz,
  kk = kk,
  trait_type = "nominal",
  method = "score",
  null_fit = res_nom_laplace$null_fit
)
```

## `vc`
If a variance component has already been estimated, it can also be supplied directly:

```r
res_nom_vc <- categorical_gwas(
  y = nominal,
  zz = zz,
  kk = kk,
  trait_type = "nominal",
  method = "score",
  vc = res_nom_laplace$null_fit$par
)
```


## `link` (ordinal only)
For ordinal traits, `CatGWAS` supports two cumulative links:

| Value | Model |
|---|---|
| `"cprobit"` | Cumulative probit (default) |
| `"clogit"` | Cumulative logit (proportional odds) |

The default is cumulative probit. The cumulative logit link can be used with:

```r
res_ord_clogit <- categorical_gwas(
  y = ordinal,
  zz = zz,
  kk = kk,
  trait_type = "ordinal",
  method = c("score", "psr"),
  link = "clogit"
)
```
The link function is used in the ordinal null model and inherited by methods that depend on the fitted null model.


## `n_cores`

Used by `method = "score"` only. `1L` is serial. `n_cores > 1` uses PSOCK workers. The null model is always fitted serially; only per-marker score calculations are parallelised.

```r
res_ord_parallel <- categorical_gwas(
  y = ordinal,
  zz = zz,
  kk = kk,
  trait_type = "ordinal",
  method = "score",
  n_cores = 4L
)
```
---

## Input

| Object | Format |
|---|---|
| `zz` | Genotype matrix, **markers × individuals** |
| `kk` | Kinship matrix, **individuals × individuals** |
| `nominal` | Factor |
| `ordinal` | Factor |

---

## Output

`categorical_gwas()` returns a list. Marker results are in `result` (one method) or `results` (several methods).

| Column | Meaning |
|---|---|
| `SNP` | Marker index |
| `Effect` | Estimated marker effect |
| `StdErr` | Standard error |
| `Score` / `Wald` | Test statistic |
| `p` | P-value |
| `iter`, `err` | Convergence, when applicable |

The fitted null, when used, is returned as `null_fit`. For ordinal scans, `link` is also returned.

---

## Citation

```text
Teng, C.-S. CatGWAS: Genome-wide association studies for categorical phenotypes.
GitHub repository: https://github.com/Jason-Teng/CatGWAS
```

---

## Author

Chin-Sheng Teng  
University of California, Riverside
