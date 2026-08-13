# CatGWAS

GWAS for **nominal** and **ordinal** categorical phenotypes.

---

## Installation

```r
install.packages("remotes")
remotes::install_github("Jason-Teng/CatGWAS")
library(CatGWAS)
```

---

## Examples

```r
data(z)
data(kk)
data(nominal)
data(ordinal)
zz <- z[1:5, ]
```

Ordinal
```r
# Ordinal, default cprobit
res_ord <- categorical_gwas(
  y = ordinal, zz = zz, kk = kk,
  trait_type = "ordinal",
  method = c("score", "psrsd")
)

# Ordinal, cumulative logit, parallel score
res_ord_clogit <- categorical_gwas(
  y = ordinal, zz = zz, kk = kk,
  trait_type = "ordinal",
  method = "score",
  link = "clogit",
  n_cores = 4L
)

# Reuse a fitted null
res_reuse <- categorical_gwas(
  y = ordinal, zz = zz, kk = kk,
  trait_type = "ordinal",
  method = c("score", "psr"),
  link = "clogit",
  null_fit = res_ord_clogit$null_fit
)
```

Nominal
```r
# Default pseudo null
res_nom <- categorical_gwas(
  y = nominal, zz = zz, kk = kk,
  trait_type = "nominal",
  method = c("score", "psr")
)
# Laplace null
res_nom_laplace <- categorical_gwas(
  y = nominal, zz = zz, kk = kk,
  trait_type = "nominal",
  method = c("score", "psr"),
  null_method = "laplace"
)
# Reuse a precomputed variance component
res_nom_vc <- categorical_gwas(
  y = nominal, zz = zz, kk = kk,
  trait_type = "nominal",
  method = c("score", "psr"),
  vc = vc
)
```

```r
res_ord$results$score
res_nom$results$score
```

---

## Functions

### `categorical_gwas()`

Main interface. All analyses go through this function.

| Argument | Description |
|---|---|
| `y` | Phenotype (factor or category-indicator matrix) |
| `zz` | Genotype matrix, **markers × individuals** |
| `kk` | Kinship matrix, **individuals × individuals** |
| `trait_type` | `"ordinal"` or `"nominal"` |
| `method` | One or more scanning methods (see below) |
| `null_method` | Null-model method for nominal; default `"pseudo"` |
| `null_fit` | Precomputed null model (optional) |
| `vc` | Precomputed variance component (optional) |
| `link` | Ordinal cumulative link (ignored for nominal); default `"cprobit"` |
| `n_cores` | Cores for parallel **score** scanning (default `1L`) |

## Arguments
### `method`

| Value | Description | Nominal | Ordinal |
|---|---|---:|---:|
| `score` | Score test using the fitted null model | Yes | Yes |
| `p3d` | Population parameters previously determined | Yes | Yes |
| `psr` | Pseudo-response marker scan | Yes | Yes |
| `psrsd` | Pseudo-response with ordinal-specific structure | No | Yes |
| `glm` | GLM without kinship correction | Yes | Yes |
| `exact` | Marker-specific model fitting | Yes | Yes |

`score`, `p3d`, `psr`, and `psrsd` need a null model (`kk` required). `glm` and `exact` do not. If `null_fit` / `vc` are not supplied, the null is fitted automatically.

### `null_method` (nominal only)

Default is `"pseudo"`. Ordinal traits always use pseudo and ignore this argument.

| Value | Description |
|---|---|
| `"pseudo"` | Pseudo-response variance-component estimation (default) |
| `"laplace"` | Laplace approximation |

### `link` (ordinal only)

| Value | Model |
|---|---|
| `"cprobit"` | Cumulative probit (default) |
| `"clogit"` | Cumulative logit (proportional odds) |

Used in the ordinal null model and in P3D, GLM, and exact. Score, PSR, and PSRSD inherit it from the null.

### `n_cores`

Used by `method = "score"` only. `1L` is serial. `n_cores > 1` uses PSOCK workers. The null model is always fitted serially; only per-marker score calculations are parallelised.

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
