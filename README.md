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

## Functions

### `categorical_gwas()`

Main interface. All analyses go through this function.

| Argument | Description |
|---|---|
| `y` | Phenotype (factor or category-indicator matrix) |
| `zz` | Genotype matrix, **markers × individuals** |
| `kk` | Kinship matrix, **individuals × individuals** |
| `x0` | Optional fixed-effect design matrix |
| `trait_type` | `"ordinal"` or `"nominal"` |
| `method` | One or more scanning methods (see below) |
| `null_method` | Null-model variance-component method (see below) |
| `null_fit` | Precomputed null model (optional) |
| `vc` | Precomputed variance component (optional) |
| `link` | Ordinal cumulative link (ignored for nominal) |
| `n_cores` | Cores for parallel **score** scanning (default `1L`) |
| `maxiter`, `minerr` | Iteration control |
| `ids`, `outdir` | Optional IDs and output directory |

#### `method`

| Value | Description | Nominal | Ordinal |
|---|---|---:|---:|
| `score` | Score test using the fitted null model | Yes | Yes |
| `p3d` | Population parameters previously determined | Yes | Yes |
| `psr` | Pseudo-response marker scan | Yes | Yes |
| `psrsd` | Pseudo-response with ordinal-specific structure | No | Yes |
| `glm` | GLM without kinship correction | Yes | Yes |
| `exact` | Marker-specific model fitting | Yes | Yes |

`score`, `p3d`, `psr`, and `psrsd` need a null model (`kk` required). `glm` and `exact` do not. If `null_fit` / `vc` are not supplied, the null is fitted automatically.

#### `null_method`

| Value | Trait | Description |
|---|---|---|
| `pseudo` | Both | Pseudo-response variance-component estimation (default) |
| `laplace` | Nominal only | Laplace approximation |

#### `link` (ordinal only)

| Value | Model |
|---|---|
| `"cprobit"` | Cumulative probit (default) |
| `"clogit"` | Cumulative logit (proportional odds) |
| `"probit"` | Alias of `"cprobit"` |
| `"logit"` | Alias of `"clogit"` |

Used in the ordinal null model and in P3D, GLM, and exact. Score, PSR, and PSRSD inherit it from the null (`ps`, `rr`).

#### `n_cores`

Used by `method = "score"` only. `1L` is serial. `n_cores > 1` uses PSOCK workers (Windows / macOS / Linux). The null model is always fitted serially; only per-marker score calculations are parallelised. Works with either ordinal link.

---

## Examples

```r
data(z)
data(kk)
data(nominal)
data(ordinal)
zz <- z[1:5, ]

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

# Nominal
res_nom <- categorical_gwas(
  y = nominal, zz = zz, kk = kk,
  trait_type = "nominal",
  method = c("score", "psr")
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

## Input

| Object | Format |
|---|---|
| `z` / `zz` | Genotype matrix, **markers × individuals** |
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
