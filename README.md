# saemMicrobiome

`saemMicrobiome` implements two zero-inflated mixed models for longitudinal
microbiome data, estimated with the Stochastic Approximation EM (SAEM)
algorithm:

- **ZIBR** (zero-inflated beta regression) for proportions or relative
  abundances — see [`fit_zibr()`].
- **ZIBBMR** (zero-inflated beta-binomial mixed regression) for counts with
  known sequencing depth — see [`fit_zibbmr()`].

Both models and their estimation algorithm were originally developed by John
Barrera:

- ZIBR: <https://github.com/jbarrera232/saem-zibr>
- ZIBBMR: <https://github.com/jbarrera232/saem-zibbmr>

This package organizes, documents and tests that implementation for general
use in longitudinal microbiome analysis.

## Installation

```r
# install.packages("remotes")
remotes::install_github("gabrielagutierrezbernal/saemMicrobiome")
```

## Example: ZIBR (proportions)

> **Note:** the example uses 300 subjects, a sample size large enough for the
> estimates to land close to the true values used in the simulation
> (`alpha = c(-0.3, 0.5)`, `beta = c(0.2, -0.4)`, `phi = 15`). With smaller
> samples, the estimates from a single fit carry more sampling noise, which
> decreases as the number of subjects grows.


``` r
library(saemMicrobiome)

# You only need to change n_subjects / n_time; n_obs is derived from them.
n_subjects <- 300
n_time <- 4
n_obs <- n_subjects * n_time

set.seed(3)
dat <- simulate_zibr_data(
  n_subjects = n_subjects, n_time = n_time,
  X = matrix(rbinom(n_obs, 1, 0.5)), Z = matrix(rbinom(n_obs, 1, 0.5)),
  alpha = c(-0.3, 0.5), beta = c(0.2, -0.4),
  sigma_alpha = 0.4, sigma_beta = 0.3, phi = 15, seed = 3
)

fit <- fit_zibr(
  y = dat$Y, id = dat$Subject, X = dat$X.1, Z = dat$Z.1,
  phi_start = 10, alpha_start = c(-0.2, 0.1), beta_start = c(0.1, 0.1),
  n_iter = 500, seed = 1, compute_fim = FALSE
)

# The estimates land close to the true values (-0.3, 0.5, 0.2, -0.4, 15)
print(fit)
#> ===== Results SAEM-ZIBR =====
#> == Logistic part: p_it ==
#>             Estimate   Type
#> Intercept -0.3452480 Random
#> X.1        0.5332735  Fixed
#> == Beta part: u_it ==
#>             Estimate   Type
#> Intercept  0.2656632 Random
#> Z.1       -0.5009376  Fixed
#> === Random-effect variances ===
#> == Logistic part ==
#>            Variance  sqrt.Var
#> Intercept 0.1999159 0.4471195
#> == Beta part ==
#>            Variance  sqrt.Var
#> Intercept 0.1034754 0.3216759
#> === Phi: 15.63579
#> === Marginal log-likelihood (importance sampling): -493.249
```

## Example: ZIBBMR (counts with sequencing depth)

> **Note:** as in the ZIBR example, 300 subjects are used, a sample size large
> enough for the estimates to land close to the true values used in the
> simulation (`alpha = c(-0.3, 0.5)`, `beta = c(0.2, -0.4)`, `phi = 15`). With
> smaller samples, the estimates from a single fit carry more sampling noise,
> which decreases as the number of subjects grows.


``` r
# Reuses n_subjects / n_time / n_obs from the previous example.
S <- rep(1000, n_obs)
set.seed(3)
dat_counts <- simulate_zibbmr_data(
  n_subjects = n_subjects, n_time = n_time, S = S,
  X = matrix(rbinom(n_obs, 1, 0.5)), Z = matrix(rbinom(n_obs, 1, 0.5)),
  alpha = c(-0.3, 0.5), beta = c(0.2, -0.4),
  sigma_alpha = 0.4, sigma_beta = 0.3, phi = 15, seed = 3
)

fit_counts <- fit_zibbmr(
  y = dat_counts$Y, S = dat_counts$TotalCounts, id = dat_counts$Subject,
  X = dat_counts$X.1, Z = dat_counts$Z.1,
  phi_start = 10, alpha_start = c(-0.2, 0.1), beta_start = c(0.1, 0.1),
  n_iter = 500, seed = 1, compute_fim = FALSE
)

# The estimates land close to the true values (-0.3, 0.5, 0.2, -0.4, 15)
print(fit_counts)
#> ===== Results SAEM-ZIBBMR =====
#> == Logistic part: p_it ==
#>             Estimate   Type
#> Intercept -0.3165230 Random
#> X.1        0.5434143  Fixed
#> == Beta-binomial part: u_it ==
#>            Estimate   Type
#> Intercept  0.227703 Random
#> Z.1       -0.444607  Fixed
#> === Random-effect variances ===
#> == Logistic part ==
#>            Variance  sqrt.Var
#> Intercept 0.1338111 0.3658019
#> == Beta-binomial part ==
#>            Variance  sqrt.Var
#> Intercept 0.1277649 0.3574422
#> === Phi: 16.97994
#> === Marginal log-likelihood (importance sampling): -4542.437
```

## Per-taxon fit and nested-model comparison


``` r
sim <- simulate_microbiome_data(n_ind = 15, n_time = 4, n_taxa = 3, seed = 1)

full <- fit_zibr_taxon(
  data = sim$proportion, taxon = "Taxon1", id = "id",
  covariates = c("time", "group"), n_iter = 50, seed = 1
)
reduced <- fit_zibr_taxon(
  data = sim$proportion, taxon = "Taxon1", id = "id",
  covariates = "time", n_iter = 50, seed = 1
)

lrt_zibr(full, reduced, df = 1)
#>    LL_full LL_reduced      LRT df   p_value
#> 1 14.80238   13.75381 2.097151  1 0.1475739
```

## More information

See `vignette("get-started", package = "saemMicrobiome")` for a more complete
introduction, including when to use ZIBR vs. ZIBBMR and how to prepare your own
data.
