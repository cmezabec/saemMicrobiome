# saemMicrobiome 0.0.2

## Changes to the ZIBBMR estimator (`fit_zibbmr()`)

* **Simulated annealing of the random-effect variances is now on by default**
  (`annealing = TRUE`). During the first half of the burn-in no variance may
  decrease by more than a factor `annealing_tau` (0.97) per iteration, and the
  random effects start with variance at least 1, as in saemix and Monolix.
  Without it, the variance of a random intercept about which each subject
  carries little information (for instance the zero-inflation intercept with
  three observations per subject, mostly zeros) could collapse to zero, an
  absorbing state of the algorithm that is not the maximum likelihood estimate.
  In simulations (Setting 3 of the accompanying paper, T = 3) the collapse
  occurred in 56-70% of the data sets and biased the zero-inflation intercept
  by about +1; with annealing the exact marginal log-likelihood at the SAEM
  estimate matches that at the Laplace estimate, and the median bias of every
  fixed effect is at most that of glmmTMB. Use `annealing = FALSE` to
  reproduce the original algorithm.

* **The default number of iterations is now 2000** (was 1000), also in
  `fit_zibbmr_taxon()` and `fit_zibbmr_taxa()`. With a covariate that is
  constant within subject and a large random-effect variance, the parameters
  were still drifting at the end of the 750-iteration burn-in; with 2000 the
  remaining bias is negligible.

* New option `mstep = "score"`: after the burn-in, the coefficients without a
  random effect and `phi` are updated by a Newton step on the score with a
  stochastic-approximation average of the Hessian (Gu and Kong, 1998) instead
  of averaging the per-iteration maximizers. It gives the same results as the
  default (`"argmax"`) in every case examined so far and is kept as an option.

* `fit_zibr()` is unchanged.

# saemMicrobiome 0.0.1

* Initial version.
