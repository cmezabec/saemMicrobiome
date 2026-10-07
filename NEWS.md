# saemMicrobiome 0.0.3

* The logistic probabilities of both model parts are kept strictly inside
  (0, 1). With very few positive counts (for instance a single positive count
  in a data set), the M-step on one simulated sample can be separated during
  the burn-in and a fixed effect can diverge; `plogis()` then returned exactly
  0 or 1, `lgamma()` returned `Inf`, and the fit stopped with an error. Only
  those two values are moved to the nearest representable probability, so fits
  in which no linear predictor goes beyond about +-37 are unchanged. Such data
  sets carry no information on the abundance part: the fit now ends at the
  boundary (very large `phi` or slopes) instead of failing, as glmmTMB does.

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

* The stochastic Fisher information is now accumulated only during the
  averaging phase. During the burn-in its running averages were overwritten at
  every iteration, so the result is identical; with the default burn-in the
  fit is about 40% faster when `compute_fim = TRUE`.

* `fit_zibr()` is unchanged.

# saemMicrobiome 0.0.1

* Initial version.
