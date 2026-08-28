# External validation of fit_zibr() against the ZIBR package of Chen and Li
# (github.com/chvlyl/ZIBR), which fits the same model by Gauss-Hermite
# quadrature rather than by SAEM. Independent code and an independent
# estimation method, so agreement is evidence that the estimator is correct
# rather than merely self-consistent.
#
# The data are `ibd` from that package, the data set used in their paper.
#
# ZIBR is not on CRAN, so it is not declared as a dependency; the test skips
# when it is absent. It is also slow (a few minutes), hence skip_on_cran().
#
# Tolerances are deliberately loose: the two methods are not expected to agree
# to machine precision, only to well within the uncertainty of the estimates.
# They are set roughly five times the differences observed when this test was
# written, so a genuine regression is caught while Monte Carlo variation is not.

test_that("fit_zibr agrees with the ZIBR package on the ibd data", {
  skip_on_cran()
  skip_if_not_installed("ZIBR")

  utils::data("ibd", package = "ZIBR", envir = environment())
  X <- matrix(ibd$Treatment, ncol = 1, dimnames = list(NULL, "Treatment"))

  ref <- ZIBR::zibr(
    logistic_cov = X, beta_cov = X, Y = ibd$Abundance,
    subject_ind = ibd$Subject, time_ind = ibd$Time, verbose = FALSE
  )

  fit <- fit_zibr(
    y = ibd$Abundance, id = ibd$Subject, X = X, Z = X, zi = TRUE,
    phi_start = 10, alpha_start = c(0.1, 0.1), beta_start = c(0.1, 0.1),
    n_iter = 1000, n_chains = 10, seed = 1, compute_fim = FALSE
  )

  # Abundance component: both coefficients are well determined here, so these
  # are the informative comparisons.
  expect_equal(as.numeric(fit$mu[3]), as.numeric(ref$beta_est_table[1, 1]),
               tolerance = 0.15)
  expect_equal(as.numeric(fit$mu[4]), as.numeric(ref$beta_est_table[2, 1]),
               tolerance = 0.15)

  # Variance components and dispersion.
  expect_equal(sqrt(fit$G[1, 1]), as.numeric(ref$logistic_s1_est), tolerance = 0.35)
  expect_equal(sqrt(fit$G[2, 2]), as.numeric(ref$beta_s2_est),     tolerance = 0.10)
  expect_equal(fit$phi,           as.numeric(ref$beta_v_est),      tolerance = 1.0)

  # The intercept of the zero-inflation component.
  expect_equal(as.numeric(fit$mu[1]), as.numeric(ref$logistic_est_table[1, 1]),
               tolerance = 0.7)

  # The marginal log-likelihood is the most stringent check: the two packages
  # compute it by entirely different means, quadrature against importance
  # sampling, and it summarises the whole fit.
  expect_equal(fit$loglik, ref$loglikelihood, tolerance = 2.0)

  # The Treatment coefficient of the zero-inflation component is deliberately
  # not compared. ZIBR itself reports p = 0.89 for it on these data, so it is
  # not distinguishable from zero and the two fits disagree on its sign while
  # agreeing on its magnitude. Testing it would only encode noise.
})
