# Regression: the covariance between the random effects of the two model parts
# (zero inflation and beta-binomial / beta part) must be estimated.
#
# Before the fix, .saem_diag_inverse() computed diag(1 / diag(G)), i.e. it
# inverted G discarding the off-diagonal terms. Since that precision is used in
# the Metropolis-Hastings acceptance ratio, the sampler targeted a model with
# independent random effects no matter what; the accumulated covariance never
# grew and rho stayed trapped near 0 (with a true rho of 0.7 it estimated
# ~0.09).

test_that(".saem_diag_inverse inverts the full matrix, not just the diagonal", {
  G <- matrix(c(0.49, 0.21, 0.21, 0.25), 2, 2)

  expect_equal(.saem_diag_inverse(G), solve(G))
  # The previous version returned this, which is NOT the inverse of G:
  expect_false(isTRUE(all.equal(.saem_diag_inverse(G), diag(1 / diag(G)))))
  # G %*% G^{-1} must give the identity.
  expect_equal(G %*% .saem_diag_inverse(G), diag(2), ignore_attr = TRUE)
})

test_that(".saem_diag_inverse is still correct with diagonal G", {
  # Backward compatibility: all published results use diagonal G, and there the
  # two versions agree exactly.
  G <- diag(c(0.49, 0.25))
  expect_equal(.saem_diag_inverse(G), diag(1 / diag(G)))

  # 1x1 case (models without a zero-inflation part).
  expect_equal(.saem_diag_inverse(matrix(0.4, 1, 1)), matrix(2.5, 1, 1))
})

test_that(".saem_diag_inverse does not fail with a singular G", {
  # This can happen in early SAEM iterations, before the variance components
  # stabilize: it must fall back to the pseudo-inverse.
  G <- matrix(c(1, 1, 1, 1), 2, 2)
  expect_no_error(inv <- .saem_diag_inverse(G))
  expect_true(all(is.finite(inv)))
})

test_that("fit_zibbmr recovers a non-zero correlation between the random effects", {
  skip_on_cran()
  skip_if_not_installed("MASS")

  set.seed(20260805)
  N <- 100; T <- 10; n <- N * T
  rho <- 0.7; s1 <- 0.7; s2 <- 0.5; phi <- 6.4

  G <- matrix(c(s1^2, rho * s1 * s2, rho * s1 * s2, s2^2), 2, 2)
  re <- MASS::mvrnorm(N, c(-0.5, -0.5), G)
  id <- rep(seq_len(N), each = T)
  x <- rep(c(rep(0, (N %/% 2) * T), rep(1, (N - N %/% 2) * T)), length.out = n)
  S <- sample(200:800, n, TRUE)
  p <- plogis(re[id, 1] + 0.5 * x)
  u <- plogis(re[id, 2] + 0.5 * x)
  Y <- rbinom(n, S, rbeta(n, u * phi, (1 - u) * phi)) * rbinom(n, 1, p)

  fit <- fit_zibbmr(
    y = Y, S = S, id = id, X = x, Z = x,
    phi_start = 18, alpha_start = c(0.1, 0.1), beta_start = c(0.1, 0.1),
    n_iter = 1000, n_chains = 5, seed = 1, compute_fim = FALSE,
    cov_random = "unstructured"
  )

  rho_hat <- fit$G[1, 2] / sqrt(fit$G[1, 1] * fit$G[2, 2])

  # With the bug, rho_hat hovered around 0.09. The 0.3 threshold is loose on
  # purpose: it separates "the algorithm sees the correlation" from "it does
  # not" without depending on the exact precision on a single dataset. Over 16
  # replicates with N=100 and T=10 the estimator is unbiased: mean rho_hat 0.746
  # against a true value of 0.70, with a Monte Carlo error of 0.043.
  expect_gt(rho_hat, 0.3)
  expect_lt(rho_hat, 1)
})

test_that("fit_zibbmr does not invent correlation when the effects are independent", {
  skip_on_cran()
  skip_if_not_installed("MASS")

  set.seed(20260806)
  N <- 100; T <- 10; n <- N * T
  s1 <- 0.7; s2 <- 0.5; phi <- 6.4

  re <- cbind(rnorm(N, -0.5, s1), rnorm(N, -0.5, s2))
  id <- rep(seq_len(N), each = T)
  x <- rep(c(rep(0, (N %/% 2) * T), rep(1, (N - N %/% 2) * T)), length.out = n)
  S <- sample(200:800, n, TRUE)
  p <- plogis(re[id, 1] + 0.5 * x)
  u <- plogis(re[id, 2] + 0.5 * x)
  Y <- rbinom(n, S, rbeta(n, u * phi, (1 - u) * phi)) * rbinom(n, 1, p)

  fit <- fit_zibbmr(
    y = Y, S = S, id = id, X = x, Z = x,
    phi_start = 18, alpha_start = c(0.1, 0.1), beta_start = c(0.1, 0.1),
    n_iter = 1000, n_chains = 5, seed = 1, compute_fim = FALSE,
    cov_random = "unstructured"
  )

  rho_hat <- fit$G[1, 2] / sqrt(fit$G[1, 1] * fit$G[2, 2])
  expect_lt(abs(rho_hat), 0.3)
})

test_that("cov_random = \"diag\" forces the covariance to zero", {
  skip_on_cran()

  set.seed(20260807)
  N <- 40; T <- 5; n <- N * T
  id <- rep(seq_len(N), each = T)
  x <- rep(c(0, 1), each = T, length.out = n)
  S <- rep(1000, n)
  dat <- simulate_zibbmr_data(
    n_subjects = N, n_time = T, S = S,
    X = matrix(x, ncol = 1), Z = matrix(x, ncol = 1),
    alpha = c(-0.3, 0.5), beta = c(0.2, -0.4),
    sigma_alpha = 0.4, sigma_beta = 0.3, phi = 15, seed = 3
  )

  fit <- fit_zibbmr(
    y = dat$Y, S = dat$TotalCounts, id = dat$Subject,
    X = matrix(x, ncol = 1), Z = matrix(x, ncol = 1),
    phi_start = 10, alpha_start = c(-0.2, 0.1), beta_start = c(0.1, 0.1),
    n_iter = 300, n_chains = 5, seed = 5, compute_fim = FALSE
  )

  # It is the default and the specification of the paper.
  expect_identical(fit$cov_random, "diag")
  expect_equal(fit$G[1, 2], 0)
  expect_equal(fit$G[2, 1], 0)
})
