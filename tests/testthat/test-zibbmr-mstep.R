test_that("mstep only accepts the two documented values", {
  dat <- simulate_zibbmr_data(
    n_subjects = 10, n_time = 3, S = rep(1000, 30),
    alpha = c(-0.3, 0.5), beta = c(0.2, -0.4),
    sigma_alpha = 0.4, sigma_beta = 0.3, phi = 15,
    X = matrix(rbinom(30, 1, 0.5)), Z = matrix(rbinom(30, 1, 0.5)), seed = 1
  )
  expect_error(
    fit_zibbmr(y = dat$Y, S = dat$TotalCounts, id = dat$Subject,
               X = dat$X.1, Z = dat$Z.1, phi_start = 10,
               alpha_start = c(-0.2, 0.1), beta_start = c(0.1, 0.1),
               n_iter = 20, seed = 1, compute_fim = FALSE, mstep = "other"),
    "should be one of"
  )
})

test_that("the score step is only used after the burn-in", {
  # With the burn-in covering every iteration but the last ones, the two M-steps
  # must follow exactly the same trajectory up to the end of the burn-in: the
  # same seed, the same random draws, the same argmax updates.
  dat <- simulate_zibbmr_data(
    n_subjects = 20, n_time = 4, S = rep(1000, 80),
    alpha = c(-0.3, 0.5), beta = c(0.2, -0.4),
    sigma_alpha = 0.4, sigma_beta = 0.3, phi = 15,
    X = matrix(rep(c(0, 1), each = 4, length.out = 80)),
    Z = matrix(rep(c(0, 1), each = 4, length.out = 80)), seed = 3
  )
  args <- list(y = dat$Y, S = dat$TotalCounts, id = dat$Subject,
               X = dat$X.1, Z = dat$Z.1, phi_start = 10,
               alpha_start = c(-0.2, 0.1), beta_start = c(0.1, 0.1),
               n_iter = 40, seed = 11, compute_fim = FALSE)
  fa <- do.call(fit_zibbmr, c(args, mstep = "argmax"))
  fs <- do.call(fit_zibbmr, c(args, mstep = "score"))
  burn_in <- floor(0.75 * 40)
  expect_equal(fa$trace[seq_len(burn_in), ], fs$trace[seq_len(burn_in), ])
  expect_identical(fs$mstep, "score")
  expect_true(all(is.finite(c(fs$mu, fs$phi, fs$loglik))))
  expect_gt(fs$phi, 0)
})

test_that("both M-steps agree when every iteration is informative", {
  # Many observations per subject and chains: the per-iteration sample is
  # precise, the bias of averaging the argmax is negligible, and the two
  # schemes must land on essentially the same estimates.
  skip_on_cran()
  n_subjects <- 40
  n_time <- 8
  n_total <- n_subjects * n_time
  X <- rep(c(0, 1), each = n_time, length.out = n_total)
  dat <- simulate_zibbmr_data(
    n_subjects = n_subjects, n_time = n_time, S = rep(1000, n_total),
    X = matrix(X, ncol = 1), Z = matrix(X, ncol = 1),
    alpha = c(-0.3, 0.5), beta = c(0.2, -0.4),
    sigma_alpha = 0.4, sigma_beta = 0.3, phi = 15, seed = 5
  )
  args <- list(y = dat$Y, S = dat$TotalCounts, id = dat$Subject,
               X = dat$X.1, Z = dat$Z.1, phi_start = 10,
               alpha_start = c(-0.2, 0.1), beta_start = c(0.1, 0.1),
               n_iter = 400, n_chains = 10, seed = 2, compute_fim = FALSE)
  fa <- do.call(fit_zibbmr, c(args, mstep = "argmax"))
  fs <- do.call(fit_zibbmr, c(args, mstep = "score"))
  expect_lt(max(abs(fa$mu - fs$mu)), 0.15)
  expect_lt(abs(log(fa$phi) - log(fs$phi)), 0.15)
})

test_that("annealing keeps each variance from falling faster than tau", {
  dat <- simulate_zibbmr_data(
    n_subjects = 20, n_time = 3, S = rep(1000, 60),
    alpha = c(-0.3, 0.5), beta = c(0.2, -0.4),
    sigma_alpha = 0.4, sigma_beta = 0.3, phi = 15,
    X = matrix(rep(c(0, 1), each = 3, length.out = 60)),
    Z = matrix(rep(c(0, 1), each = 3, length.out = 60)), seed = 4
  )
  args <- list(y = dat$Y, S = dat$TotalCounts, id = dat$Subject,
               X = dat$X.1, Z = dat$Z.1, phi_start = 10,
               alpha_start = c(-0.2, 0.1), beta_start = c(0.1, 0.1),
               n_iter = 80, seed = 9, compute_fim = FALSE)
  fit <- do.call(fit_zibbmr, c(args, annealing = TRUE, annealing_tau = 0.9))
  v <- fit$trace[, grep("^var", colnames(fit$trace)), drop = FALSE]
  n_anneal <- floor(floor(0.75 * 80) / 2)
  # From iteration 11 (first update of G) to the end of the annealing phase,
  # no variance may drop below tau times its previous value.
  ratio <- v[12:n_anneal, , drop = FALSE] / v[11:(n_anneal - 1), , drop = FALSE]
  expect_true(all(ratio >= 0.9 - 1e-12))
  expect_equal(unname(fit$annealing), c(0.9, n_anneal))
  # Annealing is the default since 0.0.2; annealing = FALSE turns it off.
  expect_false(is.null(do.call(fit_zibbmr, args)$annealing))
  expect_null(do.call(fit_zibbmr, c(args, annealing = FALSE))$annealing)
})
