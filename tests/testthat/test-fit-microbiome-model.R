test_that("fit_microbiome_model dispatches to zibbmr by default", {
  set.seed(1)
  sim <- simulate_microbiome_data(n_ind = 10, n_time = 3, n_taxa = 2, N = 500, seed = 1)

  fit <- fit_microbiome_model(
    model = "zibbmr", data = sim$count, taxon = "Taxon1",
    id = "id", total = "N", covariates = c("time", "group"),
    iter = 20, seed = 1
  )

  expect_s3_class(fit, "zibbmr_saem")
})

test_that("fit_microbiome_model dispatches to zibr when requested", {
  set.seed(1)
  sim <- simulate_microbiome_data(n_ind = 10, n_time = 3, n_taxa = 2, N = 500, seed = 1)

  fit <- fit_microbiome_model(
    model = "zibr", data = sim$proportion, taxon = "Taxon1",
    id = "id", covariates = c("time", "group"),
    iter = 20, seed = 1
  )

  expect_s3_class(fit, "zibr_saem")
})

test_that("fit_saem_microbiome is an identical alias of fit_microbiome_model", {
  sim <- simulate_microbiome_data(n_ind = 10, n_time = 3, n_taxa = 2, N = 500, seed = 1)

  via_alias <- fit_saem_microbiome(
    model = "zibbmr", data = sim$count, taxon = "Taxon1",
    id = "id", total = "N", covariates = c("time", "group"),
    iter = 20, seed = 1
  )
  via_original <- fit_microbiome_model(
    model = "zibbmr", data = sim$count, taxon = "Taxon1",
    id = "id", total = "N", covariates = c("time", "group"),
    iter = 20, seed = 1
  )

  expect_equal(via_alias$mu, via_original$mu)
  expect_equal(via_alias$loglik, via_original$loglik)
})

test_that("fit_microbiome_model requires a valid model", {
  sim <- simulate_microbiome_data(n_ind = 5, n_time = 2, n_taxa = 1, seed = 1)
  expect_error(
    fit_microbiome_model(model = "other", data = sim$count, taxon = "Taxon1"),
    "arg"
  )
})

test_that("fit_microbiome_model validates that taxon/id/total exist in data", {
  sim <- simulate_microbiome_data(n_ind = 5, n_time = 2, n_taxa = 1, seed = 1)

  expect_error(
    fit_microbiome_model(model = "zibbmr", data = sim$count, taxon = "DoesNotExist"),
    "taxon"
  )
  expect_error(
    fit_microbiome_model(model = "zibbmr", data = sim$count, taxon = "Taxon1", id = "DoesNotExist"),
    "id column"
  )
  expect_error(
    fit_microbiome_model(model = "zibbmr", data = sim$count, taxon = "Taxon1", total = "DoesNotExist"),
    "total column"
  )
})

test_that("fit_microbiome_model allows different covariates for X and Z", {
  sim <- simulate_microbiome_data(n_ind = 10, n_time = 3, n_taxa = 2, N = 500, seed = 1)

  fit <- fit_microbiome_model(
    model = "zibbmr", data = sim$count, taxon = "Taxon1",
    id = "id", total = "N", x_covariates = "time", z_covariates = "group",
    iter = 20, seed = 1
  )

  expect_s3_class(fit, "zibbmr_saem")
  expect_length(fit$mu, 4)
})

test_that("fit_microbiome_model works with zi = FALSE (fixed bug)", {
  # simulate_microbiome_data() with n_taxa = 1 degenerates (a single
  # multinomial category always gets the full total, with no variance), so
  # this test uses simulate_zibbmr_data(), which does generate real variation
  # via the beta-binomial. Before this fix, fit_microbiome_model built X the
  # same way with zi = TRUE or FALSE, and fit_zibbmr failed because it does not
  # accept X when zi = FALSE.
  S <- rep(500, 30)
  sim <- simulate_zibbmr_data(
    n_subjects = 10, n_time = 3, S = S, zi = FALSE,
    Z = matrix(rbinom(30, 1, 0.5)), beta = c(0.2, -0.3), sigma_beta = 0.3,
    phi = 15, seed = 8
  )
  dat <- data.frame(
    id = as.numeric(factor(sim$Subject, levels = unique(sim$Subject))),
    time = sim$Z.1, Taxon1 = sim$Y, N = S
  )

  fit <- fit_microbiome_model(
    model = "zibbmr", data = dat, taxon = "Taxon1",
    id = "id", total = "N", covariates = "time", zi = FALSE,
    iter = 30, seed = 1
  )

  expect_s3_class(fit, "zibbmr_saem")
  expect_length(fit$mu, 2)
})

test_that("fit_microbiome_model draws starting values when none are supplied", {
  sim <- simulate_microbiome_data(n_ind = 10, n_time = 3, n_taxa = 3, N = 500, seed = 1)

  fit1 <- fit_microbiome_model(
    model = "zibbmr", data = sim$count, taxon = "Taxon1",
    id = "id", total = "N", covariates = "time", iter = 15, seed = 11
  )
  fit2 <- fit_microbiome_model(
    model = "zibbmr", data = sim$count, taxon = "Taxon1",
    id = "id", total = "N", covariates = "time", iter = 15, seed = 99
  )

  # different seeds -> different starting values -> different traces
  expect_false(isTRUE(all.equal(fit1$trace[1, ], fit2$trace[1, ])))
})
