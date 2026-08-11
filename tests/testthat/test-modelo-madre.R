test_that("ajustar_modelo_microbioma despacha a zibbmr por defecto", {
  set.seed(1)
  sim <- simular_datos_microbioma(n_ind = 10, n_time = 3, n_taxa = 2, N = 500, seed = 1)

  fit <- ajustar_modelo_microbioma(
    modelo = "zibbmr", datos = sim$conteo, taxon = "Taxon1",
    id = "id", total = "N", covariables = c("tiempo", "grupo"),
    iter = 20, seed = 1
  )

  expect_s3_class(fit, "zibbmr_saem")
})

test_that("ajustar_modelo_microbioma despacha a zibr cuando se pide", {
  set.seed(1)
  sim <- simular_datos_microbioma(n_ind = 10, n_time = 3, n_taxa = 2, N = 500, seed = 1)

  fit <- ajustar_modelo_microbioma(
    modelo = "zibr", datos = sim$proporcion, taxon = "Taxon1",
    id = "id", covariables = c("tiempo", "grupo"),
    iter = 20, seed = 1
  )

  expect_s3_class(fit, "zibr_saem")
})

test_that("fit_saem_microbiome es un alias identico a ajustar_modelo_microbioma", {
  sim <- simular_datos_microbioma(n_ind = 10, n_time = 3, n_taxa = 2, N = 500, seed = 1)

  via_alias <- fit_saem_microbiome(
    modelo = "zibbmr", datos = sim$conteo, taxon = "Taxon1",
    id = "id", total = "N", covariables = c("tiempo", "grupo"),
    iter = 20, seed = 1
  )
  via_original <- ajustar_modelo_microbioma(
    modelo = "zibbmr", datos = sim$conteo, taxon = "Taxon1",
    id = "id", total = "N", covariables = c("tiempo", "grupo"),
    iter = 20, seed = 1
  )

  expect_equal(via_alias$mu, via_original$mu)
  expect_equal(via_alias$loglik, via_original$loglik)
})

test_that("ajustar_modelo_microbioma exige un modelo valido", {
  sim <- simular_datos_microbioma(n_ind = 5, n_time = 2, n_taxa = 1, seed = 1)
  expect_error(
    ajustar_modelo_microbioma(modelo = "otro", datos = sim$conteo, taxon = "Taxon1"),
    "arg"
  )
})

test_that("ajustar_modelo_microbioma valida que existan taxon/id/total en datos", {
  sim <- simular_datos_microbioma(n_ind = 5, n_time = 2, n_taxa = 1, seed = 1)

  expect_error(
    ajustar_modelo_microbioma(modelo = "zibbmr", datos = sim$conteo, taxon = "NoExiste"),
    "taxon"
  )
  expect_error(
    ajustar_modelo_microbioma(modelo = "zibbmr", datos = sim$conteo, taxon = "Taxon1", id = "NoExiste"),
    "columna id"
  )
  expect_error(
    ajustar_modelo_microbioma(modelo = "zibbmr", datos = sim$conteo, taxon = "Taxon1", total = "NoExiste"),
    "columna total"
  )
})

test_that("ajustar_modelo_microbioma admite covariables distintas para X y Z", {
  sim <- simular_datos_microbioma(n_ind = 10, n_time = 3, n_taxa = 2, N = 500, seed = 1)

  fit <- ajustar_modelo_microbioma(
    modelo = "zibbmr", datos = sim$conteo, taxon = "Taxon1",
    id = "id", total = "N", x_covariables = "tiempo", z_covariables = "grupo",
    iter = 20, seed = 1
  )

  expect_s3_class(fit, "zibbmr_saem")
  expect_length(fit$mu, 4)
})

test_that("ajustar_modelo_microbioma funciona con zi = FALSE (bug corregido)", {
  # simular_datos_microbioma() con n_taxa = 1 degenera (una sola categoria
  # multinomial siempre recibe el total completo, sin varianza), asi que
  # para este test se usa simulate_zibbmr_data(), que si genera variacion
  # real via beta-binomial. Antes de este arreglo, ajustar_modelo_microbioma
  # armaba X igual con zi = TRUE o FALSE, y fit_zibbmr fallaba porque no
  # acepta X cuando zi = FALSE.
  S <- rep(500, 30)
  sim <- simulate_zibbmr_data(
    n_subjects = 10, n_time = 3, S = S, zi = FALSE,
    Z = matrix(rbinom(30, 1, 0.5)), beta = c(0.2, -0.3), sigma_beta = 0.3,
    phi = 15, seed = 8
  )
  dat <- data.frame(
    id = as.numeric(factor(sim$Subject, levels = unique(sim$Subject))),
    tiempo = sim$Z.1, Taxon1 = sim$Y, N = S
  )

  fit <- ajustar_modelo_microbioma(
    modelo = "zibbmr", datos = dat, taxon = "Taxon1",
    id = "id", total = "N", covariables = "tiempo", zi = FALSE,
    iter = 30, seed = 1
  )

  expect_s3_class(fit, "zibbmr_saem")
  expect_length(fit$mu, 2)
})

test_that("ajustar_modelo_microbioma sortea valores iniciales cuando no se entregan", {
  sim <- simular_datos_microbioma(n_ind = 10, n_time = 3, n_taxa = 3, N = 500, seed = 1)

  fit1 <- ajustar_modelo_microbioma(
    modelo = "zibbmr", datos = sim$conteo, taxon = "Taxon1",
    id = "id", total = "N", covariables = "tiempo", iter = 15, seed = 11
  )
  fit2 <- ajustar_modelo_microbioma(
    modelo = "zibbmr", datos = sim$conteo, taxon = "Taxon1",
    id = "id", total = "N", covariables = "tiempo", iter = 15, seed = 99
  )

  # semillas distintas -> valores iniciales distintos -> trazas distintas
  expect_false(isTRUE(all.equal(fit1$trace[1, ], fit2$trace[1, ])))
})
