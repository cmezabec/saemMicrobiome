mock_fit <- function(loglik, class) {
  structure(list(loglik = loglik), class = class)
}

test_that("lrt_zibr calcula el estadistico y p-valor correctos", {
  full <- mock_fit(-100, "zibr_saem")
  reduced <- mock_fit(-105, "zibr_saem")

  res <- lrt_zibr(full, reduced, df = 2)

  expect_equal(res$LL_full, -100)
  expect_equal(res$LL_reduced, -105)
  expect_equal(res$LRT, 2 * (-100 - (-105)))
  expect_equal(res$p_value, stats::pchisq(res$LRT, df = 2, lower.tail = FALSE))
})

test_that("lrt_zibbmr calcula el estadistico y p-valor correctos", {
  full <- mock_fit(-200, "zibbmr_saem")
  reduced <- mock_fit(-210, "zibbmr_saem")

  res <- lrt_zibbmr(full, reduced, df = 1)

  expect_equal(res$LRT, 2 * (-200 - (-210)))
  expect_equal(res$df, 1)
})

test_that("lrt_zibr_table arma una fila por especie y marca Detected correctamente", {
  full_models <- list(sp1 = mock_fit(-100, "zibr_saem"), sp2 = mock_fit(-50, "zibr_saem"))
  reduced_models <- list(sp1 = mock_fit(-105, "zibr_saem"), sp2 = mock_fit(-50.001, "zibr_saem"))

  tab <- lrt_zibr_table(full_models, reduced_models, df = 2, alpha = 0.05)

  expect_equal(nrow(tab), 2)
  expect_equal(tab$Species, c("sp1", "sp2"))
  expect_true(tab$Detected[1])
  expect_false(tab$Detected[2])
})

test_that("lrt_zibr_table exige listas de la misma longitud", {
  expect_error(
    lrt_zibr_table(
      list(mock_fit(-1, "zibr_saem")),
      list(mock_fit(-1, "zibr_saem"), mock_fit(-2, "zibr_saem"))
    ),
    "misma longitud"
  )
})

test_that("lrt_zibr acepta cualquier objeto con metodo logLik() (no solo zibr_saem)", {
  full <- lm(mpg ~ wt, data = mtcars)
  reduced <- mock_fit(as.numeric(stats::logLik(full)) - 5, "zibr_saem")

  res <- lrt_zibr(full, reduced, df = 1)

  expect_equal(res$LL_full, as.numeric(stats::logLik(full)))
  expect_equal(res$LRT, 2 * 5)
})

test_that("zibr_results_table arma las tres comparaciones LRT", {
  species <- c("sp1", "sp2")

  tab <- zibr_results_table(
    species = species,
    mod1_full = list(mock_fit(-90, "zibr_saem"), mock_fit(-40, "zibr_saem")),
    mod1_no_preg = list(mock_fit(-100, "zibr_saem"), mock_fit(-40.001, "zibr_saem")),
    mod2_full = list(mock_fit(-90, "zibr_saem"), mock_fit(-40, "zibr_saem")),
    mod2_no_preg = list(mock_fit(-100, "zibr_saem"), mock_fit(-40.001, "zibr_saem")),
    mod2_no_inter = list(mock_fit(-100, "zibr_saem"), mock_fit(-40.001, "zibr_saem")),
    df = 2, alpha = 0.05
  )

  expect_equal(nrow(tab), 2)
  expect_equal(tab$Species, species)
  expect_true(tab$Detec_Preg1[1])
  expect_false(tab$Detec_Preg1[2])
})

test_that("zibbmr_results_table arma las tres comparaciones LRT", {
  species <- c("sp1", "sp2")

  tab <- zibbmr_results_table(
    species = species,
    mod1_full = list(mock_fit(-90, "zibbmr_saem"), mock_fit(-40, "zibbmr_saem")),
    mod1_no_preg = list(mock_fit(-100, "zibbmr_saem"), mock_fit(-40.001, "zibbmr_saem")),
    mod2_full = list(mock_fit(-90, "zibbmr_saem"), mock_fit(-40, "zibbmr_saem")),
    mod2_no_preg = list(mock_fit(-100, "zibbmr_saem"), mock_fit(-40.001, "zibbmr_saem")),
    mod2_no_inter = list(mock_fit(-100, "zibbmr_saem"), mock_fit(-40.001, "zibbmr_saem")),
    df = 2, alpha = 0.05
  )

  expect_equal(nrow(tab), 2)
  expect_true(tab$Detec_Preg1[1])
  expect_false(tab$Detec_Preg1[2])
})
