# Regresion: la covarianza entre los efectos aleatorios de las dos partes del
# modelo (inflacion de ceros y parte beta-binomial / beta) debe estimarse.
#
# Antes del arreglo, .saem_diag_inverse() calculaba diag(1 / diag(G)), o sea
# invertia G descartando los terminos fuera de la diagonal. Como esa precision
# se usa en la razon de aceptacion de Metropolis-Hastings, el muestreador tenia
# como objetivo un modelo con efectos aleatorios independientes pase lo que
# pase; la covarianza acumulada nunca crecia y rho quedaba atrapado cerca de 0
# (con rho verdadero 0.7 se estimaba ~0.09).

test_that(".saem_diag_inverse invierte la matriz completa, no solo la diagonal", {
  G <- matrix(c(0.49, 0.21, 0.21, 0.25), 2, 2)

  expect_equal(.saem_diag_inverse(G), solve(G))
  # La version anterior devolvia esto, que NO es la inversa de G:
  expect_false(isTRUE(all.equal(.saem_diag_inverse(G), diag(1 / diag(G)))))
  # G %*% G^{-1} tiene que dar la identidad.
  expect_equal(G %*% .saem_diag_inverse(G), diag(2), ignore_attr = TRUE)
})

test_that(".saem_diag_inverse sigue siendo correcta con G diagonal", {
  # Retrocompatibilidad: todos los resultados publicados usan G diagonal, y ahi
  # las dos versiones coinciden exactamente.
  G <- diag(c(0.49, 0.25))
  expect_equal(.saem_diag_inverse(G), diag(1 / diag(G)))

  # Caso 1x1 (modelos sin parte de inflacion de ceros).
  expect_equal(.saem_diag_inverse(matrix(0.4, 1, 1)), matrix(2.5, 1, 1))
})

test_that(".saem_diag_inverse no falla con una G singular", {
  # Puede ocurrir en iteraciones tempranas del SAEM, antes de que las
  # componentes de varianza se estabilicen: debe recurrir a la pseudo-inversa.
  G <- matrix(c(1, 1, 1, 1), 2, 2)
  expect_no_error(inv <- .saem_diag_inverse(G))
  expect_true(all(is.finite(inv)))
})

test_that("fit_zibbmr recupera una correlacion no nula entre los efectos aleatorios", {
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

  ajuste <- fit_zibbmr(
    y = Y, S = S, id = id, X = x, Z = x,
    phi_start = 18, alpha_start = c(0.1, 0.1), beta_start = c(0.1, 0.1),
    n_iter = 1000, n_chains = 5, seed = 1, compute_fim = FALSE,
    cov_random = "unstructured"
  )

  rho_estimado <- ajuste$G[1, 2] / sqrt(ajuste$G[1, 1] * ajuste$G[2, 2])

  # Con el bug, rho_estimado rondaba 0.09. El umbral de 0.3 es holgado a
  # proposito: separa "el algoritmo ve la correlacion" de "no la ve" sin
  # depender de la precision exacta en un solo dataset. Sobre 16 replicas con
  # N=100 y T=10 el estimador queda insesgado: rho_hat medio 0.746 contra un
  # verdadero de 0.70, con un error de Monte Carlo de 0.043.
  expect_gt(rho_estimado, 0.3)
  expect_lt(rho_estimado, 1)
})

test_that("fit_zibbmr no inventa correlacion cuando los efectos son independientes", {
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

  ajuste <- fit_zibbmr(
    y = Y, S = S, id = id, X = x, Z = x,
    phi_start = 18, alpha_start = c(0.1, 0.1), beta_start = c(0.1, 0.1),
    n_iter = 1000, n_chains = 5, seed = 1, compute_fim = FALSE,
    cov_random = "unstructured"
  )

  rho_estimado <- ajuste$G[1, 2] / sqrt(ajuste$G[1, 1] * ajuste$G[2, 2])
  expect_lt(abs(rho_estimado), 0.3)
})

test_that("cov_random = \"diag\" fuerza la covarianza a cero", {
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

  ajuste <- fit_zibbmr(
    y = dat$Y, S = dat$TotalCounts, id = dat$Subject,
    X = matrix(x, ncol = 1), Z = matrix(x, ncol = 1),
    phi_start = 10, alpha_start = c(-0.2, 0.1), beta_start = c(0.1, 0.1),
    n_iter = 300, n_chains = 5, seed = 5, compute_fim = FALSE
  )

  # Es el valor por defecto y la especificacion del articulo.
  expect_identical(ajuste$cov_random, "diag")
  expect_equal(ajuste$G[1, 2], 0)
  expect_equal(ajuste$G[2, 1], 0)
})
