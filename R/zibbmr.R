#### Implementacion limpia de SAEM-ZIBBMR ####

#### Utilidades internas ####
#### (.saem_check_packages, .saem_diag*, .saem_covariate_matrix, ####
####  .saem_replicate_design, .saem_linear_prob viven en R/utils.R, ####
####  compartidas con ZIBR) ####


#### Log-verosimilitudes condicionales usadas en el paso M ####
#### (.saem_neg_loglik_zero, la parte de inflacion de ceros, vive en ####
####  R/utils.R, compartida con ZIBR) ####

.zibbmr_neg_loglik_beta_binomial <- function(par, psi_chain, beta_random,
                                             z_design_chain, id_chain,
                                             is_positive_chain, y_chain, s_chain,
                                             n_alpha, n_beta, n_beta_random) {
  phi <- par[length(par)]
  n_rows <- nrow(psi_chain)

  if (n_beta_random != n_beta) {
    beta_fixed <- par[-length(par)]
    fixed_index <- which(!beta_random)

    psi_chain[, n_alpha + fixed_index] <- matrix(
      rep(beta_fixed, each = n_rows),
      ncol = length(beta_fixed),
      nrow = n_rows
    )
  }

  u <- .saem_linear_prob(
    psi_chain,
    n_alpha + seq_len(n_beta),
    id_chain,
    z_design_chain
  )

  loglik <- sum(
    lgamma(y_chain[is_positive_chain] + phi * u[is_positive_chain]) +
      lgamma(s_chain[is_positive_chain] - y_chain[is_positive_chain] +
               phi * (1 - u[is_positive_chain])) -
      lgamma(s_chain[is_positive_chain] + phi) +
      lgamma(phi) -
      lgamma(phi * u[is_positive_chain]) -
      lgamma(phi * (1 - u[is_positive_chain]))
  )

  -loglik
}


#### Importance sampling para log-verosimilitud marginal ####

.zibbmr_loglik_importance <- function(mu, G, phi, zi, y, s, id,
                                      x_design, z_design, n_alpha, n_beta,
                                      is_positive, is_zero,
                                      psi_mean, psi_var, random_index,
                                      n_random, n_samples = 500, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  G_inv <- .saem_diag_inverse(G)
  # prod(diag(G)) es el determinante solo si G es diagonal. Como el paso M
  # estima la matriz completa, hay que usar det(): de lo contrario la densidad
  # normal de los efectos aleatorios queda mal normalizada y la verosimilitud
  # marginal (y con ella el LRT) sale sesgada cuando la covarianza no es cero.
  G_det <- det(as.matrix(G))

  psi_array <- array(rep(psi_mean, n_samples), dim = c(dim(psi_mean), n_samples))
  sd_array <- array(rep(sqrt(psi_var), n_samples), dim = dim(psi_array))
  t_draws <- array(rt(prod(dim(psi_array)), df = 5), dim = dim(psi_array))
  psi_draws <- psi_array + sd_array * t_draws

  u_draws <- apply(
    psi_draws,
    3,
    .saem_linear_prob,
    cols = n_alpha + seq_len(n_beta),
    id = id,
    design = z_design
  )

  if (zi) {
    p_draws <- apply(
      psi_draws,
      3,
      .saem_linear_prob,
      cols = seq_len(n_alpha),
      id = id,
      design = x_design
    )
  } else {
    p_draws <- u_draws / u_draws
  }

  log_y_given_psi <- matrix(0, nrow = length(id), ncol = n_samples)

  if (zi) {
    log_y_given_psi[is_zero, ] <- log(1 - p_draws[is_zero, ])
  }

  log_y_given_psi[is_positive, ] <-
    log(p_draws[is_positive, ]) +
    lchoose(s[is_positive], y[is_positive]) +
    lbeta(
      y[is_positive] + phi * u_draws[is_positive, ],
      s[is_positive] - y[is_positive] + phi * (1 - u_draws[is_positive, ])
    ) -
    lbeta(
      phi * u_draws[is_positive, ],
      phi * (1 - u_draws[is_positive, ])
    )

  P1 <- rowsum(log_y_given_psi, id)

  mu_random <- matrix(rep(mu[random_index], n_samples), ncol = n_samples)
  random_draws <- array(
    psi_draws[, random_index, ],
    dim = c(nrow(psi_mean), n_random, n_samples)
  )

  P2 <- apply(random_draws, 1, function(x) {
    -0.5 * (
      diag(t(x - mu_random) %*% G_inv %*% (x - mu_random)) +
        log(G_det) +
        n_random * log(2 * pi)
    )
  })

  log_proposal_terms <-
    dt(
      array(t_draws[, random_index, ], dim = c(nrow(psi_mean), n_random, n_samples)),
      df = 5,
      log = TRUE
    ) -
    log(array(sd_array[, random_index, ], dim = c(nrow(psi_mean), n_random, n_samples)))

  P3 <- apply(log_proposal_terms, 3, rowSums)

  log_weights <- P1 + t(P2) - P3
  sum(log(rowMeans(exp(log_weights))))
}


#### Gradiente y Hessiano de la log-verosimilitud completa ####

.zibbmr_complete_grad <- function(mu, G, phi, zi, psi_chain,
                                  random_index, alpha_random, beta_random,
                                  n_random, x_design_chain = NULL, id_chain,
                                  is_positive_chain, is_zero_chain,
                                  n_alpha, n_beta, z_design_chain,
                                  y_chain, s_chain,
                                  n_alpha_random, n_beta_random) {
  G_diag <- .saem_diag(G)
  n_rows <- nrow(psi_chain)

  if (zi) {
    alpha <- mu[seq_len(n_alpha)]
  }
  beta <- mu[n_alpha + seq_len(n_beta)]

  psi_sum <- colSums(as.matrix(psi_chain[, random_index, drop = FALSE]))
  psi_sum2 <- colSums(as.matrix(psi_chain[, random_index, drop = FALSE]^2))

  mu_grad <- rep(0, length(mu))
  mu_grad[random_index] <- (psi_sum - n_rows * mu[random_index]) / G_diag

  G_grad <- 0.5 * (
    (psi_sum2 - 2 * mu[random_index] * psi_sum + n_rows * mu[random_index]^2) /
      G_diag^2 -
      n_rows / G_diag
  )

  alpha_grad <- NULL

  if (zi && n_alpha_random != n_alpha) {
    alpha_grad <- -numDeriv::grad(
      .saem_neg_loglik_zero,
      alpha[!alpha_random],
      psi_chain = psi_chain,
      alpha_random = alpha_random,
      x_design_chain = x_design_chain,
      id_chain = id_chain,
      is_positive_chain = is_positive_chain,
      is_zero_chain = is_zero_chain,
      n_alpha = n_alpha
    )
  }

  if (n_beta_random != n_beta) {
    beta_phi_par <- c(beta[!beta_random], phi)
  } else {
    beta_phi_par <- phi
  }

  beta_phi_grad <- -numDeriv::grad(
    .zibbmr_neg_loglik_beta_binomial,
    beta_phi_par,
    psi_chain = psi_chain,
    beta_random = beta_random,
    z_design_chain = z_design_chain,
    id_chain = id_chain,
    is_positive_chain = is_positive_chain,
    y_chain = y_chain,
    s_chain = s_chain,
    n_alpha = n_alpha,
    n_beta = n_beta,
    n_beta_random = n_beta_random
  )

  phi_grad <- beta_phi_grad[length(beta_phi_par)]

  if (n_beta_random != n_beta) {
    beta_grad <- beta_phi_grad[-length(beta_phi_par)]
  } else {
    beta_grad <- NULL
  }

  if (n_random != length(mu)) {
    mu_grad[-random_index] <- c(alpha_grad, beta_grad)
  }

  c(mu_grad, phi_grad, G_grad)
}

.zibbmr_complete_hess <- function(mu, G, phi, zi, psi_chain,
                                  random_index, alpha_random, beta_random,
                                  n_random, x_design_chain = NULL, id_chain,
                                  is_positive_chain, is_zero_chain,
                                  n_alpha, n_beta, z_design_chain,
                                  y_chain, s_chain,
                                  n_alpha_random, n_beta_random) {
  G_diag <- .saem_diag(G)
  n_rows <- nrow(psi_chain)

  if (zi) {
    alpha <- mu[seq_len(n_alpha)]
  }
  beta <- mu[n_alpha + seq_len(n_beta)]

  psi_sum <- colSums(as.matrix(psi_chain[, random_index, drop = FALSE]))
  psi_sum2 <- colSums(as.matrix(psi_chain[, random_index, drop = FALSE]^2))

  n_hess <- n_alpha + n_beta + n_random + 1
  H <- matrix(0, nrow = n_hess, ncol = n_hess)

  diag(H)[random_index] <- -n_rows / G_diag

  variance_index <- (n_hess - n_random + 1):n_hess
  diag(H)[variance_index] <-
    (1 / G_diag^2) *
    (
      0.5 * n_rows -
        (psi_sum2 - 2 * mu[random_index] * psi_sum + n_rows * mu[random_index]^2) /
        G_diag
    )

  H[random_index, variance_index] <-
    .saem_diag(-1 / G_diag^2) %*% .saem_diag(psi_sum - n_rows * mu[random_index])
  H[variance_index, random_index] <-
    .saem_diag(-1 / G_diag^2) %*% .saem_diag(psi_sum - n_rows * mu[random_index])

  if (zi && n_alpha_random != n_alpha) {
    alpha_hess <- -numDeriv::hessian(
      .saem_neg_loglik_zero,
      alpha[!alpha_random],
      psi_chain = psi_chain,
      alpha_random = alpha_random,
      x_design_chain = x_design_chain,
      id_chain = id_chain,
      is_positive_chain = is_positive_chain,
      is_zero_chain = is_zero_chain,
      n_alpha = n_alpha
    )

    alpha_fixed_index <- setdiff(seq_len(n_alpha), which(alpha_random))
    H[alpha_fixed_index, alpha_fixed_index] <- alpha_hess
  }

  if (n_beta_random != n_beta) {
    beta_phi_par <- c(beta[!beta_random], phi)
  } else {
    beta_phi_par <- phi
  }

  beta_phi_hess <- -numDeriv::hessian(
    .zibbmr_neg_loglik_beta_binomial,
    beta_phi_par,
    psi_chain = psi_chain,
    beta_random = beta_random,
    z_design_chain = z_design_chain,
    id_chain = id_chain,
    is_positive_chain = is_positive_chain,
    y_chain = y_chain,
    s_chain = s_chain,
    n_alpha = n_alpha,
    n_beta = n_beta,
    n_beta_random = n_beta_random
  )

  beta_phi_index <- setdiff(seq_len(n_beta + 1), which(beta_random)) + n_alpha
  H[beta_phi_index, beta_phi_index] <- beta_phi_hess

  H
}


#### Ajuste principal SAEM-ZIBBMR ####

#' Ajustar un modelo ZIBBMR (zero-inflated beta-binomial mixed regression) via SAEM
#'
#' Estima por Stochastic Approximation EM (SAEM) un modelo mixto
#' beta-binomial con inflacion de ceros para datos longitudinales de conteo
#' con profundidad de secuenciacion conocida (`S`), siguiendo el metodo
#' descrito en Barrera (ZIBBMR: "Stochastic EM Estimation and Inference in
#' Zero-Inflated Beta-Binomial Mixed Models for Longitudinal Count Data"). Es
#' el analogo de [fit_zibr()] para conteos: en vez de modelar una proporcion
#' observada directamente, modela el numero de lecturas `y` de un taxon sobre
#' un total `S` (profundidad de secuenciacion de la muestra) con una
#' verosimilitud beta-binomial.
#'
#' @param y Vector de conteos (numero de lecturas del taxon), `0 <= y <= S`.
#' @param S Vector con el total de lecturas (profundidad de secuenciacion) de
#'   cada observacion, misma longitud que `y`.
#' @param id Vector (o factor) que identifica al sujeto de cada observacion.
#' @param X Matriz o data frame de covariables para la parte de inflacion de
#'   ceros (parte logistica). `NULL` si `zi = FALSE`.
#' @param Z Matriz o data frame de covariables para la parte beta-binomial
#'   (magnitud condicional). `NULL` equivale a solo intercepto.
#' @param zi Logico. Si `TRUE` (por defecto) ajusta la parte de inflacion de
#'   ceros.
#' @param phi_start Valor inicial del parametro de dispersion `phi`.
#' @param alpha_start Vector de valores iniciales para los coeficientes de la
#'   parte logistica. Requerido si `zi = TRUE`.
#' @param beta_start Vector de valores iniciales para los coeficientes de la
#'   parte beta-binomial.
#' @param n_iter Numero de iteraciones del algoritmo SAEM.
#' @param n_chains Numero de cadenas MCMC paralelas usadas en el S-step.
#' @param seed Semilla aleatoria opcional.
#' @param alpha_random Vector logico que indica que coeficientes de la parte
#'   logistica son efectos aleatorios (por defecto, solo el intercepto).
#' @param beta_random Vector logico que indica que coeficientes de la parte
#'   beta-binomial son efectos aleatorios (por defecto, solo el intercepto).
#' @param n_is Numero de muestras de importance sampling para la
#'   log-verosimilitud marginal.
#' @param compute_fim Logico. Si `TRUE`, calcula la matriz de informacion de
#'   Fisher estocastica (necesaria para `vcov()`/`se()`).
#' @param cov_random Estructura de la matriz de covarianza de los efectos
#'   aleatorios. `"diag"` (por defecto) los trata como independientes entre la
#'   parte de inflacion de ceros y la parte beta-binomial: es la especificacion
#'   del articulo y la unica que las alternativas (`glmmTMB`, `gamlss`) tambien
#'   pueden ajustar. `"unstructured"` estima ademas la covarianza entre ambos.
#'
#'   Conviene usar `"unstructured"` solo con bastante informacion: con pocos
#'   sujetos u observaciones por sujeto la matriz sin restricciones degenera
#'   (la correlacion se va a \eqn{\pm 1} y una de las varianzas colapsa). Como
#'   referencia, con 40 sujetos y 4 observaciones cada uno se ha observado
#'   \eqn{\hat\rho = -0.92} sobre datos generados con \eqn{\rho = 0}.
#'
#'   Ojo: la matriz de informacion de Fisher esta derivada para el caso
#'   diagonal, asi que con `"unstructured"` los errores estandar son
#'   aproximados y no hay error estandar para la correlacion. Para contrastar
#'   \eqn{H_0:\rho=0} conviene usar el test de razon de verosimilitud sobre
#'   `logLik()` en vez de un test de Wald.
#'
#' @return Un objeto de clase `zibbmr_saem` (y `SAEM_ZIBBMR_result` por
#'   compatibilidad), con los mismos elementos que [fit_zibr()] (`mu`, `G`,
#'   `phi`, `loglik`, `trace`, `fisher_stoch`, etc.). Tiene metodos
#'   [print()], [plot()], [stats::logLik()], [stats::coef()],
#'   [stats::vcov()] y [se()].
#'
#' @seealso [fit_zibbmr_taxon()], [simulate_zibbmr_data()], [lrt_zibbmr()],
#'   [fit_zibr()] para la version en proporciones.
#'
#' @examples
#' \donttest{
#' n_subjects <- 20
#' n_time <- 4
#' n_obs <- n_subjects * n_time
#' dat <- simulate_zibbmr_data(
#'   n_subjects = n_subjects, n_time = n_time, S = rep(1000, n_obs),
#'   alpha = c(-0.3, 0.5), beta = c(0.2, -0.4),
#'   sigma_alpha = 0.4, sigma_beta = 0.3, phi = 15,
#'   X = matrix(rbinom(n_obs, 1, 0.5)), Z = matrix(rbinom(n_obs, 1, 0.5)),
#'   seed = 1
#' )
#' fit <- fit_zibbmr(
#'   y = dat$Y, S = dat$TotalCounts, id = dat$Subject,
#'   X = dat$X.1, Z = dat$Z.1,
#'   phi_start = 10, alpha_start = c(-0.2, 0.1), beta_start = c(0.1, 0.1),
#'   n_iter = 50, seed = 1, compute_fim = FALSE
#' )
#' print(fit)
#' }
#' @export
fit_zibbmr <- function(y, S, id, X = NULL, Z = NULL, zi = TRUE,
                       phi_start, alpha_start = NULL, beta_start,
                       n_iter = 1000, n_chains = 5, seed = NULL,
                       alpha_random = NULL, beta_random = NULL,
                       n_is = 500, compute_fim = TRUE,
                       cov_random = c("diag", "unstructured")) {
  .saem_check_packages(inference = compute_fim)
  cov_random <- .saem_valida_estructura(cov_random)

  if (!is.null(seed)) {
    set.seed(seed)
  }

  y <- as.numeric(y)
  S <- as.numeric(S)

  if (length(y) != length(S) || length(y) != length(id)) {
    stop("y, S e id deben tener la misma longitud.", call. = FALSE)
  }

  if (any(y < 0) || any(S < 0) || any(y > S)) {
    stop("Los conteos deben satisfacer 0 <= y <= S.", call. = FALSE)
  }

  n_total <- length(y)
  subject_id <- as.numeric(factor(id, levels = unique(id)))
  n_subjects <- length(unique(subject_id))

  if (zi) {
    x_design <- .saem_covariate_matrix(X, n_total, "X")
    n_alpha <- ncol(x_design)

    if (length(alpha_start) != n_alpha) {
      stop("alpha_start debe tener longitud igual a 1 + numero de columnas de X.", call. = FALSE)
    }

    alpha_random <- if (is.null(alpha_random)) c(TRUE, rep(FALSE, n_alpha - 1)) else alpha_random
    if (length(alpha_random) != n_alpha) {
      stop("alpha_random debe tener longitud igual a length(alpha_start).", call. = FALSE)
    }

    x_design_chain <- .saem_replicate_design(x_design, n_chains)
    alpha_labels <- colnames(x_design)
    n_alpha_random <- sum(alpha_random)
  } else {
    if (!is.null(X) || !is.null(alpha_start)) {
      stop("No entregue X ni alpha_start cuando zi = FALSE.", call. = FALSE)
    }

    y <- ifelse(y == 0, 1e-6, y)
    x_design <- NULL
    x_design_chain <- NULL
    n_alpha <- 0
    alpha_start <- NULL
    alpha_random <- NULL
    alpha_labels <- NULL
    n_alpha_random <- 0
  }

  z_design <- .saem_covariate_matrix(Z, n_total, "Z")
  n_beta <- ncol(z_design)

  if (length(beta_start) != n_beta) {
    stop("beta_start debe tener longitud igual a 1 + numero de columnas de Z.", call. = FALSE)
  }

  beta_random <- if (is.null(beta_random)) c(TRUE, rep(FALSE, n_beta - 1)) else beta_random
  if (length(beta_random) != n_beta) {
    stop("beta_random debe tener longitud igual a length(beta_start).", call. = FALSE)
  }

  beta_labels <- colnames(z_design)
  n_beta_random <- sum(beta_random)
  z_design_chain <- .saem_replicate_design(z_design, n_chains)

  is_positive <- y != 0
  is_zero <- y == 0

  n_psi <- length(c(alpha_start, beta_start))
  random_index <- which(c(alpha_random, beta_random))
  n_random <- length(random_index)

  id_rep <- rep(subject_id, n_chains)
  id_chain <- id_rep + n_subjects * (rep(seq_len(n_chains), each = n_total) - 1)

  # Grupo sujeto para colapsar las filas (sujeto, cadena) de psi_chain sobre
  # las cadenas; cada sujeto aparece exactamente n_chains veces, asi que la
  # media por sujeto es rowsum(.)/n_chains (equivalente a tapply(., mean) pero
  # sin reconstruir el factor en cada iteracion).
  subject_group_chain <- rep(seq_len(n_subjects), n_chains)

  is_positive_chain <- rep(is_positive, n_chains)
  is_zero_chain <- rep(is_zero, n_chains)
  y_chain <- rep(y, n_chains)
  S_chain <- rep(S, n_chains)

  mu <- c(alpha_start, beta_start)
  G_full <- 0.5 * .saem_diag(abs(mu))
  G <- .saem_impone_estructura(G_full[random_index, random_index, drop = FALSE], cov_random)
  phi <- phi_start

  psi_chain <- matrix(
    rep(mu, n_chains * n_subjects),
    nrow = n_chains * n_subjects,
    ncol = n_psi,
    byrow = TRUE
  )

  u_chain <- .saem_linear_prob(
    psi_chain,
    n_alpha + seq_len(n_beta),
    id_chain,
    z_design_chain
  )

  if (zi) {
    p_chain <- .saem_linear_prob(psi_chain, seq_len(n_alpha), id_chain, x_design_chain)
  } else {
    p_chain <- u_chain / u_chain
  }

  proposal_sd_uni <- 0.5 * G
  proposal_sd_uni[proposal_sd_uni < 0.5 & proposal_sd_uni > 0] <- 0.5
  proposal_sd_multi <- proposal_sd_uni

  if (compute_fim) {
    grad_avg <- rep(0, n_psi + n_random + 1)
    hess_avg <- matrix(0, nrow = n_psi + n_random + 1, ncol = n_psi + n_random + 1)
    score2_avg <- matrix(0, nrow = n_psi + n_random + 1, ncol = n_psi + n_random + 1)
    fisher_stoch <- matrix(0, nrow = n_psi + n_random + 1, ncol = n_psi + n_random + 1)
  } else {
    fisher_stoch <- NULL
  }

  psi_subject_mean <- rowsum(psi_chain, subject_group_chain) / n_chains
  psi_subject_second <- rowsum(psi_chain^2, subject_group_chain) / n_chains

  stat1 <- n_subjects * mu
  stat2 <- G_full

  trace <- NULL
  burn_in <- floor(0.75 * n_iter)

  for (iter in seq_len(n_iter)) {
    gamma <- if (iter <= burn_in) 1 else 1 / (iter - burn_in)

    mu_chain <- matrix(
      rep(mu, n_chains * n_subjects),
      nrow = n_chains * n_subjects,
      ncol = n_psi,
      byrow = TRUE
    )

    G_inv <- .saem_diag_inverse(G)
    log_ratio_data <- rep(0, n_chains * n_total)

    for (mh_iter in seq_len(4)) {
      psi_candidate <- MASS::mvrnorm(n_subjects * n_chains, mu, G_full)
      psi_candidate[, -random_index] <- psi_chain[, -random_index]

      if (zi) {
        p_candidate <- .saem_linear_prob(
          psi_candidate,
          seq_len(n_alpha),
          id_chain,
          x_design_chain
        )
        log_ratio_data[is_zero_chain] <-
          log(1 - p_chain[is_zero_chain]) - log(1 - p_candidate[is_zero_chain])
      } else {
        p_candidate <- p_chain
      }

      u_candidate <- .saem_linear_prob(
        psi_candidate,
        n_alpha + seq_len(n_beta),
        id_chain,
        z_design_chain
      )

      log_ratio_data[is_positive_chain] <-
        log(p_chain[is_positive_chain]) - log(p_candidate[is_positive_chain]) +
        lgamma(y_chain[is_positive_chain] + phi * u_chain[is_positive_chain]) +
        lgamma(S_chain[is_positive_chain] - y_chain[is_positive_chain] +
                 phi * (1 - u_chain[is_positive_chain])) -
        lgamma(y_chain[is_positive_chain] + phi * u_candidate[is_positive_chain]) -
        lgamma(S_chain[is_positive_chain] - y_chain[is_positive_chain] +
                 phi * (1 - u_candidate[is_positive_chain])) +
        lgamma(phi * u_candidate[is_positive_chain]) +
        lgamma(phi * (1 - u_candidate[is_positive_chain])) -
        lgamma(phi * u_chain[is_positive_chain]) -
        lgamma(phi * (1 - u_chain[is_positive_chain]))

      subject_log_ratio <- as.vector(rowsum(log_ratio_data, id_chain))
      accept <- subject_log_ratio < -log(runif(n_subjects * n_chains))

      psi_chain <- psi_candidate * accept + psi_chain * (!accept)

      if (zi) {
        p_chain <- .saem_linear_prob(psi_chain, seq_len(n_alpha), id_chain, x_design_chain)
      }
      u_chain <- .saem_linear_prob(psi_chain, n_alpha + seq_len(n_beta), id_chain, z_design_chain)
    }

    accepted_uni <- proposed_uni <- rep(0, n_random)

    for (mh_iter in seq_len(4)) {
      delta <- matrix(0, nrow = n_subjects * n_chains, ncol = n_random)
      delta[
        matrix(
          c(
            seq_len(n_chains * n_subjects),
            sample(seq_len(n_random), n_chains * n_subjects, replace = TRUE)
          ),
          nrow = n_chains * n_subjects
        )
      ] <- rnorm(n_chains * n_subjects)

      psi_candidate <- psi_chain
      psi_candidate[, random_index] <- psi_chain[, random_index] + delta %*% proposal_sd_uni

      if (zi) {
        p_candidate <- .saem_linear_prob(
          psi_candidate,
          seq_len(n_alpha),
          id_chain,
          x_design_chain
        )
        log_ratio_data[is_zero_chain] <-
          log(1 - p_chain[is_zero_chain]) - log(1 - p_candidate[is_zero_chain])
      } else {
        p_candidate <- p_chain
      }

      u_candidate <- .saem_linear_prob(
        psi_candidate,
        n_alpha + seq_len(n_beta),
        id_chain,
        z_design_chain
      )

      log_ratio_data[is_positive_chain] <-
        log(p_chain[is_positive_chain]) - log(p_candidate[is_positive_chain]) +
        lgamma(y_chain[is_positive_chain] + phi * u_chain[is_positive_chain]) +
        lgamma(S_chain[is_positive_chain] - y_chain[is_positive_chain] +
                 phi * (1 - u_chain[is_positive_chain])) -
        lgamma(y_chain[is_positive_chain] + phi * u_candidate[is_positive_chain]) -
        lgamma(S_chain[is_positive_chain] - y_chain[is_positive_chain] +
                 phi * (1 - u_candidate[is_positive_chain])) +
        lgamma(phi * u_candidate[is_positive_chain]) +
        lgamma(phi * (1 - u_candidate[is_positive_chain])) -
        lgamma(phi * u_chain[is_positive_chain]) -
        lgamma(phi * (1 - u_chain[is_positive_chain]))

      d_candidate <- psi_candidate[, random_index, drop = FALSE] - mu_chain[, random_index, drop = FALSE]
      d_current <- psi_chain[, random_index, drop = FALSE] - mu_chain[, random_index, drop = FALSE]

      subject_log_ratio <- as.vector(rowsum(log_ratio_data, id_chain)) +
        0.5 * (
          rowSums((d_candidate %*% G_inv) * d_candidate) -
            rowSums((d_current %*% G_inv) * d_current)
        )

      accept <- subject_log_ratio < -log(runif(n_subjects * n_chains))
      psi_chain <- psi_candidate * accept + psi_chain * (!accept)

      if (zi) {
        p_chain <- .saem_linear_prob(psi_chain, seq_len(n_alpha), id_chain, x_design_chain)
      }
      u_chain <- .saem_linear_prob(psi_chain, n_alpha + seq_len(n_beta), id_chain, z_design_chain)

      accepted_uni <- accepted_uni + colSums((delta * accept) != 0)
      proposed_uni <- proposed_uni + colSums(delta != 0)
    }

    proposal_sd_uni <- (1 + 0.5 * (accepted_uni / proposed_uni - 0.4)) * proposal_sd_uni

    accepted_multi <- 0

    for (mh_iter in seq_len(4)) {
      psi_candidate <- psi_chain
      psi_candidate[, random_index] <-
        psi_chain[, random_index] +
        matrix(rnorm(n_chains * n_subjects * n_random),
               nrow = n_chains * n_subjects,
               ncol = n_random) %*%
        proposal_sd_multi

      if (zi) {
        p_candidate <- .saem_linear_prob(
          psi_candidate,
          seq_len(n_alpha),
          id_chain,
          x_design_chain
        )
        log_ratio_data[is_zero_chain] <-
          log(1 - p_chain[is_zero_chain]) - log(1 - p_candidate[is_zero_chain])
      } else {
        p_candidate <- p_chain
      }

      u_candidate <- .saem_linear_prob(
        psi_candidate,
        n_alpha + seq_len(n_beta),
        id_chain,
        z_design_chain
      )

      log_ratio_data[is_positive_chain] <-
        log(p_chain[is_positive_chain]) - log(p_candidate[is_positive_chain]) +
        lgamma(y_chain[is_positive_chain] + phi * u_chain[is_positive_chain]) +
        lgamma(S_chain[is_positive_chain] - y_chain[is_positive_chain] +
                 phi * (1 - u_chain[is_positive_chain])) -
        lgamma(y_chain[is_positive_chain] + phi * u_candidate[is_positive_chain]) -
        lgamma(S_chain[is_positive_chain] - y_chain[is_positive_chain] +
                 phi * (1 - u_candidate[is_positive_chain])) +
        lgamma(phi * u_candidate[is_positive_chain]) +
        lgamma(phi * (1 - u_candidate[is_positive_chain])) -
        lgamma(phi * u_chain[is_positive_chain]) -
        lgamma(phi * (1 - u_chain[is_positive_chain]))

      d_candidate <- psi_candidate[, random_index, drop = FALSE] - mu_chain[, random_index, drop = FALSE]
      d_current <- psi_chain[, random_index, drop = FALSE] - mu_chain[, random_index, drop = FALSE]

      subject_log_ratio <- as.vector(rowsum(log_ratio_data, id_chain)) +
        0.5 * (
          rowSums((d_candidate %*% G_inv) * d_candidate) -
            rowSums((d_current %*% G_inv) * d_current)
        )

      accept <- subject_log_ratio < -log(runif(n_subjects * n_chains))
      psi_chain <- psi_candidate * accept + psi_chain * (!accept)

      if (zi) {
        p_chain <- .saem_linear_prob(psi_chain, seq_len(n_alpha), id_chain, x_design_chain)
      }
      u_chain <- .saem_linear_prob(psi_chain, n_alpha + seq_len(n_beta), id_chain, z_design_chain)

      accepted_multi <- accepted_multi + sum(accept)
    }

    proposal_sd_multi <-
      (1 + 0.5 * (accepted_multi / (4 * n_chains * n_subjects) - 0.4)) *
      proposal_sd_multi

    psi_subject_mean <- psi_subject_mean +
      gamma * (
        rowsum(psi_chain, subject_group_chain) / n_chains -
          psi_subject_mean
      )

    psi_subject_second <- psi_subject_second +
      gamma * (
        rowsum(psi_chain^2, subject_group_chain) / n_chains -
          psi_subject_second
      )

    stat1 <- stat1 + gamma * (colSums(psi_chain) / n_chains - stat1)
    stat2 <- stat2 + gamma * ((t(psi_chain) %*% psi_chain) / n_chains - stat2)

    if (iter > 10) {
      mu <- stat1 / n_subjects
      G_full <- stat2 / n_subjects - (stat1 %*% t(stat1)) / n_subjects^2
      G <- .saem_impone_estructura(G_full[random_index, random_index, drop = FALSE], cov_random)
      # El kernel de propuesta usa G_full, asi que tiene que respetar la misma
      # restriccion; si no, propone en direcciones que el modelo no admite.
      G_full[random_index, random_index] <- G

      beta <- mu[n_alpha + seq_len(n_beta)]

      if (zi) {
        alpha <- mu[seq_len(n_alpha)]

        if (n_alpha_random != n_alpha) {
          alpha_opt <- stats::nlminb(
            start = alpha[!alpha_random],
            objective = .saem_neg_loglik_zero,
            psi_chain = psi_chain,
            alpha_random = alpha_random,
            x_design_chain = x_design_chain,
            id_chain = id_chain,
            is_positive_chain = is_positive_chain,
            is_zero_chain = is_zero_chain,
            n_alpha = n_alpha
          )$par

          alpha[!alpha_random] <- alpha[!alpha_random] +
            gamma * (alpha_opt - alpha[!alpha_random])
        }
      } else {
        alpha <- NULL
      }

      if (n_beta_random != n_beta) {
        beta_phi_par <- c(beta[!beta_random], phi)
      } else {
        beta_phi_par <- phi
      }

      beta_phi_opt <- stats::nlminb(
        start = beta_phi_par,
        objective = .zibbmr_neg_loglik_beta_binomial,
        psi_chain = psi_chain,
        beta_random = beta_random,
        z_design_chain = z_design_chain,
        id_chain = id_chain,
        is_positive_chain = is_positive_chain,
        y_chain = y_chain,
        s_chain = S_chain,
        n_alpha = n_alpha,
        n_beta = n_beta,
        n_beta_random = n_beta_random,
        lower = c(rep(-Inf, length(beta_phi_par) - 1), 0.0001)
      )$par

      beta_phi_par <- beta_phi_par + gamma * (beta_phi_opt - beta_phi_par)

      phi <- beta_phi_par[length(beta_phi_par)]

      if (n_beta_random != n_beta) {
        beta[!beta_random] <- beta_phi_par[-length(beta_phi_par)]
      }

      mu <- c(alpha, beta)
      G_full[, -random_index] <- 0
      G_full[-random_index, ] <- 0

      psi_chain[, -random_index] <- matrix(
        rep(mu[-random_index], each = n_chains * n_subjects),
        ncol = n_psi - n_random,
        nrow = n_chains * n_subjects
      )

      if (compute_fim) {
        grad_current <- .zibbmr_complete_grad(
          mu, G, phi, zi, psi_chain, random_index, alpha_random,
          beta_random, n_random, x_design_chain, id_chain,
          is_positive_chain, is_zero_chain, n_alpha, n_beta,
          z_design_chain, y_chain, S_chain, n_alpha_random, n_beta_random
        )

        hess_current <- .zibbmr_complete_hess(
          mu, G, phi, zi, psi_chain, random_index, alpha_random,
          beta_random, n_random, x_design_chain, id_chain,
          is_positive_chain, is_zero_chain, n_alpha, n_beta,
          z_design_chain, y_chain, S_chain, n_alpha_random, n_beta_random
        )

        score2_current <- matrix(0, nrow = n_psi + n_random + 1, ncol = n_psi + n_random + 1)

        for (chain in seq_len(n_chains)) {
          row_index <- n_subjects * (chain - 1) + seq_len(n_subjects)

          grad_chain <- .zibbmr_complete_grad(
            mu, G, phi, zi, psi_chain[row_index, , drop = FALSE],
            random_index, alpha_random, beta_random, n_random,
            x_design, subject_id, is_positive, is_zero,
            n_alpha, n_beta, z_design, y, S,
            n_alpha_random, n_beta_random
          )

          score2_current <- score2_current + grad_chain %*% t(grad_chain)
        }

        grad_avg <- grad_avg + gamma * (grad_current - grad_avg)
        hess_avg <- hess_avg + gamma * (hess_current - hess_avg)
        score2_avg <- score2_avg + gamma * (score2_current - score2_avg)

        fisher_stoch <-
          (hess_avg + score2_avg) / n_chains -
          (grad_avg %*% t(grad_avg)) / n_chains^2
      }
    }

    trace <- rbind(trace, c(mu, diag(G), phi))

    if (zi) {
      colnames(trace) <- c(
        paste("alpha", seq_len(n_alpha), sep = "."),
        paste("beta", seq_len(n_beta), sep = "."),
        paste("var", seq_len(n_random), sep = "."),
        "phi"
      )
    } else {
      colnames(trace) <- c(
        paste("beta", seq_len(n_beta), sep = "."),
        paste("var", seq_len(n_random), sep = "."),
        "phi"
      )
    }
  }

  psi_mean <- rowsum(psi_chain, subject_group_chain) / n_chains
  psi_var <- psi_subject_second - psi_subject_mean^2
  psi_var[, -random_index] <- 0

  loglik <- .zibbmr_loglik_importance(
    mu = mu,
    G = G,
    phi = phi,
    zi = zi,
    y = y,
    s = S,
    id = subject_id,
    x_design = x_design,
    z_design = z_design,
    n_alpha = n_alpha,
    n_beta = n_beta,
    is_positive = is_positive,
    is_zero = is_zero,
    psi_mean = psi_mean,
    psi_var = psi_var,
    random_index = random_index,
    n_random = n_random,
    n_samples = n_is,
    seed = seed
  )

  out <- list(
    mu = mu,
    G = G,
    cov_random = cov_random,
    phi = phi,
    psi_mean = psi_mean,
    psi_var = psi_var,
    loglik = loglik,
    zi = zi,
    trace = trace,
    alpha_labels = alpha_labels,
    beta_labels = beta_labels,
    n_alpha = n_alpha,
    n_beta = n_beta,
    alpha_random = alpha_random,
    beta_random = beta_random,
    random_index = random_index,
    fisher_stoch = fisher_stoch,
    nobs = n_total,
    # Datos originales, para poder calcular predicciones y residuos en plot().
    # No afectan la estimacion; solo se guardan para los graficos. Incluye S
    # (total de lecturas por muestra), que ZIBBMR necesita para las predicciones.
    data = list(y = y, S = S, x_design = x_design, z_design = z_design,
                subject_id = subject_id),
    call = match.call()
  )

  out$MU <- out$mu
  out$V <- out$phi
  out$psi.mean <- out$psi_mean
  out$psi.var <- out$psi_var
  out$graph <- out$trace
  out$labs.X <- out$alpha_labels
  out$labs.Z <- out$beta_labels
  out$nxcov <- out$n_alpha
  out$nzcov <- out$n_beta
  out$ind.a.aleat <- out$alpha_random
  out$ind.b.aleat <- out$beta_random
  out["FIM.stoch"] <- list(out$fisher_stoch)

  class(out) <- c("zibbmr_saem", "SAEM_ZIBBMR_result")
  out
}


#### Simulacion de datos ZIBBMR ####

#' Simular datos longitudinales de conteo para un modelo ZIBBMR
#'
#' Genera un data frame de datos longitudinales compatibles con
#' [fit_zibbmr()]: un conteo de lecturas `Y` sobre un total `S`, con
#' inflacion de ceros opcional, efectos fijos y un intercepto (u otros
#' coeficientes) aleatorios por sujeto.
#'
#' @param n_subjects Numero de sujetos.
#' @param n_time Numero de observaciones (tiempos) por sujeto.
#' @param S Vector con el total de lecturas de cada observacion, de longitud
#'   `n_subjects * n_time`.
#' @param zi Logico. Si `TRUE` (por defecto), simula presencia/ausencia con
#'   un modelo logistico (`X`/`alpha`) antes de simular el conteo.
#' @param X Matriz o data frame de covariables para la parte de inflacion de
#'   ceros. `NULL` si `zi = FALSE`.
#' @param Z Matriz o data frame de covariables para la parte beta-binomial.
#' @param alpha Vector de coeficientes verdaderos de la parte logistica.
#'   Requerido si `zi = TRUE`.
#' @param beta Vector de coeficientes verdaderos de la parte beta-binomial.
#' @param sigma_alpha Desviacion estandar del intercepto aleatorio de la
#'   parte logistica. Requerido si `zi = TRUE`.
#' @param sigma_beta Desviacion estandar del intercepto aleatorio de la parte
#'   beta-binomial.
#' @param phi Parametro de dispersion de la distribucion beta-binomial.
#' @param seed Semilla aleatoria opcional.
#'
#' @return Un data frame con columnas `Subject`, `Time`, `Y` (conteo
#'   simulado), `TotalCounts` (igual a `S`) y las covariables usadas.
#' @seealso [fit_zibbmr()]
#' @examples
#' n_subjects <- 10
#' n_time <- 3
#' n_obs <- n_subjects * n_time
#' dat <- simulate_zibbmr_data(
#'   n_subjects = n_subjects, n_time = n_time, S = rep(1000, n_obs),
#'   alpha = c(-0.3, 0.5), beta = c(0.2, -0.4),
#'   sigma_alpha = 0.4, sigma_beta = 0.3, phi = 15,
#'   X = matrix(rbinom(n_obs, 1, 0.5)), Z = matrix(rbinom(n_obs, 1, 0.5)),
#'   seed = 1
#' )
#' head(dat)
#' @export
simulate_zibbmr_data <- function(n_subjects, n_time, S, zi = TRUE,
                                 X = NULL, Z = NULL, alpha = NULL, beta,
                                 sigma_alpha = NULL, sigma_beta,
                                 phi, seed = NULL) {
  .saem_check_packages(inference = FALSE)

  if (!is.null(seed)) {
    set.seed(seed)
  }

  n_total <- n_subjects * n_time

  if (length(S) != n_total) {
    stop("S debe tener longitud n_subjects * n_time.", call. = FALSE)
  }

  if (zi) {
    x_design <- .saem_covariate_matrix(X, n_total, "X")
    n_alpha <- ncol(x_design)

    if (length(alpha) != n_alpha) {
      stop("alpha debe tener longitud igual a 1 + numero de columnas de X.", call. = FALSE)
    }

    random_cols <- c(1, n_alpha + 1)
  } else {
    if (!is.null(X) || !is.null(alpha)) {
      stop("No entregue X ni alpha cuando zi = FALSE.", call. = FALSE)
    }

    x_design <- NULL
    n_alpha <- 0
    random_cols <- 1
  }

  z_design <- .saem_covariate_matrix(Z, n_total, "Z")
  n_beta <- ncol(z_design)

  if (length(beta) != n_beta) {
    stop("beta debe tener longitud igual a 1 + numero de columnas de Z.", call. = FALSE)
  }

  id <- rep(seq_len(n_subjects), each = n_time)
  mu <- c(alpha, beta)

  if (zi) {
    G <- .saem_diag(c(sigma_alpha^2, sigma_beta^2))
  } else {
    G <- .saem_diag(sigma_beta^2)
  }

  psi <- matrix(
    rep(mu, n_subjects),
    ncol = n_alpha + n_beta,
    nrow = n_subjects,
    byrow = TRUE
  )

  psi[, random_cols] <- MASS::mvrnorm(n = n_subjects, mu = mu[random_cols], Sigma = G)

  u <- .saem_linear_prob(psi, n_alpha + seq_len(n_beta), id, z_design)

  if (zi) {
    p <- .saem_linear_prob(psi, seq_len(n_alpha), id, x_design)
  } else {
    p <- u / u
  }

  w <- apply(phi * cbind(u, 1 - u), 1, function(par) {
    stats::rbeta(1, par[1], par[2])
  })

  y_positive <- apply(cbind(S, w), 1, function(par) {
    stats::rbinom(1, par[1], par[2])
  })

  present <- stats::rbinom(n_total, size = 1, prob = p)
  y <- present * y_positive

  if (zi) {
    covariates <- data.frame(x_design[, -1, drop = FALSE], z_design[, -1, drop = FALSE])
    colnames(covariates) <- c(colnames(x_design)[-1], colnames(z_design)[-1])
  } else {
    covariates <- data.frame(z_design[, -1, drop = FALSE])
    colnames(covariates) <- colnames(z_design)[-1]
  }

  data.frame(
    Subject = paste("Subject", id, sep = "."),
    Time = rep(seq_len(n_time), n_subjects),
    Y = y,
    TotalCounts = S,
    covariates,
    check.names = FALSE
  )
}


#### Metodos basicos ####

#' Imprimir un ajuste ZIBBMR
#'
#' @param x Un objeto `zibbmr_saem`, resultado de [fit_zibbmr()].
#' @param ... No usado, por compatibilidad con el generico [print()].
#' @return `x`, de forma invisible.
#' @export
print.zibbmr_saem <- function(x, ...) {
  .saem_print(x, model_label = "SAEM-ZIBBMR", beta_label = "Parte beta-binomial")
}

#' Graficos de un ajuste ZIBBMR
#'
#' Analogo a [plot.zibr_saem()]; genera distintos graficos segun `which`:
#' \describe{
#'   \item{`"convergencia"`}{(por defecto) traza iteracion-a-iteracion de los
#'     parametros, para revisar la convergencia del algoritmo SAEM.}
#'   \item{`"coeficientes"`}{coeficientes estimados con su intervalo de
#'     confianza al 95%. Requiere `compute_fim = TRUE` para mostrar los
#'     intervalos.}
#'   \item{`"aleatorios"`}{distribucion entre sujetos de los efectos aleatorios
#'     estimados; la linea roja marca la media poblacional.}
#'   \item{`"ajuste"`}{observados frente a predichos de la parte continua, en las
#'     observaciones positivas (el conteo esperado dado presencia, `u * S`), con
#'     la recta `y = x` de referencia.}
#'   \item{`"residuos"`}{residuos de la parte continua (observado menos predicho
#'     individual, en observaciones positivas): su dispersion contra el valor
#'     predicho y su distribucion.}
#' }
#' Los graficos `"ajuste"` y `"residuos"` usan los datos originales que el
#' ajuste guarda.
#'
#' @param x Un objeto `zibbmr_saem`, resultado de [fit_zibbmr()].
#' @param which Tipo de grafico: `"convergencia"`, `"coeficientes"`,
#'   `"aleatorios"`, `"ajuste"` o `"residuos"`.
#' @param ... Argumentos adicionales (no usados por ahora).
#' @return `x`, de forma invisible. Se llama por su efecto secundario de
#'   graficar.
#' @export
plot.zibbmr_saem <- function(x, which = c("convergencia", "coeficientes",
                                          "aleatorios", "ajuste", "residuos"), ...) {
  .saem_plot(x, which = match.arg(which), beta_label = "beta-binomial", ...)
}

#' Log-verosimilitud marginal de un ajuste ZIBBMR
#'
#' @param object Un objeto `zibbmr_saem`, resultado de [fit_zibbmr()].
#' @param ... No usado, por compatibilidad con el generico [stats::logLik()].
#' @return Un objeto `logLik` con la log-verosimilitud marginal estimada por
#'   importance sampling, con atributos `df` y `nobs`.
#' @export
logLik.zibbmr_saem <- function(object, ...) {
  .saem_logLik(object)
}

#' Coeficientes estimados de un ajuste ZIBBMR
#'
#' @param object Un objeto `zibbmr_saem`, resultado de [fit_zibbmr()].
#' @param ... No usado, por compatibilidad con el generico [stats::coef()].
#' @return Vector numerico `mu` con los coeficientes de la parte logistica
#'   seguidos de los de la parte beta-binomial.
#' @export
coef.zibbmr_saem <- function(object, ...) {
  .saem_coef(object)
}

#' Matriz de varianza-covarianza de un ajuste ZIBBMR
#'
#' @param object Un objeto `zibbmr_saem` ajustado con `compute_fim = TRUE`.
#' @param ... No usado, por compatibilidad con el generico [stats::vcov()].
#' @return Una matriz de varianza-covarianza.
#' @export
vcov.zibbmr_saem <- function(object, ...) {
  .saem_vcov(object)
}

#' @rdname se
#' @export
se.zibbmr_saem <- function(object, ...) {
  .saem_se(object)
}


#### Funciones limpias para analisis por taxon ####

#' Ajustar ZIBBMR para un taxon de un data frame
#'
#' Envoltorio de [fit_zibbmr()] pensado para trabajar directamente sobre un
#' data frame de microbioma en formato largo con conteos por taxon y una
#' columna de profundidad de secuenciacion total.
#'
#' @param data Data frame con una fila por observacion, incluyendo la columna
#'   del taxon (conteos), la columna `total`, la columna de id y las
#'   covariables.
#' @param taxon Nombre de la columna en `data` con el conteo del taxon a
#'   modelar.
#' @param covariates Vector de nombres de columnas a usar como covariables en
#'   ambas partes del modelo. Se ignora si se entregan
#'   `x_covariates`/`z_covariates` por separado.
#' @param x_covariates Nombres de columnas para la parte de inflacion de
#'   ceros (por defecto, igual a `covariates`).
#' @param z_covariates Nombres de columnas para la parte beta-binomial (por
#'   defecto, igual a `covariates`).
#' @param total Nombre de la columna en `data` con el total de lecturas de
#'   cada observacion.
#' @param id Nombre de la columna en `data` que identifica al sujeto.
#' @param zi Logico, ver [fit_zibbmr()].
#' @param phi_start Valor inicial de `phi`. Si es `NULL`, se sortea con
#'   `runif(1, 10, 20)`.
#' @param alpha_start Valores iniciales de la parte logistica. Si es `NULL`,
#'   se sortean con `runif(., -0.1, 0.1)`.
#' @param beta_start Valores iniciales de la parte beta-binomial. Si es
#'   `NULL`, se sortean con `runif(., -0.1, 0.1)`.
#' @param seed Semilla aleatoria.
#' @param n_iter Numero de iteraciones SAEM.
#' @param n_chains Numero de cadenas MCMC.
#' @param compute_fim Logico, ver [fit_zibbmr()].
#' @param ... Argumentos adicionales pasados a [fit_zibbmr()].
#'
#' @return Un objeto `zibbmr_saem`, igual que [fit_zibbmr()].
#' @seealso [fit_zibbmr()], [fit_zibbmr_taxa()]
#' @export
fit_zibbmr_taxon <- function(data, taxon, covariates = NULL,
                             x_covariates = covariates,
                             z_covariates = covariates,
                             total, id,
                             zi = TRUE, phi_start = NULL,
                             alpha_start = NULL, beta_start = NULL,
                             seed = 232, n_iter = 1000, n_chains = 5,
                             compute_fim = FALSE, ...) {
  if (!taxon %in% names(data)) {
    stop("El taxon indicado no existe en data.", call. = FALSE)
  }
  if (!total %in% names(data)) {
    stop("La columna total no existe en data.", call. = FALSE)
  }
  if (!id %in% names(data)) {
    stop("La columna id no existe en data.", call. = FALSE)
  }
  if (is.null(x_covariates)) {
    x_covariates <- character(0)
  }
  if (is.null(z_covariates)) {
    z_covariates <- character(0)
  }

  if (!all(x_covariates %in% names(data))) {
    stop("Al menos una covariable de x_covariates no existe en data.", call. = FALSE)
  }
  if (!all(z_covariates %in% names(data))) {
    stop("Al menos una covariable de z_covariates no existe en data.", call. = FALSE)
  }

  n_x_covariates <- length(x_covariates)
  n_z_covariates <- length(z_covariates)

  if (!is.null(seed)) {
    set.seed(seed)
  }

  if (is.null(phi_start)) {
    phi_start <- runif(1, 10, 20)
  }
  if (zi && is.null(alpha_start)) {
    alpha_start <- runif(n_x_covariates + 1, -0.1, 0.1)
  }
  if (is.null(beta_start)) {
    beta_start <- runif(n_z_covariates + 1, -0.1, 0.1)
  }

  X <- if (!zi || length(x_covariates) == 0) NULL else data[, x_covariates, drop = FALSE]
  Z <- if (length(z_covariates) == 0) NULL else data[, z_covariates, drop = FALSE]

  fit_zibbmr(
    y = data[[taxon]],
    S = data[[total]],
    id = data[[id]],
    X = X,
    Z = Z,
    zi = zi,
    phi_start = phi_start,
    alpha_start = alpha_start,
    beta_start = beta_start,
    n_iter = n_iter,
    n_chains = n_chains,
    seed = seed,
    compute_fim = compute_fim,
    ...
  )
}

#' Ajustar ZIBBMR para varios taxones de un data frame
#'
#' Aplica [fit_zibbmr_taxon()] a cada elemento de `taxa`, con la misma
#' configuracion de covariables e iteraciones para todos.
#'
#' @param data Data frame con una fila por observacion.
#' @param taxa Vector de nombres de columnas (taxones) a ajustar.
#' @param covariates,x_covariates,z_covariates Ver [fit_zibbmr_taxon()].
#' @param total Nombre de la columna con el total de lecturas.
#' @param id Nombre de la columna que identifica al sujeto.
#' @param zi Logico, ver [fit_zibbmr()].
#' @param seed Semilla aleatoria (se reutiliza para cada taxon).
#' @param n_iter Numero de iteraciones SAEM.
#' @param n_chains Numero de cadenas MCMC.
#' @param compute_fim Logico, ver [fit_zibbmr()].
#' @param ... Argumentos adicionales pasados a [fit_zibbmr_taxon()].
#'
#' @return Una lista de objetos `zibbmr_saem`, nombrada segun `taxa`.
#' @seealso [fit_zibbmr_taxon()]
#' @export
fit_zibbmr_taxa <- function(data, taxa, covariates = NULL,
                            x_covariates = covariates,
                            z_covariates = covariates,
                            total, id,
                            zi = TRUE, seed = 232, n_iter = 1000,
                            n_chains = 5, compute_fim = FALSE, ...) {
  fits <- lapply(taxa, function(taxon) {
    fit_zibbmr_taxon(
      data = data,
      taxon = taxon,
      x_covariates = x_covariates,
      z_covariates = z_covariates,
      total = total,
      id = id,
      zi = zi,
      seed = seed,
      n_iter = n_iter,
      n_chains = n_chains,
      compute_fim = compute_fim,
      ...
    )
  })

  names(fits) <- taxa
  fits
}

#' Prueba de razon de verosimilitudes entre dos ajustes ZIBBMR anidados
#'
#' Analogo de [lrt_zibr()] para modelos ZIBBMR.
#'
#' @param full Modelo completo: un objeto `zibbmr_saem` o cualquier objeto
#'   con metodo [stats::logLik()].
#' @param reduced Modelo reducido (anidado en `full`), mismo tipo que `full`.
#' @param df Grados de libertad de la prueba.
#'
#' @return Un data frame de una fila con `LL_full`, `LL_reduced`, `LRT`, `df`
#'   y `p_value`.
#' @seealso [lrt_zibbmr_table()], [lrt_zibr()]
#' @export
lrt_zibbmr <- function(full, reduced, df = 2) {
  .saem_lrt(full, reduced, df = df)
}

#' Tabla de pruebas de razon de verosimilitudes para varios taxones ZIBBMR
#'
#' Analogo de [lrt_zibr_table()] para modelos ZIBBMR.
#'
#' @param full_models Lista de modelos completos (uno por taxon).
#' @param reduced_models Lista de modelos reducidos (uno por taxon, mismo
#'   orden que `full_models`).
#' @param species Vector de nombres/etiquetas para cada taxon.
#' @param df Grados de libertad de la prueba.
#' @param alpha Nivel de significancia usado para marcar `Detected`.
#'
#' @return Un data frame con una fila por taxon, ver [lrt_zibbmr()].
#' @seealso [lrt_zibbmr()]
#' @export
lrt_zibbmr_table <- function(full_models, reduced_models, species = names(full_models),
                             df = 2, alpha = 0.05) {
  .saem_lrt_table(
    full_models, reduced_models,
    species = species, df = df, alpha = alpha, lrt_fn = lrt_zibbmr
  )
}

#' Tabla resumen de tres comparaciones LRT tipicas para ZIBBMR
#'
#' Analogo de [zibr_results_table()] para modelos ZIBBMR.
#'
#' @param species Vector de nombres/etiquetas para cada taxon.
#' @param mod1_full,mod1_no_preg Listas de modelos ZIBBMR para la primera
#'   comparacion, uno por taxon.
#' @param mod2_full,mod2_no_preg,mod2_no_inter Listas de modelos ZIBBMR para
#'   la segunda comparacion y la prueba de interaccion, uno por taxon.
#' @param df Grados de libertad usados en las tres pruebas.
#' @param alpha Nivel de significancia usado para las columnas `Detec_*`.
#'
#' @return Un data frame con una fila por taxon.
#' @seealso [zibr_results_table()]
#' @export
zibbmr_results_table <- function(species,
                                 mod1_full, mod1_no_preg,
                                 mod2_full, mod2_no_preg, mod2_no_inter,
                                 df = 2, alpha = 0.05) {
  .saem_results_table(
    species, mod1_full, mod1_no_preg, mod2_full, mod2_no_preg, mod2_no_inter,
    df = df, alpha = alpha
  )
}

#' Preparar datos tipo Romero para ZIBBMR
#'
#' Analogo de [prepare_romero_zibr()], pero conserva los conteos crudos por
#' taxon (en vez de convertirlos a abundancia relativa), listos para usarse
#' con [fit_zibbmr()]/[fit_zibbmr_taxon()] junto con la columna
#' `Total.Read.Counts` como profundidad de secuenciacion.
#'
#' @param romero Una lista con elementos `SampleData` y `OTU`, ver
#'   [prepare_romero_zibr()].
#' @param taxa_out Indices de columnas de taxones a excluir explicitamente
#'   tras el filtro por proporcion de ceros.
#' @param zero_range Vector de largo 2 con el rango `[min, max]` de
#'   proporcion de ceros permitido para retener un taxon (calculado sobre los
#'   conteos crudos).
#'
#' @return Una lista con `data` (covariables + conteos de los taxones
#'   retenidos), `taxa`, `covariates`, `counts`, `taxa_removed` y
#'   `zero_range`.
#' @seealso [prepare_romero_zibr()] para la version en proporciones (ZIBR).
#' @export
prepare_romero_zibbmr <- function(romero, taxa_out = c(31, 49, 50, 60),
                                  zero_range = c(0.1, 0.9)) {
  if (!all(c("SampleData", "OTU") %in% names(romero))) {
    stop("romero debe contener los elementos SampleData y OTU.", call. = FALSE)
  }

  sample_data <- romero$SampleData
  otu <- romero$OTU

  keep_age <- stats::complete.cases(sample_data$Age)
  sample_data <- sample_data[keep_age, , drop = FALSE]
  otu <- otu[keep_age, , drop = FALSE]

  subject_id <- as.numeric(factor(sample_data$Subect_ID, levels = unique(sample_data$Subect_ID)))
  age_sc <- as.numeric(scale(
    sample_data$Age,
    center = min(sample_data$Age),
    scale = max(sample_data$Age) - min(sample_data$Age)
  ))

  month <- ifelse(
    sample_data$pregnant == 1,
    7 * sample_data$GA_Days / 30,
    sample_data$GA_Days / 30
  )

  time <- as.numeric(scale(
    month,
    center = min(month),
    scale = max(month) - min(month)
  ))

  covariates <- data.frame(
    Subect_ID = sample_data$Subect_ID,
    ID = subject_id,
    Time = time,
    pregnant = sample_data$pregnant,
    AGE_SC = age_sc,
    Time_Preg = time * sample_data$pregnant,
    Total.Read.Counts = sample_data$Total.Read.Counts,
    Month = month
  )

  zero_prop <- vapply(otu, function(x) mean(x == 0), numeric(1))
  otu_filtered <- otu[, zero_prop >= zero_range[1] & zero_prop <= zero_range[2], drop = FALSE]

  taxa_def <- setdiff(seq_len(ncol(otu_filtered)), taxa_out)
  taxa_names <- colnames(otu_filtered)[taxa_def]

  list(
    data = cbind(covariates, otu_filtered[, taxa_def, drop = FALSE]),
    taxa = taxa_names,
    covariates = covariates,
    counts = otu_filtered[, taxa_def, drop = FALSE],
    taxa_removed = taxa_out,
    zero_range = zero_range
  )
}


#### Alias de compatibilidad con nombres del codigo original ####

#' Alias historico de fit_zibbmr con la firma del codigo original de Barrera
#'
#' Envoltorio de compatibilidad hacia atras que expone [fit_zibbmr()] con los
#' mismos nombres de argumento que el script original `saem_zibbmr()`
#' (`jbarrera232/saem-zibbmr`). Se mantiene para no romper analisis
#' existentes; el codigo nuevo deberia usar [fit_zibbmr()] directamente.
#'
#' @param Y,X,Z,S,index,zi,v0,a0,b0,seed,iter,ncad,a.fix,b.fix Ver los
#'   argumentos equivalentes de [fit_zibbmr()]: `Y = y`, `index = id`,
#'   `v0 = phi_start`, `a0 = alpha_start`, `b0 = beta_start`, `iter = n_iter`,
#'   `ncad = n_chains`, `a.fix`/`b.fix` equivalen a `alpha_random`/
#'   `beta_random` (`a.fix == 0` marca las posiciones aleatorias).
#' @param compute_fim Ver [fit_zibbmr()].
#'
#' @return Un objeto `zibbmr_saem`, igual que [fit_zibbmr()].
#' @seealso [fit_zibbmr()]
#' @export
saem_zibbmr_clean <- function(Y, X = NULL, Z = NULL, S, index, zi = TRUE,
                              v0, a0 = NULL, b0, seed, iter, ncad = 5,
                              a.fix = NULL, b.fix = NULL, compute_fim = TRUE) {
  fit_zibbmr(
    y = Y,
    S = S,
    id = index,
    X = X,
    Z = Z,
    zi = zi,
    phi_start = v0,
    alpha_start = a0,
    beta_start = b0,
    n_iter = iter,
    n_chains = ncad,
    seed = seed,
    alpha_random = if (is.null(a.fix)) NULL else a.fix == 0,
    beta_random = if (is.null(b.fix)) NULL else b.fix == 0,
    compute_fim = compute_fim
  )
}

#' Alias historico de simulate_zibbmr_data con la firma del codigo original
#'
#' Envoltorio de compatibilidad hacia atras que expone
#' [simulate_zibbmr_data()] con los nombres de argumento del script original
#' de Barrera.
#'
#' @param n.ind,n.obs.ind,s.tot,zi,X,Z,alpha,beta,s1,s2,v,seed Ver los
#'   argumentos equivalentes de [simulate_zibbmr_data()]: `n.ind =
#'   n_subjects`, `n.obs.ind = n_time`, `s.tot = S`, `s1 = sigma_alpha`,
#'   `s2 = sigma_beta`, `v = phi`.
#'
#' @return Un data frame, igual que [simulate_zibbmr_data()].
#' @seealso [simulate_zibbmr_data()]
#' @export
sim_zibbmr_data_clean <- function(n.ind, n.obs.ind, s.tot, zi = TRUE,
                                  X = NULL, Z = NULL, alpha = NULL, beta,
                                  s1 = NULL, s2, v, seed) {
  simulate_zibbmr_data(
    n_subjects = n.ind,
    n_time = n.obs.ind,
    S = s.tot,
    zi = zi,
    X = X,
    Z = Z,
    alpha = alpha,
    beta = beta,
    sigma_alpha = s1,
    sigma_beta = s2,
    phi = v,
    seed = seed
  )
}
