#### Clean SAEM-ZIBR implementation ####
#### Basis for later R-package development ####

#### Internal utilities ####
#### (.saem_check_packages, .saem_diag*, .saem_covariate_matrix, ####
####  .saem_replicate_design, .saem_linear_prob live in R/utils.R, ####
####  shared with ZIBBMR) ####

.zibr_validate_y <- function(y, zi, eps = 1e-6) {
  y <- as.numeric(y)

  if (any(!is.finite(y))) {
    stop("Y contains non-finite values.", call. = FALSE)
  }
  if (any(y < 0 | y > 1)) {
    stop("ZIBR requires proportions in the interval [0, 1].", call. = FALSE)
  }
  if (any(y == 1)) {
    warning("Y contains values exactly equal to 1; they are replaced by 1 - eps to avoid log(0).", call. = FALSE)
    y[y == 1] <- 1 - eps
  }
  if (!zi) {
    y[y == 0] <- eps
  }

  y
}


#### Conditional log-likelihoods used in the M-step ####
#### (.saem_neg_loglik_zero, the zero-inflation part, lives in ####
####  R/utils.R, shared with ZIBBMR) ####

.zibr_neg_loglik_beta <- function(par, psi_chain, beta_random,
                                  z_design_chain, id_chain,
                                  is_positive_chain, y_chain,
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
    lgamma(phi) -
      lgamma(phi * u[is_positive_chain]) -
      lgamma((1 - u[is_positive_chain]) * phi) +
      phi * u[is_positive_chain] * log(y_chain[is_positive_chain]) +
      phi * (1 - u[is_positive_chain]) * log(1 - y_chain[is_positive_chain])
  )

  -loglik
}


#### Importance sampling for the marginal log-likelihood ####

.zibr_loglik_importance <- function(mu, G, phi, zi, y, id,
                                    x_design, z_design, n_alpha, n_beta,
                                    is_positive, is_zero,
                                    psi_mean, psi_var, random_index,
                                    n_random, n_samples = 500, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  G_inv <- .saem_diag_inverse(G)
  # See the note in zibbmr.R: prod(diag(G)) is only the determinant if G is
  # diagonal, and the M-step estimates the full matrix.
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
    log(y[is_positive]) * (phi * u_draws[is_positive, ] - 1) +
    log(1 - y[is_positive]) * (phi * (1 - u_draws[is_positive, ]) - 1) -
    lbeta(phi * u_draws[is_positive, ], phi * (1 - u_draws[is_positive, ]))

  P1 <- rowsum(log_y_given_psi, id)

  mu_random <- matrix(rep(mu[random_index], n_samples), nrow = n_random, ncol = n_samples)
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


#### Gradient and Hessian of the complete log-likelihood ####

.zibr_complete_grad <- function(mu, G, phi, zi, psi_chain,
                                random_index, alpha_random, beta_random,
                                n_random, x_design_chain = NULL, id_chain,
                                is_positive_chain, is_zero_chain,
                                n_alpha, n_beta, z_design_chain,
                                y_chain, n_alpha_random, n_beta_random) {
  G_diag <- .saem_diag(G)
  n_rows <- nrow(psi_chain)

  if (zi) {
    alpha <- mu[seq_len(n_alpha)]
  }
  beta <- mu[n_alpha + seq_len(n_beta)]

  psi_random <- as.matrix(psi_chain[, random_index, drop = FALSE])
  psi_sum <- colSums(psi_random)
  psi_sum2 <- colSums(psi_random^2)

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
    .zibr_neg_loglik_beta,
    beta_phi_par,
    psi_chain = psi_chain,
    beta_random = beta_random,
    z_design_chain = z_design_chain,
    id_chain = id_chain,
    is_positive_chain = is_positive_chain,
    y_chain = y_chain,
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

.zibr_complete_hess <- function(mu, G, phi, zi, psi_chain,
                                random_index, alpha_random, beta_random,
                                n_random, x_design_chain = NULL, id_chain,
                                is_positive_chain, is_zero_chain,
                                n_alpha, n_beta, z_design_chain,
                                y_chain, n_alpha_random, n_beta_random) {
  G_diag <- .saem_diag(G)
  n_rows <- nrow(psi_chain)

  if (zi) {
    alpha <- mu[seq_len(n_alpha)]
  }
  beta <- mu[n_alpha + seq_len(n_beta)]

  psi_random <- as.matrix(psi_chain[, random_index, drop = FALSE])
  psi_sum <- colSums(psi_random)
  psi_sum2 <- colSums(psi_random^2)

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
    .zibr_neg_loglik_beta,
    beta_phi_par,
    psi_chain = psi_chain,
    beta_random = beta_random,
    z_design_chain = z_design_chain,
    id_chain = id_chain,
    is_positive_chain = is_positive_chain,
    y_chain = y_chain,
    n_alpha = n_alpha,
    n_beta = n_beta,
    n_beta_random = n_beta_random
  )

  beta_phi_index <- setdiff(seq_len(n_beta + 1), which(beta_random)) + n_alpha
  H[beta_phi_index, beta_phi_index] <- beta_phi_hess

  H
}


#### Main SAEM-ZIBR fit ####

#' Fit a ZIBR (zero-inflated beta regression) model via SAEM
#'
#' Estimates, by Stochastic Approximation EM (SAEM), a zero-inflated beta
#' mixed-regression model for a response bounded in `[0, 1)` (for example, the
#' relative abundance of a microbiome taxon), following the method described in
#' Barrera (ZIBR: "A stochastic method to estimate a zero-inflated two-part
#' mixed model for human microbiome data"). The model has two parts: a logistic
#' part for the presence probability (`X`/`alpha`) and a beta part for the
#' magnitude conditional on being present (`Z`/`beta`, `phi`), both with a
#' per-subject random intercept.
#'
#' @param y Numeric vector of proportions in `[0, 1)` (the response).
#' @param id Vector (or factor) identifying the subject of each observation in
#'   `y`. Must have the same length as `y`.
#' @param X Matrix or data frame of covariates for the zero-inflation part
#'   (logistic part). `NULL` if `zi = FALSE`.
#' @param Z Matrix or data frame of covariates for the beta part (conditional
#'   magnitude). `NULL` is equivalent to intercept only.
#' @param zi Logical. If `TRUE` (default) fits the zero-inflation part; if
#'   `FALSE`, assumes `y` has no structural zeros.
#' @param phi_start Starting value for the dispersion parameter `phi` of the
#'   beta part.
#' @param alpha_start Vector of starting values for the coefficients of the
#'   logistic part (intercept + columns of `X`). Required if `zi = TRUE`.
#' @param beta_start Vector of starting values for the coefficients of the beta
#'   part (intercept + columns of `Z`).
#' @param n_iter Number of iterations of the SAEM algorithm.
#' @param n_chains Number of parallel MCMC chains used in the stochastic
#'   simulation step (S-step).
#' @param seed Optional random seed.
#' @param alpha_random Logical vector indicating which coefficients of the
#'   logistic part are random effects (by default, only the intercept).
#' @param beta_random Logical vector indicating which coefficients of the beta
#'   part are random effects (by default, only the intercept).
#' @param n_is Number of importance-sampling draws used to estimate the
#'   marginal log-likelihood at the end of the fit.
#' @param cov_random Structure of the random-effect covariance: `"diag"`
#'   (default, independent between the two model parts) or `"unstructured"`
#'   (also estimates the covariance). See [fit_zibbmr()] for the caveats about
#'   `"unstructured"`.
#' @param compute_fim Logical. If `TRUE`, computes the stochastic Fisher
#'   information matrix (needed for `vcov()`/`se()`).
#' @param eps Small value used to avoid `log(0)` when `y` contains values
#'   exactly equal to 0 (with `zi = FALSE`) or exactly equal to 1.
#'
#' @return An object of class `zibr_saem` (and `SAEM_ZIBR_result` for
#'   compatibility), a list containing, among others, the elements `mu` (alpha
#'   and beta concatenated), `G` (random-effect variance), `phi`, `loglik`,
#'   `trace` and `fisher_stoch`. It has [print()], [plot()], [stats::logLik()],
#'   [stats::coef()], [stats::vcov()] and [se()] methods.
#'
#' @seealso [fit_zibr_taxon()] to fit directly on a data-frame column,
#'   [simulate_zibr_data()] to generate test data, [lrt_zibr()] to compare
#'   nested models.
#'
#' @examples
#' \donttest{
#' n_subjects <- 20
#' n_time <- 4
#' n_obs <- n_subjects * n_time
#' dat <- simulate_zibr_data(
#'   n_subjects = n_subjects, n_time = n_time, alpha = c(-0.3, 0.5),
#'   beta = c(0.2, -0.4), sigma_alpha = 0.4, sigma_beta = 0.3, phi = 15,
#'   X = matrix(rbinom(n_obs, 1, 0.5)), Z = matrix(rbinom(n_obs, 1, 0.5)),
#'   seed = 1
#' )
#' fit <- fit_zibr(
#'   y = dat$Y, id = dat$Subject, X = dat$X.1, Z = dat$Z.1,
#'   phi_start = 10, alpha_start = c(-0.2, 0.1), beta_start = c(0.1, 0.1),
#'   n_iter = 50, seed = 1, compute_fim = FALSE
#' )
#' print(fit)
#' }
#' @export
fit_zibr <- function(y, id, X = NULL, Z = NULL, zi = TRUE,
                     phi_start, alpha_start = NULL, beta_start,
                     n_iter = 500, n_chains = 5, seed = NULL,
                     alpha_random = NULL, beta_random = NULL,
                     n_is = 500, compute_fim = TRUE, eps = 1e-6,
                     cov_random = c("diag", "unstructured")) {
  .saem_check_packages(inference = compute_fim)
  cov_random <- .saem_validate_structure(cov_random)

  if (!is.null(seed)) {
    set.seed(seed)
  }

  y <- .zibr_validate_y(y, zi = zi, eps = eps)

  if (length(y) != length(id)) {
    stop("y and id must have the same length.", call. = FALSE)
  }

  n_total <- length(y)
  subject_id <- as.numeric(factor(id, levels = unique(id)))
  n_subjects <- length(unique(subject_id))

  if (zi) {
    x_design <- .saem_covariate_matrix(X, n_total, "X")
    n_alpha <- ncol(x_design)

    if (length(alpha_start) != n_alpha) {
      stop("alpha_start must have length equal to 1 + number of columns of X.", call. = FALSE)
    }

    alpha_random <- if (is.null(alpha_random)) c(TRUE, rep(FALSE, n_alpha - 1)) else alpha_random
    if (length(alpha_random) != n_alpha) {
      stop("alpha_random must have length equal to length(alpha_start).", call. = FALSE)
    }

    x_design_chain <- .saem_replicate_design(x_design, n_chains)
    alpha_labels <- colnames(x_design)
    n_alpha_random <- sum(alpha_random)
  } else {
    if (!is.null(X) || !is.null(alpha_start)) {
      stop("Do not supply X or alpha_start when zi = FALSE.", call. = FALSE)
    }

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
    stop("beta_start must have length equal to 1 + number of columns of Z.", call. = FALSE)
  }

  beta_random <- if (is.null(beta_random)) c(TRUE, rep(FALSE, n_beta - 1)) else beta_random
  if (length(beta_random) != n_beta) {
    stop("beta_random must have length equal to length(beta_start).", call. = FALSE)
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

  # Subject group to collapse the (subject, chain) rows of psi_chain over the
  # chains; each subject appears exactly n_chains times, so the per-subject
  # mean is rowsum(.)/n_chains (equivalent to tapply(., mean) but without
  # rebuilding the factor at each iteration).
  subject_group_chain <- rep(seq_len(n_subjects), n_chains)

  is_positive_chain <- rep(is_positive, n_chains)
  is_zero_chain <- rep(is_zero, n_chains)
  y_chain <- rep(y, n_chains)

  mu <- c(alpha_start, beta_start)
  G_full <- 0.5 * .saem_diag(abs(mu))
  G <- .saem_impose_structure(G_full[random_index, random_index, drop = FALSE], cov_random)
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
        lgamma(phi * u_candidate[is_positive_chain]) +
        lgamma(phi * (1 - u_candidate[is_positive_chain])) +
        log(p_chain[is_positive_chain]) -
        lgamma(phi * u_chain[is_positive_chain]) -
        lgamma(phi * (1 - u_chain[is_positive_chain])) -
        log(p_candidate[is_positive_chain]) +
        phi * (u_chain[is_positive_chain] - u_candidate[is_positive_chain]) *
        (log(y_chain[is_positive_chain]) - log(1 - y_chain[is_positive_chain]))

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
        lgamma(phi * u_candidate[is_positive_chain]) +
        lgamma(phi * (1 - u_candidate[is_positive_chain])) +
        log(p_chain[is_positive_chain]) -
        lgamma(phi * u_chain[is_positive_chain]) -
        lgamma(phi * (1 - u_chain[is_positive_chain])) -
        log(p_candidate[is_positive_chain]) +
        phi * (u_chain[is_positive_chain] - u_candidate[is_positive_chain]) *
        (log(y_chain[is_positive_chain]) - log(1 - y_chain[is_positive_chain]))

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

    proposal_sd_uni <- (1 + 0.4 * (accepted_uni / proposed_uni - 0.4)) * proposal_sd_uni

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
        lgamma(phi * u_candidate[is_positive_chain]) +
        lgamma(phi * (1 - u_candidate[is_positive_chain])) +
        log(p_chain[is_positive_chain]) -
        lgamma(phi * u_chain[is_positive_chain]) -
        lgamma(phi * (1 - u_chain[is_positive_chain])) -
        log(p_candidate[is_positive_chain]) +
        phi * (u_chain[is_positive_chain] - u_candidate[is_positive_chain]) *
        (log(y_chain[is_positive_chain]) - log(1 - y_chain[is_positive_chain]))

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
      (1 + 0.4 * (accepted_multi / (4 * n_chains * n_subjects) - 0.4)) *
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
      G <- .saem_impose_structure(G_full[random_index, random_index, drop = FALSE], cov_random)
      # The proposal kernel uses G_full: it must respect the same restriction,
      # otherwise it proposes in directions the model does not allow.
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
        objective = .zibr_neg_loglik_beta,
        psi_chain = psi_chain,
        beta_random = beta_random,
        z_design_chain = z_design_chain,
        id_chain = id_chain,
        is_positive_chain = is_positive_chain,
        y_chain = y_chain,
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
        grad_current <- .zibr_complete_grad(
          mu, G, phi, zi, psi_chain, random_index, alpha_random,
          beta_random, n_random, x_design_chain, id_chain,
          is_positive_chain, is_zero_chain, n_alpha, n_beta,
          z_design_chain, y_chain, n_alpha_random, n_beta_random
        )

        hess_current <- .zibr_complete_hess(
          mu, G, phi, zi, psi_chain, random_index, alpha_random,
          beta_random, n_random, x_design_chain, id_chain,
          is_positive_chain, is_zero_chain, n_alpha, n_beta,
          z_design_chain, y_chain, n_alpha_random, n_beta_random
        )

        score2_current <- matrix(0, nrow = n_psi + n_random + 1, ncol = n_psi + n_random + 1)

        for (chain in seq_len(n_chains)) {
          row_index <- n_subjects * (chain - 1) + seq_len(n_subjects)

          grad_chain <- .zibr_complete_grad(
            mu, G, phi, zi, psi_chain[row_index, , drop = FALSE],
            random_index, alpha_random, beta_random, n_random,
            x_design, subject_id, is_positive, is_zero,
            n_alpha, n_beta, z_design, y,
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

    trace <- rbind(trace, c(mu, .saem_diag(G), phi))

    if (zi) {
      colnames(trace) <- c(
        paste("A", seq_len(n_alpha), sep = "."),
        paste("B", seq_len(n_beta), sep = "."),
        paste("SIGMA", seq_len(n_random), sep = "."),
        "V"
      )
    } else {
      colnames(trace) <- c(
        paste("B", seq_len(n_beta), sep = "."),
        paste("SIGMA", seq_len(n_random), sep = "."),
        "V"
      )
    }
  }

  psi_mean <- rowsum(psi_chain, subject_group_chain) / n_chains
  psi_var <- psi_subject_second - psi_subject_mean^2
  psi_var[, -random_index] <- 0

  loglik <- .zibr_loglik_importance(
    mu = mu,
    G = G,
    phi = phi,
    zi = zi,
    y = y,
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
    # Original data, to be able to compute predictions and residuals in plot().
    # They do not affect estimation; they are only stored for the plots.
    data = list(y = y, x_design = x_design, z_design = z_design,
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

  class(out) <- c("zibr_saem", "SAEM_ZIBR_result")
  out
}


#### ZIBR data simulation ####

#' Simulate longitudinal data for a ZIBR model
#'
#' Generates a data frame of longitudinal data compatible with [fit_zibr()]:
#' a continuous response bounded in `[0, 1]`, with optional zero inflation,
#' fixed effects and a per-subject random intercept (or other coefficients)
#' with a multivariate normal distribution.
#'
#' @param n_subjects Number of subjects.
#' @param n_time Number of observations (time points) per subject.
#' @param zi Logical. If `TRUE` (default), simulates presence/absence with a
#'   logistic model (`X`/`alpha`) before simulating the beta magnitude.
#' @param X Matrix or data frame of covariates for the zero-inflation part.
#'   `NULL` if `zi = FALSE`.
#' @param Z Matrix or data frame of covariates for the beta part.
#' @param alpha Vector of true coefficients of the logistic part (intercept +
#'   columns of `X`). Required if `zi = TRUE`.
#' @param beta Vector of true coefficients of the beta part (intercept +
#'   columns of `Z`).
#' @param sigma_alpha Standard deviation of the random intercept of the
#'   logistic part. Required if `zi = TRUE`.
#' @param sigma_beta Standard deviation of the random intercept of the beta
#'   part.
#' @param phi Dispersion parameter of the beta distribution.
#' @param seed Optional random seed.
#'
#' @return A data frame with columns `Subject`, `Time`, `Y` (the simulated
#'   response) and the `X`/`Z` covariates used.
#'
#' @seealso [fit_zibr()]
#' @examples
#' n_subjects <- 10
#' n_time <- 3
#' n_obs <- n_subjects * n_time
#' dat <- simulate_zibr_data(
#'   n_subjects = n_subjects, n_time = n_time, alpha = c(-0.3, 0.5),
#'   beta = c(0.2, -0.4), sigma_alpha = 0.4, sigma_beta = 0.3, phi = 15,
#'   X = matrix(rbinom(n_obs, 1, 0.5)), Z = matrix(rbinom(n_obs, 1, 0.5)),
#'   seed = 1
#' )
#' head(dat)
#' @export
simulate_zibr_data <- function(n_subjects, n_time, zi = TRUE,
                               X = NULL, Z = NULL, alpha = NULL, beta,
                               sigma_alpha = NULL, sigma_beta,
                               phi, seed = NULL) {
  .saem_check_packages(inference = FALSE)

  if (!is.null(seed)) {
    set.seed(seed)
  }

  n_total <- n_subjects * n_time

  if (zi) {
    x_design <- .saem_covariate_matrix(X, n_total, "X")
    n_alpha <- ncol(x_design)

    if (length(alpha) != n_alpha) {
      stop("alpha must have length equal to 1 + number of columns of X.", call. = FALSE)
    }

    random_cols <- c(1, n_alpha + 1)
  } else {
    if (!is.null(X) || !is.null(alpha)) {
      stop("Do not supply X or alpha when zi = FALSE.", call. = FALSE)
    }

    x_design <- NULL
    n_alpha <- 0
    random_cols <- 1
  }

  z_design <- .saem_covariate_matrix(Z, n_total, "Z")
  n_beta <- ncol(z_design)

  if (length(beta) != n_beta) {
    stop("beta must have length equal to 1 + number of columns of Z.", call. = FALSE)
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

  present <- stats::rbinom(n_total, size = 1, prob = p)
  y_positive <- stats::rbeta(n_total, shape1 = u * phi, shape2 = (1 - u) * phi)
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
    covariates,
    check.names = FALSE
  )
}


#### Basic methods ####

#' Print a ZIBR fit
#'
#' @param x A `zibr_saem` object, the result of [fit_zibr()].
#' @param ... Unused, for compatibility with the [print()] generic.
#' @return `x`, invisibly.
#' @export
print.zibr_saem <- function(x, ...) {
  .saem_print(x, model_label = "SAEM-ZIBR", beta_label = "Beta part")
}

#' Plots of a ZIBR fit
#'
#' Produces different diagnostic/result plots depending on `which`:
#' \describe{
#'   \item{`"convergence"`}{(default) iteration-by-iteration trace of the
#'     parameters, to check the convergence of the SAEM algorithm.}
#'   \item{`"coefficients"`}{estimated coefficients with their 95% confidence
#'     interval (forest-plot style). Requires fitting with
#'     `compute_fim = TRUE` to show the intervals.}
#'   \item{`"random"`}{between-subject distribution of the estimated random
#'     effects (one per subject); the red line marks the population mean.}
#'   \item{`"fit"`}{observed vs. predicted for the continuous part, on the
#'     positive observations (where the taxon is present), with the reference
#'     line `y = x`. It focuses on the continuous part because in a
#'     zero-inflated model the marginal prediction mixes the mass at zero.}
#'   \item{`"residuals"`}{residuals of the continuous part (observed minus
#'     individual predicted, on positive observations): their spread against
#'     the predicted value and their distribution.}
#' }
#' The `"fit"` and `"residuals"` plots use the original data that the fit
#' stores; they are only available for fits made with a version of the package
#' that stores them.
#'
#' @param x A `zibr_saem` object, the result of [fit_zibr()].
#' @param which Type of plot: `"convergence"`, `"coefficients"`, `"random"`,
#'   `"fit"` or `"residuals"`.
#' @param ... Additional arguments (unused for now).
#' @return `x`, invisibly. Called for its side effect of plotting.
#' @export
plot.zibr_saem <- function(x, which = c("convergence", "coefficients",
                                        "random", "fit", "residuals"), ...) {
  .saem_plot(x, which = match.arg(which), beta_label = "beta", ...)
}

#' Marginal log-likelihood of a ZIBR fit
#'
#' @param object A `zibr_saem` object, the result of [fit_zibr()].
#' @param ... Unused, for compatibility with the [stats::logLik()] generic.
#' @return A `logLik` object with the marginal log-likelihood estimated by
#'   importance sampling, with attributes `df` (degrees of freedom) and `nobs`.
#' @export
logLik.zibr_saem <- function(object, ...) {
  .saem_logLik(object)
}

#' Estimated coefficients of a ZIBR fit
#'
#' @param object A `zibr_saem` object, the result of [fit_zibr()].
#' @param ... Unused, for compatibility with the [stats::coef()] generic.
#' @return Numeric vector `mu` with the coefficients of the logistic part
#'   followed by those of the beta part.
#' @export
coef.zibr_saem <- function(object, ...) {
  .saem_coef(object)
}

#' Variance-covariance matrix of a ZIBR fit
#'
#' Computes the inverse of the stochastic Fisher information matrix
#' (`fisher_stoch`), estimated during the SAEM fit.
#'
#' @param object A `zibr_saem` object fitted with `compute_fim = TRUE`.
#' @param ... Unused, for compatibility with the [stats::vcov()] generic.
#' @return A variance-covariance matrix.
#' @export
vcov.zibr_saem <- function(object, ...) {
  .saem_vcov(object)
}

#' @rdname se
#' @export
se.zibr_saem <- function(object, ...) {
  .saem_se(object)
}


#### Clean functions for per-taxon analysis ####

#' Fit ZIBR for one taxon of a data frame
#'
#' Wrapper around [fit_zibr()] meant to work directly on a long-format
#' microbiome data frame (one row per subject-time, one column per taxon).
#' Extracts the taxon column and the covariates by name, generates reasonable
#' random starting values if none are supplied, and calls [fit_zibr()].
#'
#' @param data Data frame with one row per observation, including the taxon
#'   column, the id column and the covariates.
#' @param taxon Name of the column in `data` with the taxon proportion to be
#'   modelled.
#' @param covariates Vector of column names to use as covariates in both the
#'   logistic part and the beta part. Ignored if `x_covariates`/`z_covariates`
#'   are supplied separately.
#' @param x_covariates Column names for the zero-inflation part (by default,
#'   equal to `covariates`).
#' @param z_covariates Column names for the beta part (by default, equal to
#'   `covariates`).
#' @param id Name of the column in `data` that identifies the subject.
#' @param zi Logical, see [fit_zibr()].
#' @param phi_start Starting value of `phi`. If `NULL`, it is drawn with
#'   `runif(1, 10, 20)`.
#' @param alpha_start Starting values of the logistic part. If `NULL`, they are
#'   drawn with `runif(., -0.1, 0.1)`.
#' @param beta_start Starting values of the beta part. If `NULL`, they are
#'   drawn with `runif(., -0.1, 0.1)`.
#' @param seed Random seed, used both for the starting values and for the SAEM
#'   fit.
#' @param n_iter Number of SAEM iterations.
#' @param n_chains Number of MCMC chains.
#' @param compute_fim Logical, see [fit_zibr()].
#' @param ... Additional arguments passed to [fit_zibr()].
#'
#' @return A `zibr_saem` object, same as [fit_zibr()].
#' @seealso [fit_zibr()], [fit_zibr_taxa()]
#' @export
fit_zibr_taxon <- function(data, taxon, covariates = NULL,
                           x_covariates = covariates,
                           z_covariates = covariates,
                           id, zi = TRUE, phi_start = NULL,
                           alpha_start = NULL, beta_start = NULL,
                           seed = 232, n_iter = 500, n_chains = 5,
                           compute_fim = FALSE, ...) {
  if (!taxon %in% names(data)) {
    stop("The specified taxon does not exist in data.", call. = FALSE)
  }
  if (!id %in% names(data)) {
    stop("The id column does not exist in data.", call. = FALSE)
  }

  if (is.null(x_covariates)) {
    x_covariates <- character(0)
  }
  if (is.null(z_covariates)) {
    z_covariates <- character(0)
  }

  if (!all(x_covariates %in% names(data))) {
    stop("At least one covariate in x_covariates does not exist in data.", call. = FALSE)
  }
  if (!all(z_covariates %in% names(data))) {
    stop("At least one covariate in z_covariates does not exist in data.", call. = FALSE)
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

  fit_zibr(
    y = data[[taxon]],
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

#' Fit ZIBR for several taxa of a data frame
#'
#' Applies [fit_zibr_taxon()] to each element of `taxa`, with the same
#' covariate and iteration configuration for all of them.
#'
#' @param data Data frame with one row per observation.
#' @param taxa Vector of column names (taxa) to fit.
#' @param covariates,x_covariates,z_covariates See [fit_zibr_taxon()].
#' @param id Name of the column that identifies the subject.
#' @param zi Logical, see [fit_zibr()].
#' @param seed Random seed (reused for each taxon).
#' @param n_iter Number of SAEM iterations.
#' @param n_chains Number of MCMC chains.
#' @param compute_fim Logical, see [fit_zibr()].
#' @param ... Additional arguments passed to [fit_zibr_taxon()].
#'
#' @return A list of `zibr_saem` objects, named according to `taxa`.
#' @seealso [fit_zibr_taxon()]
#' @export
fit_zibr_taxa <- function(data, taxa, covariates = NULL,
                          x_covariates = covariates,
                          z_covariates = covariates,
                          id, zi = TRUE, seed = 232,
                          n_iter = 500, n_chains = 5,
                          compute_fim = FALSE, ...) {
  fits <- lapply(taxa, function(taxon) {
    fit_zibr_taxon(
      data = data,
      taxon = taxon,
      x_covariates = x_covariates,
      z_covariates = z_covariates,
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

#' Likelihood-ratio test between two nested ZIBR fits
#'
#' Compares two nested ZIBR models (for example, with and without a covariate)
#' via a likelihood-ratio test using the marginal log-likelihood (importance
#' sampling) of each fit.
#'
#' @param full Full model: a `zibr_saem` object or any object with a
#'   [stats::logLik()] method.
#' @param reduced Reduced model (nested in `full`), same type as `full`.
#' @param df Degrees of freedom of the test (number of extra parameters in
#'   `full` relative to `reduced`).
#'
#' @return A one-row data frame with `LL_full`, `LL_reduced`, `LRT` (the
#'   statistic), `df` and `p_value`.
#' @seealso [lrt_zibr_table()]
#' @export
lrt_zibr <- function(full, reduced, df = 2) {
  .saem_lrt(full, reduced, df = df)
}

#' Table of likelihood-ratio tests for several ZIBR taxa
#'
#' Applies [lrt_zibr()] pairing each element of `full_models` with the
#' corresponding one in `reduced_models`, and builds a table with one row per
#' taxon.
#'
#' @param full_models List of full models (one per taxon).
#' @param reduced_models List of reduced models (one per taxon, same order as
#'   `full_models`).
#' @param species Vector of names/labels for each taxon (by default,
#'   `names(full_models)`).
#' @param df Degrees of freedom of the test, see [lrt_zibr()].
#' @param alpha Significance level used to flag `Detected`.
#'
#' @return A data frame with one row per taxon: `Species`, the columns of
#'   [lrt_zibr()], and `Detected` (logical, `p_value < alpha`).
#' @seealso [lrt_zibr()]
#' @export
lrt_zibr_table <- function(full_models, reduced_models, species = names(full_models),
                           df = 2, alpha = 0.05) {
  .saem_lrt_table(
    full_models, reduced_models,
    species = species, df = df, alpha = alpha, lrt_fn = lrt_zibr
  )
}

#' Summary table of three typical LRT comparisons for ZIBR
#'
#' Builds a results table with two likelihood-ratio tests on a main effect
#' (`mod1` vs. `mod2`, e.g. pregnancy) and an interaction test (`mod2` with and
#' without an interaction term), replicating the reporting format used in the
#' Romero-type data analysis.
#'
#' @param species Vector of names/labels for each taxon.
#' @param mod1_full,mod1_no_preg Lists of ZIBR models (with and without the
#'   main effect) for the first comparison, one per taxon.
#' @param mod2_full,mod2_no_preg,mod2_no_inter Lists of ZIBR models for the
#'   second comparison (main effect) and the interaction test, one per taxon.
#' @param df Degrees of freedom used in the three tests.
#' @param alpha Significance level used for the `Detec_*` columns.
#'
#' @return A data frame with one row per taxon, with the log-likelihoods,
#'   p-values and detections of the three comparisons.
#' @export
zibr_results_table <- function(species,
                               mod1_full, mod1_no_preg,
                               mod2_full, mod2_no_preg, mod2_no_inter,
                               df = 2, alpha = 0.05) {
  .saem_results_table(
    species, mod1_full, mod1_no_preg, mod2_full, mod2_no_preg, mod2_no_inter,
    df = df, alpha = alpha
  )
}


#### Romero data preparation for ZIBR ####

#' Prepare Romero-type data for ZIBR
#'
#' Takes a list with microbiome data in the format of the Romero study (with
#' elements `SampleData` and `OTU`) and builds standard covariates (scaled
#' gestational time, scaled age, time-pregnancy interaction) together with
#' relative abundances filtered by zero proportion, ready to use with
#' [fit_zibr()]/[fit_zibr_taxon()].
#'
#' @param romero A list with elements `SampleData` (per-sample covariates,
#'   including `Age`, `Subect_ID`, `pregnant`, `GA_Days`, `Total.Read.Counts`)
#'   and `OTU` (matrix or data frame of counts per taxon, same number of rows
#'   as `SampleData`).
#' @param taxa_out Indices of taxon columns to exclude explicitly after the
#'   zero-proportion filter.
#' @param zero_range Length-2 vector with the `[min, max]` range of zero
#'   proportion allowed to retain a taxon.
#'
#' @return A list with `data` (covariates + relative abundances of the
#'   retained taxa), `taxa` (names of those taxa), `covariates`, `abundances`,
#'   `taxa_removed` and `zero_range`.
#' @seealso [prepare_romero_zibbmr()] for the count version (ZIBBMR).
#' @export
prepare_romero_zibr <- function(romero, taxa_out = c(31, 49, 50, 60),
                                zero_range = c(0.1, 0.9)) {
  if (!all(c("SampleData", "OTU") %in% names(romero))) {
    stop("romero must contain the elements SampleData and OTU.", call. = FALSE)
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

  rel_abund <- as.data.frame(sweep(otu, 1, sample_data$Total.Read.Counts, FUN = "/"))
  zero_prop <- vapply(rel_abund, function(x) mean(x == 0), numeric(1))
  taxa_filtered <- rel_abund[, zero_prop >= zero_range[1] & zero_prop <= zero_range[2], drop = FALSE]

  taxa_def <- setdiff(seq_len(ncol(taxa_filtered)), taxa_out)
  taxa_names <- colnames(taxa_filtered)[taxa_def]

  list(
    data = cbind(covariates, taxa_filtered[, taxa_def, drop = FALSE]),
    taxa = taxa_names,
    covariates = covariates,
    abundances = taxa_filtered[, taxa_def, drop = FALSE],
    taxa_removed = taxa_out,
    zero_range = zero_range
  )
}


#### Compatibility aliases with the original code's names ####

#' Historical alias of fit_zibr with the signature of Barrera's original code
#'
#' Backward-compatibility wrapper that exposes [fit_zibr()] with the same
#' argument names as the original `saem_zibr()` script
#' (`jbarrera232/saem-zibr`). Kept so as not to break existing analyses; new
#' code should use [fit_zibr()] directly.
#'
#' @param Y,X,Z,index,zi,v0,a0,b0,seed,iter,ncad,a.fix,b.fix See the equivalent
#'   arguments of [fit_zibr()]: `Y = y`, `index = id`, `v0 = phi_start`,
#'   `a0 = alpha_start`, `b0 = beta_start`, `iter = n_iter`, `ncad = n_chains`;
#'   `a.fix`/`b.fix` correspond to `alpha_random`/`beta_random` (`a.fix == 0`
#'   marks the random positions).
#' @param compute_fim See [fit_zibr()].
#'
#' @return A `zibr_saem` object, same as [fit_zibr()].
#' @seealso [fit_zibr()]
#' @export
saem_zibr_clean <- function(Y, X = NULL, Z = NULL, index, zi = TRUE,
                            v0, a0 = NULL, b0, seed, iter = 500, ncad = 5,
                            a.fix = NULL, b.fix = NULL, compute_fim = TRUE) {
  fit_zibr(
    y = Y,
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

#' Historical alias of simulate_zibr_data with the original code's signature
#'
#' Backward-compatibility wrapper that exposes [simulate_zibr_data()] with the
#' argument names of Barrera's original script.
#'
#' @param n.ind,n.obs.ind,zi,X,Z,alpha,beta,s1,s2,v,seed See the equivalent
#'   arguments of [simulate_zibr_data()]: `n.ind = n_subjects`,
#'   `n.obs.ind = n_time`, `s1 = sigma_alpha`, `s2 = sigma_beta`, `v = phi`.
#'
#' @return A data frame, same as [simulate_zibr_data()].
#' @seealso [simulate_zibr_data()]
#' @export
sim_zibr_data_clean <- function(n.ind, n.obs.ind, zi = TRUE,
                                X = NULL, Z = NULL, alpha = NULL, beta,
                                s1 = NULL, s2, v, seed) {
  simulate_zibr_data(
    n_subjects = n.ind,
    n_time = n.obs.ind,
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
