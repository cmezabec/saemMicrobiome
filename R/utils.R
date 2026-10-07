#### Internal utilities shared by ZIBR and ZIBBMR ####
#### These functions do not depend on each model's likelihood: ####
#### design-matrix construction, logistic linear predictor, ####
#### MCMC chain replication and required-package checks. ####

.saem_check_packages <- function(inference = TRUE) {
  required <- "MASS"
  if (inference) {
    required <- c(required, "numDeriv")
  }

  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]

  if (length(missing) > 0) {
    stop(
      "Missing required packages: ",
      paste(missing, collapse = ", "),
      ". Install them before fitting the model.",
      call. = FALSE
    )
  }
}

.saem_diag <- function(x) {
  if (length(x) == 1) {
    matrix(x, nrow = 1, ncol = 1)
  } else {
    diag(x)
  }
}

# Inverse of the random-effect covariance matrix.
#
# The previous version computed diag(1 / diag(G)), i.e. it inverted G
# discarding the off-diagonal terms. With diagonal G the result is the same,
# but the SAEM M-step estimates the FULL matrix (see the computation of
# `G_full` in zibbmr.R and zibr.R), so as soon as the random effects of the
# two model parts are correlated the precision came out wrong.
#
# This affected the Metropolis-Hastings acceptance ratio, which evaluated the
# prior density as if the random effects were independent: the traces came out
# uncorrelated, the accumulated covariance did not grow and G went back to
# being almost diagonal at the next iteration. Result: a self-reinforcing fixed
# point with rho ~ 0 (with a true rho of 0.7 it estimated ~0.09). It also
# affected the marginal likelihood via importance sampling, and with it the LRT.
.saem_diag_inverse <- function(G) {
  G <- as.matrix(G)

  inverse <- try(solve(G), silent = TRUE)

  # If G is singular or ill-conditioned (which can happen in early SAEM
  # iterations, before the variance components stabilize) we fall back to the
  # pseudo-inverse instead of failing.
  if (inherits(inverse, "try-error")) {
    eig <- eigen(G, symmetric = TRUE)
    positive <- eig$values >
      max(1e-10, max(eig$values) * 1e-10)

    if (!any(positive)) {
      stop("The random-effect covariance matrix is singular.",
           call. = FALSE)
    }

    inverse <- eig$vectors[, positive, drop = FALSE] %*%
      diag(1 / eig$values[positive],
           nrow = sum(positive)) %*%
      t(eig$vectors[, positive, drop = FALSE])
  }

  inverse
}

# Imposes the requested structure on the random-effect covariance.
#
# "diag"          random effects independent between the two model parts. It
#                 is the specification of the paper and the default: with few
#                 subjects or observations per subject, the unrestricted
#                 matrix degenerates (rho -> +-1 and collapse of one of the
#                 variances).
# "unstructured"  free covariance. Enables estimating the correlation between
#                 the random effect of the zero-inflation part and that of the
#                 abundance part, which neither glmmTMB nor gamlss can
#                 represent. Requires substantially more information.
.saem_impose_structure <- function(G, structure) {
  if (identical(structure, "diag")) {
    return(.saem_diag(diag(as.matrix(G))))
  }
  as.matrix(G)
}

.saem_validate_structure <- function(structure) {
  match.arg(structure, c("diag", "unstructured"))
}

.saem_covariate_matrix <- function(x, n, prefix) {
  if (is.null(x)) {
    mat <- matrix(nrow = n, ncol = 0)
  } else {
    mat <- as.matrix(x)
    if (nrow(mat) != n) {
      stop("The covariate matrix does not have the same number of rows as Y.", call. = FALSE)
    }
  }

  if (ncol(mat) == 0) {
    out <- matrix(1, nrow = n, ncol = 1)
    colnames(out) <- "Intercept"
    return(out)
  }

  if (is.null(colnames(mat))) {
    colnames(mat) <- paste(prefix, seq_len(ncol(mat)), sep = ".")
  }

  cbind(Intercept = 1, mat)
}

.saem_replicate_design <- function(design, n_chains) {
  do.call("rbind", replicate(n_chains, design, simplify = FALSE))
}

.saem_linear_prob <- function(psi, cols, id, design) {
  # The linear predictor is computed in C++ (avoids materializing psi[id, cols]
  # and the rowSums); plogis() is applied in R (vectorized) so the result is
  # byte-identical to the pure-R version.
  eta <- saem_linear_eta_cpp(psi, as.integer(cols), as.integer(id), design)
  # plogis() returns exactly 0 or 1 once |eta| exceeds about 745 or 37, which
  # happens when a fixed effect is pushed far out during the burn-in (with few
  # positive counts the M-step on a single simulated sample can be separated).
  # A probability of exactly 0 or 1 makes lgamma(phi * u) or
  # lgamma(phi * (1 - u)) infinite, and Inf - Inf = NaN then reaches the
  # Metropolis-Hastings ratio, the chains and the sufficient statistics. Only
  # those two values are moved to the nearest representable probability;
  # everything strictly inside (0, 1) is left untouched.
  p <- plogis(eta)
  p[p == 0] <- .Machine$double.xmin
  p[p == 1] <- 1 - .Machine$double.neg.eps
  p
}


#### Zero-inflation part: identical for ZIBR and ZIBBMR ####
#### (the presence/absence probability does not depend on whether the ####
#### positive part is beta or beta-binomial) ####

.saem_neg_loglik_zero <- function(alpha_fixed, psi_chain, alpha_random,
                                  x_design_chain, id_chain, is_positive_chain,
                                  is_zero_chain, n_alpha) {
  fixed_index <- which(!alpha_random)
  n_rows <- nrow(psi_chain)

  psi_chain[, fixed_index] <- matrix(
    rep(alpha_fixed, each = n_rows),
    ncol = length(alpha_fixed),
    nrow = n_rows
  )

  p <- .saem_linear_prob(psi_chain, seq_len(n_alpha), id_chain, x_design_chain)
  loglik <- sum(log(1 - p[is_zero_chain])) + sum(log(p[is_positive_chain]))

  -loglik
}


#### Log-likelihood extraction and nested-model comparison (LRT) ####
#### Common to both models: they only depend on the object having $loglik ####
#### or a generic logLik() method. ####

.saem_extract_loglik <- function(model) {
  if (inherits(model, c("zibr_saem", "zibbmr_saem"))) {
    return(model$loglik)
  }

  as.numeric(stats::logLik(model))
}

.saem_lrt <- function(full, reduced, df = 2) {
  ll_full <- .saem_extract_loglik(full)
  ll_reduced <- .saem_extract_loglik(reduced)
  statistic <- 2 * (ll_full - ll_reduced)

  data.frame(
    LL_full = ll_full,
    LL_reduced = ll_reduced,
    LRT = statistic,
    df = df,
    p_value = stats::pchisq(statistic, df = df, lower.tail = FALSE)
  )
}

.saem_lrt_table <- function(full_models, reduced_models, species, df, alpha, lrt_fn) {
  if (length(full_models) != length(reduced_models)) {
    stop("full_models and reduced_models must have the same length.", call. = FALSE)
  }

  if (is.null(species) || length(species) == 0) {
    species <- seq_along(full_models)
  }

  out <- do.call(
    rbind,
    Map(function(full, reduced, sp) {
      res <- lrt_fn(full, reduced, df = df)
      data.frame(Species = sp, res, row.names = NULL)
    }, full_models, reduced_models, species)
  )

  out$Detected <- out$p_value < alpha
  rownames(out) <- NULL
  out
}

.saem_results_table <- function(species,
                                mod1_full, mod1_no_preg,
                                mod2_full, mod2_no_preg, mod2_no_inter,
                                df, alpha) {
  LL1.0 <- vapply(mod1_full, .saem_extract_loglik, numeric(1))
  LL1.1 <- vapply(mod1_no_preg, .saem_extract_loglik, numeric(1))
  LL2.0 <- vapply(mod2_full, .saem_extract_loglik, numeric(1))
  LL2.1 <- vapply(mod2_no_preg, .saem_extract_loglik, numeric(1))
  LL2.2 <- vapply(mod2_no_inter, .saem_extract_loglik, numeric(1))

  pval_Preg1 <- stats::pchisq(2 * (LL1.0 - LL1.1), df = df, lower.tail = FALSE)
  pval_Preg2 <- stats::pchisq(2 * (LL2.0 - LL2.1), df = df, lower.tail = FALSE)
  pval_Inter <- stats::pchisq(2 * (LL2.0 - LL2.2), df = df, lower.tail = FALSE)

  data.frame(
    Species = species,
    LL1.0 = LL1.0,
    LL1.1 = LL1.1,
    LL2.0 = LL2.0,
    LL2.1 = LL2.1,
    LL2.2 = LL2.2,
    pval_Preg1 = pval_Preg1,
    pval_Preg2 = pval_Preg2,
    pval_Inter = pval_Inter,
    Detec_Preg1 = pval_Preg1 < alpha,
    Detec_Preg2 = pval_Preg2 < alpha,
    Detec_Inter = pval_Inter < alpha,
    row.names = NULL
  )
}


#### S3 methods: the mechanics of printing/plotting/extracting coefficients ####
#### are identical for zibr_saem and zibbmr_saem; only the text labels change ####

.saem_print <- function(x, model_label, beta_label) {
  cat("===== Results ", model_label, " =====\n", sep = "")

  if (x$zi) {
    alpha <- x$mu[seq_len(x$n_alpha)]
    alpha_tab <- data.frame(
      Estimate = alpha,
      Type = ifelse(x$alpha_random, "Random", "Fixed"),
      Variance = 0,
      sqrt.Var = 0,
      row.names = x$alpha_labels
    )

    n_alpha_random <- sum(x$alpha_random)
    if (n_alpha_random > 0) {
      alpha_tab[x$alpha_random, "Variance"] <- diag(x$G)[seq_len(n_alpha_random)]
      alpha_tab[, "sqrt.Var"] <- sqrt(alpha_tab[, "Variance"])
    }

    cat("== Logistic part: p_it ==\n")
    print(alpha_tab[, c("Estimate", "Type")])
  } else {
    n_alpha_random <- 0
  }

  beta <- x$mu[x$n_alpha + seq_len(x$n_beta)]
  beta_tab <- data.frame(
    Estimate = beta,
    Type = ifelse(x$beta_random, "Random", "Fixed"),
    Variance = 0,
    sqrt.Var = 0,
    row.names = x$beta_labels
  )

  n_beta_random <- sum(x$beta_random)
  if (n_beta_random > 0) {
    beta_tab[x$beta_random, "Variance"] <- diag(x$G)[n_alpha_random + seq_len(n_beta_random)]
    beta_tab[, "sqrt.Var"] <- sqrt(beta_tab[, "Variance"])
  }

  cat("== ", beta_label, ": u_it ==\n", sep = "")
  print(beta_tab[, c("Estimate", "Type")])

  cat("=== Random-effect variances ===\n")
  if (x$zi && n_alpha_random > 0) {
    cat("== Logistic part ==\n")
    print(alpha_tab[x$alpha_random, c("Variance", "sqrt.Var"), drop = FALSE])
  }
  if (n_beta_random > 0) {
    cat("== ", beta_label, " ==\n", sep = "")
    print(beta_tab[x$beta_random, c("Variance", "sqrt.Var"), drop = FALSE])
  }

  cat("=== Phi: ", x$phi, "\n", sep = "")
  cat("=== Marginal log-likelihood (importance sampling): ", x$loglik, "\n", sep = "")

  invisible(x)
}

## Labels of the fixed-effect parameters (logistic part + beta/beta-binomial
## part), in the same order as x$mu.
.saem_param_labels <- function(x, beta_label = "beta") {
  labs <- character(0)
  if (isTRUE(x$zi) && x$n_alpha > 0) {
    labs <- paste0("logistic: ", x$alpha_labels)
  }
  c(labs, paste0(beta_label, ": ", x$beta_labels))
}

## Plot 1: convergence trace (parameters across the iterations).
.saem_plot_trace <- function(x, ...) {
  trace <- x$trace
  n_iter <- nrow(trace)
  burn_in <- floor(0.75 * n_iter)
  n_panels <- ncol(trace)
  n_rows <- ceiling(n_panels / 3)

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)

  graphics::par(mfrow = c(n_rows, 3))
  for (j in seq_len(n_panels)) {
    graphics::plot(
      seq_len(n_iter),
      trace[, j],
      type = "l",
      xlab = "Iteration",
      ylab = "Value",
      main = colnames(trace)[j]
    )
    graphics::abline(v = burn_in, lty = 2)
  }

  invisible(x)
}

## Plot 2: estimated coefficients with 95% confidence interval (forest-plot
## style). Needs the fit with compute_fim = TRUE for the CIs.
.saem_plot_coef <- function(x, beta_label = "beta") {
  est <- x$mu
  labs <- .saem_param_labels(x, beta_label)
  se_all <- tryCatch(suppressWarnings(.saem_se(x)), error = function(e) NULL)
  se_coef <- if (!is.null(se_all) && length(se_all) >= length(est)) {
    se_all[seq_along(est)]
  } else {
    rep(NA_real_, length(est))
  }
  lo <- est - 1.96 * se_coef
  hi <- est + 1.96 * se_coef
  n <- length(est)

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::par(mar = c(4, 12, 3, 1))

  xr <- range(c(est, lo, hi, 0), na.rm = TRUE)
  graphics::plot(est, seq_len(n), xlim = xr, ylim = c(0.5, n + 0.5), yaxt = "n",
                 xlab = "Estimate (95% CI)", ylab = "", pch = 19,
                 main = "Estimated coefficients")
  graphics::axis(2, at = seq_len(n), labels = labs, las = 1, cex.axis = 0.85)
  graphics::abline(v = 0, lty = 2, col = "gray50")
  ok <- !is.na(lo)
  if (any(ok)) graphics::segments(lo[ok], which(ok), hi[ok], which(ok), lwd = 2)
  if (!all(ok)) {
    note <- if (is.null(x$fisher_stoch)) {
      "No CI: refit with compute_fim = TRUE"
    } else {
      "CI not available for some coefficient (ill-conditioned information matrix)"
    }
    graphics::mtext(note, side = 1, line = 2.5, cex = 0.75, col = "gray40")
  }
  invisible(x)
}

## Plot 3: between-subject distribution of the estimated random effects (one
## per subject). The red line marks the population mean.
.saem_plot_random <- function(x, beta_label = "beta") {
  ri <- x$random_index
  if (length(ri) == 0) {
    message("The fit has no random effects to plot.")
    return(invisible(x))
  }
  labs <- .saem_param_labels(x, beta_label)[ri]
  vals <- x$psi_mean[, ri, drop = FALSE]
  k <- ncol(vals)

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::par(mfrow = c(1, k))
  for (j in seq_len(k)) {
    graphics::hist(vals[, j], main = labs[j], xlab = "Value per subject",
                   col = "#92c5de", border = "white")
    graphics::abline(v = x$mu[ri[j]], col = "red", lwd = 2)
  }
  invisible(x)
}

## Predictions of the CONTINUOUS PART of the model (the magnitude given that
## the taxon is present), using the original data that the fit stores in
## `x$data`. It focuses on the continuous part -not the marginal E[Y] = p * u-
## because in a zero-inflated model the marginal prediction mixes the point
## mass at zero with the continuous part and does not read like a classic
## observed-vs-predicted plot. The conditional mean given presence is:
##   ZIBR   : u   (the mean of the beta part, in [0, 1))
##   ZIBBMR : u * S (the expected count given presence, with S = total reads)
## It also returns `is_positive` to restrict to the observations where the
## taxon is present, and the population (only `mu`) and individual (per-subject
## effects, `psi_mean`) versions.
.saem_predict <- function(x) {
  d <- x$data
  if (is.null(d)) {
    stop("This fit did not store the original data (it was created with an ",
         "older version of the package). Refit the model to use these plots.",
         call. = FALSE)
  }
  id <- d$subject_id
  n_subjects <- nrow(x$psi_mean)
  beta_cols <- x$n_alpha + seq_len(x$n_beta)

  # population psi: all rows equal to mu (no per-subject deviations)
  psi_pop <- matrix(x$mu, nrow = n_subjects, ncol = length(x$mu), byrow = TRUE)

  u_ind <- .saem_linear_prob(x$psi_mean, beta_cols, id, d$z_design)
  u_pop <- .saem_linear_prob(psi_pop,    beta_cols, id, d$z_design)

  mult <- if (!is.null(d$S)) d$S else 1  # ZIBBMR: total reads per sample

  list(
    observed    = d$y,
    is_positive = d$y != 0,
    pred_ind    = u_ind * mult,
    pred_pop    = u_pop * mult
  )
}

## Plot 4: observed vs. predicted for the continuous part, using only the
## positive observations (where the taxon is present). Shows the population and
## individual predictions; the red y = x line marks a perfect fit.
.saem_plot_fit <- function(x) {
  pr <- .saem_predict(x)
  pos <- pr$is_positive
  if (!any(pos)) {
    message("There are no positive observations to plot.")
    return(invisible(x))
  }
  obs <- pr$observed[pos]; pi <- pr$pred_ind[pos]; pp <- pr$pred_pop[pos]
  rng <- range(c(obs, pi, pp), na.rm = TRUE)

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)

  graphics::plot(pp, obs, xlim = rng, ylim = rng,
                 xlab = "Predicted (continuous part)", ylab = "Observed",
                 main = "Observed vs. predicted\n(positive observations)",
                 pch = 1, col = grDevices::adjustcolor("gray40", 0.5))
  graphics::points(pi, obs, pch = 19, col = grDevices::adjustcolor("#2166ac", 0.5))
  graphics::abline(0, 1, col = "red", lwd = 2)
  graphics::legend("topleft", bty = "n",
                   pch = c(1, 19, NA), lty = c(NA, NA, 1), lwd = c(NA, NA, 2),
                   col = c("gray40", "#2166ac", "red"),
                   legend = c("population", "individual", "y = x"), cex = 0.85)
  invisible(x)
}

## Plot 5: residuals of the continuous part (observed - individual predicted),
## on the positive observations. Two panels: residuals against the predicted
## value, and their distribution. The red reference marks 0.
.saem_plot_resid <- function(x) {
  pr <- .saem_predict(x)
  pos <- pr$is_positive
  if (!any(pos)) {
    message("There are no positive observations to plot.")
    return(invisible(x))
  }
  pred <- pr$pred_ind[pos]
  resid <- pr$observed[pos] - pred

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::par(mfrow = c(1, 2))

  graphics::plot(pred, resid, xlab = "Predicted (continuous part)",
                 ylab = "Residual (obs - pred)", main = "Residuals vs. predicted",
                 pch = 19, col = grDevices::adjustcolor("black", 0.4))
  graphics::abline(h = 0, col = "red", lwd = 2)

  graphics::hist(resid, main = "Distribution of residuals", xlab = "Residual",
                 col = "#92c5de", border = "white")
  graphics::abline(v = 0, col = "red", lwd = 2)
  invisible(x)
}

## Plot dispatcher used by plot.zibr_saem / plot.zibbmr_saem.
.saem_plot <- function(x, which = c("convergence", "coefficients", "random",
                                    "fit", "residuals"),
                       beta_label = "beta", ...) {
  which <- match.arg(which)
  switch(which,
    convergence  = .saem_plot_trace(x, ...),
    coefficients = .saem_plot_coef(x, beta_label),
    random       = .saem_plot_random(x, beta_label),
    fit          = .saem_plot_fit(x),
    residuals    = .saem_plot_resid(x))
}

.saem_logLik <- function(object) {
  value <- object$loglik
  attr(value, "df") <- length(object$mu) + 1 + length(diag(object$G))
  attr(value, "nobs") <- object$nobs
  class(value) <- "logLik"
  value
}

.saem_coef <- function(object) {
  object$mu
}

.saem_vcov <- function(object) {
  if (is.null(object$fisher_stoch)) {
    stop("The fit does not contain a FIM matrix. Refit with compute_fim = TRUE.", call. = FALSE)
  }

  -solve(object$fisher_stoch)
}

.saem_se <- function(object) {
  sqrt(diag(.saem_vcov(object)))
}
