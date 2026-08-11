#### Utilidades internas compartidas por ZIBR y ZIBBMR ####
#### Estas funciones no dependen de la verosimilitud de cada modelo: ####
#### construccion de matrices de diseno, predictor lineal logistico, ####
#### replicacion de cadenas MCMC y verificacion de paquetes requeridos. ####

.saem_check_packages <- function(inference = TRUE) {
  required <- "MASS"
  if (inference) {
    required <- c(required, "numDeriv")
  }

  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]

  if (length(missing) > 0) {
    stop(
      "Faltan paquetes requeridos: ",
      paste(missing, collapse = ", "),
      ". Instala antes de ajustar el modelo.",
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

# Inversa de la matriz de covarianza de los efectos aleatorios.
#
# La version anterior calculaba diag(1 / diag(G)), es decir invertia G
# descartando los terminos fuera de la diagonal. Con G diagonal el resultado es
# el mismo, pero el paso M del SAEM estima la matriz COMPLETA (ver el calculo de
# `G_full` en zibbmr.R y zibr.R), de modo que en cuanto los efectos aleatorios
# de las dos partes del modelo estan correlacionados la precision quedaba mal.
#
# Esto afectaba la razon de aceptacion de Metropolis-Hastings, que evaluaba la
# densidad previa como si los efectos aleatorios fueran independientes: las
# trazas salian sin correlacion, la covarianza acumulada no crecia y G volvia a
# ser casi diagonal en la iteracion siguiente. Resultado: un punto fijo
# autoreforzado con rho ~ 0 (con rho verdadero 0.7 se estimaba ~0.09). Tambien
# afectaba la verosimilitud marginal por importance sampling, y con ella el LRT.
.saem_diag_inverse <- function(G) {
  G <- as.matrix(G)

  inversa <- try(solve(G), silent = TRUE)

  # Si G es singular o esta mal condicionada (puede pasar en iteraciones
  # tempranas del SAEM, antes de que las componentes de varianza se estabilicen)
  # se recurre a la pseudo-inversa en vez de fallar.
  if (inherits(inversa, "try-error")) {
    descomposicion <- eigen(G, symmetric = TRUE)
    positivos <- descomposicion$values >
      max(1e-10, max(descomposicion$values) * 1e-10)

    if (!any(positivos)) {
      stop("La matriz de covarianza de los efectos aleatorios es singular.",
           call. = FALSE)
    }

    inversa <- descomposicion$vectors[, positivos, drop = FALSE] %*%
      diag(1 / descomposicion$values[positivos],
           nrow = sum(positivos)) %*%
      t(descomposicion$vectors[, positivos, drop = FALSE])
  }

  inversa
}

# Impone la estructura pedida a la covarianza de los efectos aleatorios.
#
# "diag"          efectos aleatorios independientes entre las dos partes del
#                 modelo. Es la especificacion del articulo y el valor por
#                 defecto: con pocos sujetos u observaciones por sujeto, la
#                 matriz sin restricciones degenera (rho -> +-1 y colapso de una
#                 de las varianzas).
# "unstructured"  covarianza libre. Habilita estimar la correlacion entre el
#                 efecto aleatorio de la parte de inflacion de ceros y el de la
#                 parte de abundancia, que ni glmmTMB ni gamlss pueden
#                 representar. Requiere bastante mas informacion.
.saem_impone_estructura <- function(G, estructura) {
  if (identical(estructura, "diag")) {
    return(.saem_diag(diag(as.matrix(G))))
  }
  as.matrix(G)
}

.saem_valida_estructura <- function(estructura) {
  match.arg(estructura, c("diag", "unstructured"))
}

.saem_covariate_matrix <- function(x, n, prefix) {
  if (is.null(x)) {
    mat <- matrix(nrow = n, ncol = 0)
  } else {
    mat <- as.matrix(x)
    if (nrow(mat) != n) {
      stop("La matriz de covariables no tiene el mismo numero de filas que Y.", call. = FALSE)
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
  # El predictor lineal se calcula en C++ (evita materializar psi[id, cols] y
  # el rowSums); plogis() se aplica en R (vectorizado) para que el resultado
  # sea byte-identico al de la version en R puro.
  eta <- saem_linear_eta_cpp(psi, as.integer(cols), as.integer(id), design)
  plogis(eta)
}


#### Parte de inflacion de ceros: identica para ZIBR y ZIBBMR ####
#### (la probabilidad de presencia/ausencia no depende de si la parte ####
#### positiva es beta o beta-binomial) ####

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


#### Extraccion de log-verosimilitud y comparacion de modelos anidados (LRT) ####
#### Comun a ambos modelos: solo dependen de que el objeto tenga $loglik ####
#### o un metodo logLik() generico. ####

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
    stop("full_models y reduced_models deben tener la misma longitud.", call. = FALSE)
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


#### Metodos S3: la mecanica de imprimir/graficar/extraer coeficientes es ####
#### identica para zibr_saem y zibbmr_saem; solo cambian etiquetas de texto ####

.saem_print <- function(x, model_label, beta_label) {
  cat("===== Resultados ", model_label, " =====\n", sep = "")

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

    cat("== Parte logistica: p_it ==\n")
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

  cat("=== Varianzas de efectos aleatorios ===\n")
  if (x$zi && n_alpha_random > 0) {
    cat("== Parte logistica ==\n")
    print(alpha_tab[x$alpha_random, c("Variance", "sqrt.Var"), drop = FALSE])
  }
  if (n_beta_random > 0) {
    cat("== ", beta_label, " ==\n", sep = "")
    print(beta_tab[x$beta_random, c("Variance", "sqrt.Var"), drop = FALSE])
  }

  cat("=== Phi: ", x$phi, "\n", sep = "")
  cat("=== Log-verosimilitud marginal (importance sampling): ", x$loglik, "\n", sep = "")

  invisible(x)
}

## Etiquetas de los parametros de efectos fijos (parte logistica + parte
## beta/beta-binomial), en el mismo orden que x$mu.
.saem_param_labels <- function(x, beta_label = "beta") {
  labs <- character(0)
  if (isTRUE(x$zi) && x$n_alpha > 0) {
    labs <- paste0("logistica: ", x$alpha_labels)
  }
  c(labs, paste0(beta_label, ": ", x$beta_labels))
}

## Grafico 1: traza de convergencia (parametros a lo largo de las iteraciones).
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
      xlab = "Iteracion",
      ylab = "Valor",
      main = colnames(trace)[j]
    )
    graphics::abline(v = burn_in, lty = 2)
  }

  invisible(x)
}

## Grafico 2: coeficientes estimados con intervalo de confianza al 95%
## (tipo "forest plot"). Necesita el ajuste con compute_fim = TRUE para los IC.
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
                 xlab = "Estimado (IC 95%)", ylab = "", pch = 19,
                 main = "Coeficientes estimados")
  graphics::axis(2, at = seq_len(n), labels = labs, las = 1, cex.axis = 0.85)
  graphics::abline(v = 0, lty = 2, col = "gray50")
  ok <- !is.na(lo)
  if (any(ok)) graphics::segments(lo[ok], which(ok), hi[ok], which(ok), lwd = 2)
  if (!all(ok)) {
    aviso <- if (is.null(x$fisher_stoch)) {
      "Sin IC: reajusta con compute_fim = TRUE"
    } else {
      "IC no disponible para algun coeficiente (matriz de informacion mal condicionada)"
    }
    graphics::mtext(aviso, side = 1, line = 2.5, cex = 0.75, col = "gray40")
  }
  invisible(x)
}

## Grafico 3: distribucion entre sujetos de los efectos aleatorios estimados
## (uno por sujeto). La linea roja marca la media poblacional.
.saem_plot_random <- function(x, beta_label = "beta") {
  ri <- x$random_index
  if (length(ri) == 0) {
    message("El ajuste no tiene efectos aleatorios que graficar.")
    return(invisible(x))
  }
  labs <- .saem_param_labels(x, beta_label)[ri]
  vals <- x$psi_mean[, ri, drop = FALSE]
  k <- ncol(vals)

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::par(mfrow = c(1, k))
  for (j in seq_len(k)) {
    graphics::hist(vals[, j], main = labs[j], xlab = "Valor por sujeto",
                   col = "#92c5de", border = "white")
    graphics::abline(v = x$mu[ri[j]], col = "red", lwd = 2)
  }
  invisible(x)
}

## Predicciones de la PARTE CONTINUA del modelo (la magnitud dado que el taxon
## esta presente), usando los datos originales que el ajuste guarda en `x$data`.
## Se enfoca en la parte continua -no en la marginal E[Y] = p * u- porque en un
## modelo con inflacion de ceros la prediccion marginal mezcla la masa de
## probabilidad en cero con la parte continua y no se lee como un diagrama
## observado-vs-predicho clasico. La media condicional dado presencia es:
##   ZIBR   : u   (la media de la parte beta, en [0, 1))
##   ZIBBMR : u * S (el conteo esperado dado presencia, con S = total de lecturas)
## Devuelve tambien `is_positive` para restringir a las observaciones donde el
## taxon esta presente, y las versiones poblacional (solo `mu`) e individual
## (efectos por sujeto, `psi_mean`).
.saem_predict <- function(x) {
  d <- x$data
  if (is.null(d)) {
    stop("Este ajuste no guardo los datos originales (fue creado con una version ",
         "anterior del paquete). Vuelve a ajustar el modelo para usar estos graficos.",
         call. = FALSE)
  }
  id <- d$subject_id
  n_subjects <- nrow(x$psi_mean)
  beta_cols <- x$n_alpha + seq_len(x$n_beta)

  # psi poblacional: todas las filas iguales a mu (sin desviaciones por sujeto)
  psi_pop <- matrix(x$mu, nrow = n_subjects, ncol = length(x$mu), byrow = TRUE)

  u_ind <- .saem_linear_prob(x$psi_mean, beta_cols, id, d$z_design)
  u_pop <- .saem_linear_prob(psi_pop,    beta_cols, id, d$z_design)

  mult <- if (!is.null(d$S)) d$S else 1  # ZIBBMR: total de lecturas por muestra

  list(
    observed    = d$y,
    is_positive = d$y != 0,
    pred_ind    = u_ind * mult,
    pred_pop    = u_pop * mult
  )
}

## Grafico 4: observados vs. predichos de la parte continua, usando solo las
## observaciones positivas (donde el taxon esta presente). Muestra la prediccion
## poblacional e individual; la recta roja y = x marca el ajuste perfecto.
.saem_plot_fit <- function(x) {
  pr <- .saem_predict(x)
  pos <- pr$is_positive
  if (!any(pos)) {
    message("No hay observaciones positivas que graficar.")
    return(invisible(x))
  }
  obs <- pr$observed[pos]; pi <- pr$pred_ind[pos]; pp <- pr$pred_pop[pos]
  rng <- range(c(obs, pi, pp), na.rm = TRUE)

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)

  graphics::plot(pp, obs, xlim = rng, ylim = rng,
                 xlab = "Predicho (parte continua)", ylab = "Observado",
                 main = "Observados vs. predichos\n(observaciones positivas)",
                 pch = 1, col = grDevices::adjustcolor("gray40", 0.5))
  graphics::points(pi, obs, pch = 19, col = grDevices::adjustcolor("#2166ac", 0.5))
  graphics::abline(0, 1, col = "red", lwd = 2)
  graphics::legend("topleft", bty = "n",
                   pch = c(1, 19, NA), lty = c(NA, NA, 1), lwd = c(NA, NA, 2),
                   col = c("gray40", "#2166ac", "red"),
                   legend = c("poblacional", "individual", "y = x"), cex = 0.85)
  invisible(x)
}

## Grafico 5: residuos de la parte continua (observado - predicho individual),
## en las observaciones positivas. Dos paneles: residuos contra el valor
## predicho, y su distribucion. La referencia roja marca el 0.
.saem_plot_resid <- function(x) {
  pr <- .saem_predict(x)
  pos <- pr$is_positive
  if (!any(pos)) {
    message("No hay observaciones positivas que graficar.")
    return(invisible(x))
  }
  pred <- pr$pred_ind[pos]
  resid <- pr$observed[pos] - pred

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::par(mfrow = c(1, 2))

  graphics::plot(pred, resid, xlab = "Predicho (parte continua)",
                 ylab = "Residuo (obs - pred)", main = "Residuos vs. predicho",
                 pch = 19, col = grDevices::adjustcolor("black", 0.4))
  graphics::abline(h = 0, col = "red", lwd = 2)

  graphics::hist(resid, main = "Distribucion de residuos", xlab = "Residuo",
                 col = "#92c5de", border = "white")
  graphics::abline(v = 0, col = "red", lwd = 2)
  invisible(x)
}

## Despachador de graficos usado por plot.zibr_saem / plot.zibbmr_saem.
.saem_plot <- function(x, which = c("convergencia", "coeficientes", "aleatorios",
                                    "ajuste", "residuos"),
                       beta_label = "beta", ...) {
  which <- match.arg(which)
  switch(which,
    convergencia = .saem_plot_trace(x, ...),
    coeficientes = .saem_plot_coef(x, beta_label),
    aleatorios   = .saem_plot_random(x, beta_label),
    ajuste       = .saem_plot_fit(x),
    residuos     = .saem_plot_resid(x))
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
    stop("El ajuste no contiene matriz FIM. Reajusta con compute_fim = TRUE.", call. = FALSE)
  }

  -solve(object$fisher_stoch)
}

.saem_se <- function(object) {
  sqrt(diag(.saem_vcov(object)))
}
