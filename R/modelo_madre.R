#' Ajustar un modelo SAEM de microbioma (ZIBR o ZIBBMR) desde un data frame
#'
#' Funcion "madre" que despacha a [saem_zibr_clean()] o [saem_zibbmr_clean()]
#' segun `modelo`, extrayendo la respuesta, el total de lecturas (si aplica)
#' y las covariables directamente de un data frame por nombre de columna, y
#' generando valores iniciales por defecto razonables si no se entregan.
#' Pensada para ajustar rapidamente un taxon sin tener que armar a mano las
#' matrices `X`/`Z` y los vectores de valores iniciales que piden
#' [fit_zibr()]/[fit_zibbmr()]. Sigue el mismo patron que
#' [fit_zibr_taxon()]/[fit_zibbmr_taxon()].
#'
#' `fit_saem_microbiome()` es un alias identico, con un nombre que sigue la
#' convencion `fit_*` del resto del paquete.
#'
#' @param modelo Cual modelo ajustar: `"zibbmr"` (por defecto, conteos con
#'   profundidad de secuenciacion) o `"zibr"` (proporciones).
#' @param datos Data frame con una fila por observacion.
#' @param taxon Nombre de la columna en `datos` con la respuesta (conteo o
#'   proporcion, segun `modelo`) a modelar.
#' @param id Nombre de la columna en `datos` que identifica al sujeto.
#' @param total Nombre de la columna en `datos` con el total de lecturas.
#'   Solo se usa si `modelo = "zibbmr"`.
#' @param covariables Vector de nombres de columnas en `datos` a usar como
#'   covariables tanto en la parte logistica como en la parte
#'   beta/beta-binomial. Se ignora si se entregan `x_covariables`/
#'   `z_covariables` por separado.
#' @param x_covariables Nombres de columnas para la parte de inflacion de
#'   ceros (por defecto, igual a `covariables`).
#' @param z_covariables Nombres de columnas para la parte beta/beta-binomial
#'   (por defecto, igual a `covariables`).
#' @param zi Logico, ver [fit_zibr()]/[fit_zibbmr()].
#' @param phi_start Valor inicial de dispersion. Si es `NULL` (por defecto),
#'   se sortea con `runif(1, 10, 20)`.
#' @param alpha_start Valores iniciales de la parte logistica. Si es `NULL`
#'   (por defecto), se sortean con `runif(., -0.1, 0.1)`. Se ignora si
#'   `zi = FALSE`.
#' @param beta_start Valores iniciales de la parte beta/beta-binomial. Si es
#'   `NULL` (por defecto), se sortean con `runif(., -0.1, 0.1)`.
#' @param iter Numero de iteraciones SAEM.
#' @param ncad Numero de cadenas MCMC.
#' @param compute_fim Logico. Si `TRUE`, calcula la matriz de informacion de
#'   Fisher estocastica (necesaria para `vcov()`/`se()`).
#' @param seed Semilla aleatoria, usada tanto para los valores iniciales
#'   (cuando se sortean) como para el ajuste SAEM.
#'
#' @return Un objeto `zibr_saem` o `zibbmr_saem`, segun `modelo`.
#' @seealso [fit_zibr_taxon()], [fit_zibbmr_taxon()] para envoltorios
#'   equivalentes con la firma de argumentos moderna (`n_iter`, `n_chains`).
#' @export
ajustar_modelo_microbioma <- function(modelo = c("zibbmr", "zibr"),
                                      datos,
                                      taxon,
                                      id = "id",
                                      total = "N",
                                      covariables = c("tiempo", "grupo"),
                                      x_covariables = covariables,
                                      z_covariables = covariables,
                                      zi = TRUE,
                                      phi_start = NULL,
                                      alpha_start = NULL,
                                      beta_start = NULL,
                                      iter = 200,
                                      ncad = 5,
                                      compute_fim = FALSE,
                                      seed = 1) {
  modelo <- match.arg(modelo)

  if (!taxon %in% names(datos)) {
    stop("El taxon indicado no existe en datos.", call. = FALSE)
  }
  if (!id %in% names(datos)) {
    stop("La columna id no existe en datos.", call. = FALSE)
  }
  if (modelo == "zibbmr" && !total %in% names(datos)) {
    stop("La columna total no existe en datos.", call. = FALSE)
  }

  if (is.null(x_covariables)) {
    x_covariables <- character(0)
  }
  if (is.null(z_covariables)) {
    z_covariables <- character(0)
  }

  if (!all(x_covariables %in% names(datos))) {
    stop("Al menos una covariable de x_covariables no existe en datos.", call. = FALSE)
  }
  if (!all(z_covariables %in% names(datos))) {
    stop("Al menos una covariable de z_covariables no existe en datos.", call. = FALSE)
  }

  n_x_covariables <- length(x_covariables)
  n_z_covariables <- length(z_covariables)

  if (!is.null(seed)) {
    set.seed(seed)
  }

  if (is.null(phi_start)) {
    phi_start <- stats::runif(1, 10, 20)
  }
  if (zi && is.null(alpha_start)) {
    alpha_start <- stats::runif(n_x_covariables + 1, -0.1, 0.1)
  }
  if (is.null(beta_start)) {
    beta_start <- stats::runif(n_z_covariables + 1, -0.1, 0.1)
  }

  Y <- datos[[taxon]]
  index <- datos[[id]]

  X <- if (!zi || n_x_covariables == 0) NULL else as.matrix(datos[, x_covariables, drop = FALSE])
  Z <- if (n_z_covariables == 0) NULL else as.matrix(datos[, z_covariables, drop = FALSE])

  if (modelo == "zibr") {
    ajuste <- saem_zibr_clean(
      Y = Y,
      X = X,
      Z = Z,
      index = index,
      zi = zi,
      v0 = phi_start,
      a0 = alpha_start,
      b0 = beta_start,
      iter = iter,
      ncad = ncad,
      seed = seed,
      compute_fim = compute_fim
    )
  }

  if (modelo == "zibbmr") {
    ajuste <- saem_zibbmr_clean(
      Y = Y,
      S = datos[[total]],
      X = X,
      Z = Z,
      index = index,
      zi = zi,
      v0 = phi_start,
      a0 = alpha_start,
      b0 = beta_start,
      iter = iter,
      ncad = ncad,
      seed = seed,
      compute_fim = compute_fim
    )
  }

  ajuste
}

#' @rdname ajustar_modelo_microbioma
#' @export
fit_saem_microbiome <- function(modelo = c("zibbmr", "zibr"),
                                datos,
                                taxon,
                                id = "id",
                                total = "N",
                                covariables = c("tiempo", "grupo"),
                                x_covariables = covariables,
                                z_covariables = covariables,
                                zi = TRUE,
                                phi_start = NULL,
                                alpha_start = NULL,
                                beta_start = NULL,
                                iter = 200,
                                ncad = 5,
                                compute_fim = FALSE,
                                seed = 1) {
  ajustar_modelo_microbioma(
    modelo = modelo,
    datos = datos,
    taxon = taxon,
    id = id,
    total = total,
    x_covariables = x_covariables,
    z_covariables = z_covariables,
    zi = zi,
    phi_start = phi_start,
    alpha_start = alpha_start,
    beta_start = beta_start,
    iter = iter,
    ncad = ncad,
    compute_fim = compute_fim,
    seed = seed
  )
}
