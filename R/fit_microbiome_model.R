#' Fit a SAEM microbiome model (ZIBR or ZIBBMR) from a data frame
#'
#' "Mother" function that dispatches to [saem_zibr_clean()] or
#' [saem_zibbmr_clean()] depending on `model`, extracting the response, the
#' total reads (if applicable) and the covariates directly from a data frame
#' by column name, and generating reasonable default starting values if none
#' are supplied. Meant to quickly fit a single taxon without having to build
#' the `X`/`Z` matrices and the starting-value vectors that
#' [fit_zibr()]/[fit_zibbmr()] require by hand. Follows the same pattern as
#' [fit_zibr_taxon()]/[fit_zibbmr_taxon()].
#'
#' `fit_saem_microbiome()` is an identical alias, with a name that follows the
#' `fit_*` convention used across the package.
#'
#' @param model Which model to fit: `"zibbmr"` (default, counts with
#'   sequencing depth) or `"zibr"` (proportions).
#' @param data Data frame with one row per observation.
#' @param taxon Name of the column in `data` with the response (count or
#'   proportion, depending on `model`) to be modelled.
#' @param id Name of the column in `data` that identifies the subject.
#' @param total Name of the column in `data` with the total reads. Only used
#'   if `model = "zibbmr"`.
#' @param covariates Vector of column names in `data` to use as covariates in
#'   both the logistic part and the beta/beta-binomial part. Ignored if
#'   `x_covariates`/`z_covariates` are supplied separately.
#' @param x_covariates Column names for the zero-inflation part (by default,
#'   equal to `covariates`).
#' @param z_covariates Column names for the beta/beta-binomial part (by
#'   default, equal to `covariates`).
#' @param zi Logical, see [fit_zibr()]/[fit_zibbmr()].
#' @param phi_start Starting value for the dispersion. If `NULL` (default), it
#'   is drawn with `runif(1, 10, 20)`.
#' @param alpha_start Starting values for the logistic part. If `NULL`
#'   (default), they are drawn with `runif(., -0.1, 0.1)`. Ignored if
#'   `zi = FALSE`.
#' @param beta_start Starting values for the beta/beta-binomial part. If
#'   `NULL` (default), they are drawn with `runif(., -0.1, 0.1)`.
#' @param iter Number of SAEM iterations.
#' @param ncad Number of MCMC chains.
#' @param compute_fim Logical. If `TRUE`, computes the stochastic Fisher
#'   information matrix (needed for `vcov()`/`se()`).
#' @param seed Random seed, used both for the starting values (when drawn) and
#'   for the SAEM fit.
#'
#' @return A `zibr_saem` or `zibbmr_saem` object, depending on `model`.
#' @seealso [fit_zibr_taxon()], [fit_zibbmr_taxon()] for equivalent wrappers
#'   with the modern argument signature (`n_iter`, `n_chains`).
#' @export
fit_microbiome_model <- function(model = c("zibbmr", "zibr"),
                                 data,
                                 taxon,
                                 id = "id",
                                 total = "N",
                                 covariates = c("time", "group"),
                                 x_covariates = covariates,
                                 z_covariates = covariates,
                                 zi = TRUE,
                                 phi_start = NULL,
                                 alpha_start = NULL,
                                 beta_start = NULL,
                                 iter = 200,
                                 ncad = 5,
                                 compute_fim = FALSE,
                                 seed = 1) {
  model <- match.arg(model)

  if (!taxon %in% names(data)) {
    stop("The specified taxon does not exist in data.", call. = FALSE)
  }
  if (!id %in% names(data)) {
    stop("The id column does not exist in data.", call. = FALSE)
  }
  if (model == "zibbmr" && !total %in% names(data)) {
    stop("The total column does not exist in data.", call. = FALSE)
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
    phi_start <- stats::runif(1, 10, 20)
  }
  if (zi && is.null(alpha_start)) {
    alpha_start <- stats::runif(n_x_covariates + 1, -0.1, 0.1)
  }
  if (is.null(beta_start)) {
    beta_start <- stats::runif(n_z_covariates + 1, -0.1, 0.1)
  }

  Y <- data[[taxon]]
  index <- data[[id]]

  X <- if (!zi || n_x_covariates == 0) NULL else as.matrix(data[, x_covariates, drop = FALSE])
  Z <- if (n_z_covariates == 0) NULL else as.matrix(data[, z_covariates, drop = FALSE])

  if (model == "zibr") {
    fit <- saem_zibr_clean(
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

  if (model == "zibbmr") {
    fit <- saem_zibbmr_clean(
      Y = Y,
      S = data[[total]],
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

  fit
}

#' @rdname fit_microbiome_model
#' @export
fit_saem_microbiome <- function(model = c("zibbmr", "zibr"),
                                data,
                                taxon,
                                id = "id",
                                total = "N",
                                covariates = c("time", "group"),
                                x_covariates = covariates,
                                z_covariates = covariates,
                                zi = TRUE,
                                phi_start = NULL,
                                alpha_start = NULL,
                                beta_start = NULL,
                                iter = 200,
                                ncad = 5,
                                compute_fim = FALSE,
                                seed = 1) {
  fit_microbiome_model(
    model = model,
    data = data,
    taxon = taxon,
    id = id,
    total = total,
    x_covariates = x_covariates,
    z_covariates = z_covariates,
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
