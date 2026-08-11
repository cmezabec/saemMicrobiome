#### Package-specific S3 generics ####

#' Standard errors of a SAEM fit
#'
#' Generic that computes the standard errors of the parameters of a fitted
#' model, from the stochastic Fisher information matrix (`fisher_stoch`)
#' stored in the fit object. Requires the model to have been fitted with
#' `compute_fim = TRUE`.
#'
#' @param object An object fitted by [fit_zibr()] or [fit_zibbmr()] (classes
#'   `zibr_saem` or `zibbmr_saem`).
#' @param ... Additional arguments passed to specific methods.
#'
#' @return A numeric vector with the standard errors, in the same order as
#'   `coef(object)` followed by the dispersion and the random-effect
#'   variances.
#'
#' @seealso [fit_zibr()], [fit_zibbmr()], [stats::vcov()]
#' @export
se <- function(object, ...) {
  UseMethod("se")
}
