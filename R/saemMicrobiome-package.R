#' saemMicrobiome: SAEM models for longitudinal microbiome data
#'
#' Tools to fit, simulate and compare zero-inflated mixed models for
#' longitudinal microbiome data:
#'
#' - **ZIBR** (zero-inflated beta regression), for proportions or relative
#'   abundances: see [fit_zibr()].
#' - **ZIBBMR** (zero-inflated beta-binomial mixed regression), for counts
#'   with known sequencing depth: see [fit_zibbmr()].
#'
#' Both are estimated with the Stochastic Approximation EM (SAEM) algorithm,
#' following the methodology developed by John Barrera.
#'
#' @keywords internal
#' @importFrom stats rnorm runif rt dt plogis rbinom vcov
#' @useDynLib saemMicrobiome, .registration = TRUE
#' @importFrom Rcpp sourceCpp
"_PACKAGE"
