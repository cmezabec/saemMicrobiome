#' saemMicrobiome: modelos SAEM para microbioma longitudinal
#'
#' Herramientas para ajustar, simular y comparar modelos mixtos con
#' inflacion de ceros para datos longitudinales de microbioma:
#'
#' - **ZIBR** (zero-inflated beta regression), para proporciones o
#'   abundancias relativas: ver [fit_zibr()].
#' - **ZIBBMR** (zero-inflated beta-binomial mixed regression), para
#'   conteos con profundidad de secuenciacion conocida: ver [fit_zibbmr()].
#'
#' Ambos se estiman con el algoritmo Stochastic Approximation EM (SAEM),
#' siguiendo la metodologia desarrollada por John Barrera.
#'
#' @keywords internal
#' @importFrom stats rnorm runif rt dt plogis rbinom vcov
#' @useDynLib saemMicrobiome, .registration = TRUE
#' @importFrom Rcpp sourceCpp
"_PACKAGE"
