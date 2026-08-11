#### Genericos S3 propios del paquete ####

#' Errores estandar de un ajuste SAEM
#'
#' Generico que calcula los errores estandar de los parametros de un modelo
#' ajustado, a partir de la matriz de informacion de Fisher estocastica
#' (`fisher_stoch`) almacenada en el objeto de ajuste. Requiere que el modelo
#' se haya ajustado con `compute_fim = TRUE`.
#'
#' @param object Un objeto ajustado por [fit_zibr()] o [fit_zibbmr()] (clases
#'   `zibr_saem` o `zibbmr_saem`).
#' @param ... Argumentos adicionales pasados a metodos especificos.
#'
#' @return Un vector numerico con los errores estandar, en el mismo orden que
#'   `coef(object)` seguido de la dispersion y las varianzas de efectos
#'   aleatorios.
#'
#' @seealso [fit_zibr()], [fit_zibbmr()], [stats::vcov()]
#' @export
se <- function(object, ...) {
  UseMethod("se")
}
