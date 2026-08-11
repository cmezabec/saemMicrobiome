#' Simular un conjunto de datos de microbioma con varios taxones
#'
#' Genera conteos multinomiales por sujeto-tiempo para varios taxones
#' simultaneamente, pensado como dataset de juguete rapido para probar el
#' flujo completo del paquete (no usa un modelo ZIBR/ZIBBMR generador; para
#' simular datos consistentes con esos modelos usar [simulate_zibr_data()] o
#' [simulate_zibbmr_data()]).
#'
#' @param n_ind Numero de sujetos.
#' @param n_time Numero de observaciones (tiempos) por sujeto.
#' @param n_taxa Numero de taxones a simular.
#' @param N Total de lecturas por observacion (profundidad de secuenciacion).
#' @param seed Semilla aleatoria.
#'
#' @return Una lista con `conteo` (data frame de conteos por taxon, con
#'   columnas `id`, `tiempo`, `grupo`, `N` y una columna por taxon),
#'   `proporcion` (igual pero con las columnas de taxon convertidas a
#'   proporciones) y `taxa` (nombres de las columnas de taxon).
#' @seealso [simulate_zibr_data()], [simulate_zibbmr_data()]
#' @examples
#' sim <- simular_datos_microbioma(n_ind = 5, n_time = 3, n_taxa = 4)
#' head(sim$conteo)
#' @export
simular_datos_microbioma <- function(n_ind = 8, n_time = 3, n_taxa = 5,
                                     N = 1000, seed = 123) {
  set.seed(seed)

  n <- n_ind * n_time
  taxa <- paste0("Taxon", seq_len(n_taxa))

  conteos <- matrix(0, nrow = n, ncol = n_taxa)

  for (i in seq_len(n)) {
    p <- stats::rgamma(n_taxa, shape = 1.2, rate = 1)
    p <- p / sum(p)

    conteos[i, ] <- as.vector(stats::rmultinom(1, size = N, prob = p))
  }

  colnames(conteos) <- taxa

  datos_conteo <- data.frame(
    id = rep(seq_len(n_ind), each = n_time),
    tiempo = rep(seq_len(n_time), times = n_ind),
    grupo = rep(rep(c(0, 1), length.out = n_ind), each = n_time),
    N = N,
    conteos,
    check.names = FALSE
  )

  datos_proporcion <- datos_conteo
  datos_proporcion[, taxa] <- datos_conteo[, taxa, drop = FALSE] /
    rowSums(datos_conteo[, taxa, drop = FALSE])

  list(
    conteo = datos_conteo,
    proporcion = datos_proporcion,
    taxa = taxa
  )
}
