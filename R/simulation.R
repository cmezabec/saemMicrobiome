#' Simulate a microbiome data set with several taxa
#'
#' Generates multinomial counts per subject-time for several taxa
#' simultaneously, meant as a quick toy data set to test the full package
#' workflow (it does not use a ZIBR/ZIBBMR generating model; to simulate data
#' consistent with those models use [simulate_zibr_data()] or
#' [simulate_zibbmr_data()]).
#'
#' @param n_ind Number of subjects.
#' @param n_time Number of observations (time points) per subject.
#' @param n_taxa Number of taxa to simulate.
#' @param N Total reads per observation (sequencing depth).
#' @param seed Random seed.
#'
#' @return A list with `count` (data frame of counts per taxon, with columns
#'   `id`, `time`, `group`, `N` and one column per taxon), `proportion`
#'   (the same but with the taxon columns converted to proportions) and
#'   `taxa` (names of the taxon columns).
#' @seealso [simulate_zibr_data()], [simulate_zibbmr_data()]
#' @examples
#' sim <- simulate_microbiome_data(n_ind = 5, n_time = 3, n_taxa = 4)
#' head(sim$count)
#' @export
simulate_microbiome_data <- function(n_ind = 8, n_time = 3, n_taxa = 5,
                                     N = 1000, seed = 123) {
  set.seed(seed)

  n <- n_ind * n_time
  taxa <- paste0("Taxon", seq_len(n_taxa))

  counts <- matrix(0, nrow = n, ncol = n_taxa)

  for (i in seq_len(n)) {
    p <- stats::rgamma(n_taxa, shape = 1.2, rate = 1)
    p <- p / sum(p)

    counts[i, ] <- as.vector(stats::rmultinom(1, size = N, prob = p))
  }

  colnames(counts) <- taxa

  count_data <- data.frame(
    id = rep(seq_len(n_ind), each = n_time),
    time = rep(seq_len(n_time), times = n_ind),
    group = rep(rep(c(0, 1), length.out = n_ind), each = n_time),
    N = N,
    counts,
    check.names = FALSE
  )

  proportion_data <- count_data
  proportion_data[, taxa] <- count_data[, taxa, drop = FALSE] /
    rowSums(count_data[, taxa, drop = FALSE])

  list(
    count = count_data,
    proportion = proportion_data,
    taxa = taxa
  )
}
