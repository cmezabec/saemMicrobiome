test_that("simular_datos_microbioma devuelve conteos y proporciones consistentes", {
  sim <- simular_datos_microbioma(n_ind = 6, n_time = 3, n_taxa = 4, N = 200, seed = 1)

  expect_named(sim, c("conteo", "proporcion", "taxa"))
  expect_equal(nrow(sim$conteo), 18)
  expect_length(sim$taxa, 4)
  expect_equal(unname(rowSums(sim$conteo[, sim$taxa])), rep(200, 18))
  expect_equal(unname(rowSums(sim$proporcion[, sim$taxa])), rep(1, 18), tolerance = 1e-8)
})

make_romero_fixture <- function() {
  n <- 20
  otu <- data.frame(
    T1 = rep(0, n),
    T2 = rep(c(0, 5), length.out = n),
    T3 = rep(3, n),
    T4 = c(rep(0, 2), rep(4, n - 2))
  )
  sample_data <- data.frame(
    Subect_ID = rep(seq_len(n / 2), each = 2),
    Age = seq(20, 39, length.out = n),
    pregnant = rep(c(0, 1), length.out = n),
    GA_Days = seq(50, 280, length.out = n),
    Total.Read.Counts = rowSums(otu) + 10
  )
  list(SampleData = sample_data, OTU = otu)
}

test_that("prepare_romero_zibr filtra por proporcion de ceros y devuelve abundancias relativas", {
  romero <- make_romero_fixture()
  out <- prepare_romero_zibr(romero, taxa_out = integer(0), zero_range = c(0.1, 0.9))

  expect_true(all(c("T2", "T4") %in% out$taxa))
  expect_false("T1" %in% out$taxa)
  expect_false("T3" %in% out$taxa)
  expect_true(all(out$abundances[[1]] <= 1))
})

test_that("prepare_romero_zibbmr conserva los conteos crudos", {
  romero <- make_romero_fixture()
  out <- prepare_romero_zibbmr(romero, taxa_out = integer(0), zero_range = c(0.1, 0.9))

  expect_true(all(c("T2", "T4") %in% out$taxa))
  expect_equal(out$counts$T2, romero$OTU$T2)
})

test_that("prepare_romero_zibr exige SampleData y OTU", {
  expect_error(
    prepare_romero_zibr(list(SampleData = data.frame())),
    "SampleData y OTU"
  )
})
