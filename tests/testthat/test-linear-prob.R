test_that(".saem_linear_prob never returns exactly 0 or 1", {
  # Linear predictors far enough out for plogis() to round to 0 or 1, as when a
  # fixed effect diverges during the burn-in on data with very few positives.
  eta <- c(-800, -30, 0, 30, 800)
  p <- saemMicrobiome:::.saem_linear_prob(matrix(eta, ncol = 1), 1L, seq_along(eta),
                                          matrix(1, length(eta), 1))
  expect_true(all(p > 0 & p < 1))
  expect_true(all(is.finite(lgamma(10 * p)) & is.finite(lgamma(10 * (1 - p)))))
  # Values that plogis() leaves strictly inside (0, 1) are untouched
  expect_identical(p[2:4], plogis(eta[2:4]))
})
