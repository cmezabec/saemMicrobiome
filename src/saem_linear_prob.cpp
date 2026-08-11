#include <Rcpp.h>
using namespace Rcpp;

// Compiled version of .saem_linear_prob(): computes the logistic linear
// predictor per observation and applies plogis() to it. It is the most called
// function of the SAEM engine (~6 times per iteration). In a single pass it
// fuses the "gather" of psi rows by subject-chain (id), the product with the
// design matrix and the row-wise sum, avoiding materializing the intermediate
// matrices that the pure-R version created.
//
// It is deterministic (no RNG). It returns ONLY the linear predictor eta; the
// logistic plogis(eta) is applied afterwards in R (vectorized), so the result
// is byte-identical to the pure-R version (same stats::plogis() routine) while
// also avoiding materializing in R the intermediate matrix psi[id, cols] and
// the rowSums.
//
// psi:    (n_subjects*n_chains) x n_psi matrix with the per-subject-chain effects
// cols:   column indices of psi to use (1-based, as in R)
// id:     row index of psi for each observation (1-based, length M)
// design: M x length(cols) matrix with the covariates (includes intercept)
//
// [[Rcpp::export]]
NumericVector saem_linear_eta_cpp(NumericMatrix psi,
                                  IntegerVector cols,
                                  IntegerVector id,
                                  NumericMatrix design) {
  const int m = design.nrow();
  const int k = cols.size();
  NumericVector eta(m);

  // 0-based indices
  std::vector<int> col0(k);
  for (int j = 0; j < k; ++j) col0[j] = cols[j] - 1;

  for (int i = 0; i < m; ++i) {
    const int row = id[i] - 1;
    double s = 0.0;
    for (int j = 0; j < k; ++j) {
      s += psi(row, col0[j]) * design(i, j);
    }
    eta[i] = s;
  }

  return eta;
}
