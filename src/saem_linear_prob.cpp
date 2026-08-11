#include <Rcpp.h>
using namespace Rcpp;

// Version compilada de .saem_linear_prob(): calcula el predictor lineal
// logistico por observacion y le aplica plogis(). Es la funcion mas llamada
// del motor SAEM (~6 veces por iteracion). Fusiona en un solo recorrido el
// "gather" de filas de psi por sujeto-cadena (id), el producto por la matriz
// de diseno y la suma por fila, evitando materializar las matrices
// intermedias que creaba la version en R puro.
//
// Es deterministica (no usa RNG). Devuelve SOLO el predictor lineal eta; el
// logistico plogis(eta) se aplica despues en R (vectorizado), de modo que el
// resultado es byte-identico al de la version en R puro (misma rutina
// stats::plogis()) y a la vez se evita materializar en R la matriz
// intermedia psi[id, cols] y el rowSums.
//
// psi:    matriz (n_subjects*n_chains) x n_psi con los efectos por sujeto-cadena
// cols:   indices de columna de psi a usar (base 1, como en R)
// id:     indice de fila de psi para cada observacion (base 1, longitud M)
// design: matriz M x length(cols) con las covariables (incluye intercepto)
//
// [[Rcpp::export]]
NumericVector saem_linear_eta_cpp(NumericMatrix psi,
                                  IntegerVector cols,
                                  IntegerVector id,
                                  NumericMatrix design) {
  const int m = design.nrow();
  const int k = cols.size();
  NumericVector eta(m);

  // indices base 0
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
