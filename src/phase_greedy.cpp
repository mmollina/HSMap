// src/phase_greedy.cpp
// Greedy coordinate ascent for the signed-graph phase objective used by
// phase_from_pairwise(). This is an internal performance replacement for the R
// reference `.pf_greedy()`; it implements the identical algorithm (same objective,
// same deterministic best-improvement rule and tie-breaking) with incremental score
// updates, so it converges to the same local optimum without the R-loop overhead.
//
// Objective maximized over x in {+1,-1}^n:  sum_{i<j} J_ij x_i x_j.
// Per pass: g = J x; then repeatedly flip the marker with the largest positive gain
// delta_i = -2 x_i g_i (first index on ties, matching R's which.max), updating g
// incrementally (g -= 2 * x_old * J[,best]) after each accepted flip.
//
// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
using namespace Rcpp;

//' Greedy coordinate ascent for the signed phase graph (internal)
//' @keywords internal
// [[Rcpp::export]]
List pf_greedy_cpp(NumericMatrix J, IntegerVector x_init, int max_passes, double tol) {
  const int n = x_init.size();
  std::vector<double> x(n), g(n, 0.0);
  for (int i = 0; i < n; ++i) x[i] = (double) x_init[i];   // +1 / -1

  int nflips = 0, pass = 0;
  bool improve = true;

  // g = J x, computed ONCE (column-major friendly) and then maintained incrementally:
  // after each accepted flip of column `best`, g is updated in O(n) rather than being
  // recomputed, and the full objective is never recomputed inside the loop.
  for (int j = 0; j < n; ++j) {
    const double xj = x[j];
    if (xj == 0.0) continue;
    for (int i = 0; i < n; ++i) g[i] += J(i, j) * xj;
  }

  while (improve && pass < max_passes) {
    ++pass;
    improve = false;

    const int inner_cap = 5 * n;
    for (int it = 0; it < inner_cap; ++it) {
      // best = argmax_i (-2 x_i g_i), first index on ties (matches which.max)
      int best = 0;
      double bestval = -2.0 * x[0] * g[0];
      for (int i = 1; i < n; ++i) {
        const double d = -2.0 * x[i] * g[i];
        if (d > bestval) { bestval = d; best = i; }
      }
      if (bestval <= tol) break;

      improve = true; ++nflips;
      const double xo = x[best];
      x[best] = -xo;
      // g -= 2 * xo * J[,best]  (column best is contiguous)
      const double f = 2.0 * xo;
      for (int i = 0; i < n; ++i) g[i] -= f * J(i, best);
    }
  }

  // objective = sum_{i<j} J_ij x_i x_j
  double obj = 0.0;
  for (int i = 0; i < n; ++i) {
    const double xi = x[i];
    for (int j = i + 1; j < n; ++j) obj += J(i, j) * xi * x[j];
  }

  IntegerVector xout(n);
  for (int i = 0; i < n; ++i) xout[i] = (int) x[i];
  return List::create(_["x"] = xout, _["objective"] = obj,
                      _["n_flips"] = nflips, _["iters"] = pass,
                      _["converged"] = !improve);
}
