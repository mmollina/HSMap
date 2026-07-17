// src/hmm_fullsib.cpp
// Oracle-phase four-state full-sib HMM (known genotyped mother AND sire).
// See dev/known_sire_design.md for the model.
//
// Hidden state Z = (maternal homolog hm, paternal homolog hp), hm,hp in {0,1};
// encoded s = 2*hm + hp in {0,1,2,3}. Initial 1/4 each. Transition is the Kronecker
// product T = Tm(rm) (x) Tp(rp) (phase-free). Phase enters ONLY the emission through
// the phased homolog allele labels Am (2xT) and Ap (2xT): dosage d = Am[hm,k]+Ap[hp,k];
// P(y|s) = 1-eps if y==d, eps/2 otherwise, 1 if y missing. No paternal q.
//
// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <cmath>
#include <vector>
using namespace Rcpp;

static inline double fs_emit(int y, int d, double eps) {
  if (y == NA_INTEGER) return 1.0;
  return (y == d) ? (1.0 - eps) : (0.5 * eps);
}

// Fill the 4-vector of emissions at marker k for one offspring genotype row.
static inline void fs_emis_col(int y, const IntegerMatrix& Am, const IntegerMatrix& Ap,
                               int k, double eps, double e[4]) {
  for (int hm = 0; hm < 2; ++hm) {
    int am = Am(hm, k);
    for (int hp = 0; hp < 2; ++hp) {
      int d = am + Ap(hp, k);
      e[2 * hm + hp] = fs_emit(y, d, eps);
    }
  }
}

// Core forward-backward for one offspring. Returns the observed-data log-likelihood
// contribution; if counts != nullptr, accumulates maternal/paternal switch and total
// counts per interval; if gamma_out != nullptr (length 4T), writes marginal posteriors.
static double fs_fb_one(const IntegerVector& y, const IntegerMatrix& Am,
                        const IntegerMatrix& Ap, const NumericVector& rm,
                        const NumericVector& rp, double eps, int T,
                        double* m_sw, double* p_sw, double* tot,
                        double* gamma_out) {
  std::vector<double> alpha(4 * T), E(4 * T), c(T);
  // emissions
  for (int k = 0; k < T; ++k) { double e[4]; fs_emis_col(y[k], Am, Ap, k, eps, e);
    for (int s = 0; s < 4; ++s) E[4 * k + s] = e[s]; }
  // forward (scaled)
  double c0 = 0.0;
  for (int s = 0; s < 4; ++s) { alpha[s] = 0.25 * E[s]; c0 += alpha[s]; }
  if (c0 <= 0) c0 = 1e-300;
  for (int s = 0; s < 4; ++s) alpha[s] /= c0;
  c[0] = c0;
  double loglik = std::log(c0);
  for (int k = 1; k < T; ++k) {
    double rmv = rm[k - 1], rpv = rp[k - 1];
    double ck = 0.0;
    for (int sp = 0; sp < 4; ++sp) {
      int hmp = sp >> 1, hpp = sp & 1;
      double acc = 0.0;
      for (int s = 0; s < 4; ++s) {
        int hm = s >> 1, hp = s & 1;
        double tm = (hm == hmp) ? (1.0 - rmv) : rmv;
        double tp = (hp == hpp) ? (1.0 - rpv) : rpv;
        acc += alpha[4 * (k - 1) + s] * tm * tp;
      }
      alpha[4 * k + sp] = acc * E[4 * k + sp];
      ck += alpha[4 * k + sp];
    }
    if (ck <= 0) ck = 1e-300;
    for (int sp = 0; sp < 4; ++sp) alpha[4 * k + sp] /= ck;
    c[k] = ck;
    loglik += std::log(ck);
  }
  if (m_sw == nullptr && gamma_out == nullptr) return loglik;

  // backward (scaled)
  std::vector<double> beta(4 * T);
  for (int s = 0; s < 4; ++s) beta[4 * (T - 1) + s] = 1.0;
  for (int k = T - 2; k >= 0; --k) {
    double rmv = rm[k], rpv = rp[k];
    for (int s = 0; s < 4; ++s) {
      int hm = s >> 1, hp = s & 1;
      double acc = 0.0;
      for (int sp = 0; sp < 4; ++sp) {
        int hmp = sp >> 1, hpp = sp & 1;
        double tm = (hm == hmp) ? (1.0 - rmv) : rmv;
        double tp = (hp == hpp) ? (1.0 - rpv) : rpv;
        acc += tm * tp * E[4 * (k + 1) + sp] * beta[4 * (k + 1) + sp];
      }
      beta[4 * k + s] = acc / c[k + 1];
    }
  }
  // gamma (marginal posteriors) if requested
  if (gamma_out != nullptr) {
    for (int k = 0; k < T; ++k) {
      double z = 0.0; double g[4];
      for (int s = 0; s < 4; ++s) { g[s] = alpha[4 * k + s] * beta[4 * k + s]; z += g[s]; }
      if (z <= 0) z = 1e-300;
      for (int s = 0; s < 4; ++s) gamma_out[4 * k + s] = g[s] / z;
    }
  }
  // xi -> switch/total counts
  if (m_sw != nullptr) {
    for (int k = 0; k < T - 1; ++k) {
      double rmv = rm[k], rpv = rp[k];
      for (int s = 0; s < 4; ++s) {
        int hm = s >> 1, hp = s & 1;
        double a = alpha[4 * k + s];
        for (int sp = 0; sp < 4; ++sp) {
          int hmp = sp >> 1, hpp = sp & 1;
          double tm = (hm == hmp) ? (1.0 - rmv) : rmv;
          double tp = (hp == hpp) ? (1.0 - rpv) : rpv;
          double xi = a * tm * tp * E[4 * (k + 1) + sp] * beta[4 * (k + 1) + sp] / c[k + 1];
          tot[k] += xi;
          if (hm != hmp) m_sw[k] += xi;
          if (hp != hpp) p_sw[k] += xi;
        }
      }
    }
  }
  return loglik;
}

//' Full-sib four-state forward log-likelihood (summed over offspring)
//' @keywords internal
// [[Rcpp::export]]
double fs_loglik_cpp(IntegerMatrix G, IntegerMatrix Am, IntegerMatrix Ap,
                     NumericVector rm, NumericVector rp, double epsilon) {
  int n = G.nrow(), T = G.ncol();
  double ll = 0.0;
  for (int i = 0; i < n; ++i) {
    IntegerVector y = G(i, _);
    ll += fs_fb_one(y, Am, Ap, rm, rp, epsilon, T, nullptr, nullptr, nullptr, nullptr);
  }
  return ll;
}

//' Full-sib E-step: observed log-likelihood and expected maternal/paternal switch
//' and total counts per interval.
//' @keywords internal
// [[Rcpp::export]]
List fs_estep_cpp(IntegerMatrix G, IntegerMatrix Am, IntegerMatrix Ap,
                  NumericVector rm, NumericVector rp, double epsilon,
                  bool return_gamma = false) {
  int n = G.nrow(), T = G.ncol();
  std::vector<double> m_sw(T - 1, 0.0), p_sw(T - 1, 0.0), tot(T - 1, 0.0);
  double ll = 0.0;
  NumericMatrix gamma(return_gamma ? n : 0, return_gamma ? 4 * T : 0);
  for (int i = 0; i < n; ++i) {
    IntegerVector y = G(i, _);
    std::vector<double> g;
    double* gptr = nullptr;
    if (return_gamma) { g.assign(4 * T, 0.0); gptr = g.data(); }
    ll += fs_fb_one(y, Am, Ap, rm, rp, epsilon, T,
                    m_sw.data(), p_sw.data(), tot.data(), gptr);
    if (return_gamma) for (int j = 0; j < 4 * T; ++j) gamma(i, j) = g[j];
  }
  NumericVector m_switch(T - 1), p_switch(T - 1), total(T - 1);
  for (int k = 0; k < T - 1; ++k) { m_switch[k] = m_sw[k]; p_switch[k] = p_sw[k]; total[k] = tot[k]; }
  List out = List::create(_["loglik"] = ll, _["m_switch"] = m_switch,
                          _["p_switch"] = p_switch, _["total"] = total);
  if (return_gamma) out["gamma"] = gamma;
  return out;
}
