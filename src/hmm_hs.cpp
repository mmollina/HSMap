// src/hmm_hs.cpp
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppParallel)]]

#include <Rcpp.h>
#include <RcppParallel.h>
#include <array>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <string>

using namespace Rcpp;
using namespace RcppParallel;

// ------------------------------------------------------------------
// Utilities
// ------------------------------------------------------------------

// small utilities
inline double clamp01(double x, double lo = 1e-9, double hi = 1.0 - 1e-9) {
  return x < lo ? lo : (x > hi ? hi : x);
}

// P(y | maternal state h in {0,1}, paternal g in {1..3}=AA,Aa,aa,
// maternal genotype Mt in {0,1,2}=aa,Aa,AA) with error epsilon
inline double py_hgM_obs(int y, int h, int g, int Mt, double epsilon) {
  if (y == NA_INTEGER) return 1.0; // missing observation
  double p = 0.0;
  if (Mt == 1) { // maternal Aa
    if (h == 1) { // mom transmits 'A'
      if      (y == 2) p = (g==1 ? 1.0 : (g==2 ? 0.5 : 0.0));
      else if (y == 1) p = (g==3 ? 1.0 : (g==2 ? 0.5 : 0.0));
      else              p = 0.0;
    } else {       // mom transmits 'a'
      if      (y == 0) p = (g==3 ? 1.0 : (g==2 ? 0.5 : 0.0));
      else if (y == 1) p = (g==1 ? 1.0 : (g==2 ? 0.5 : 0.0));
      else              p = 0.0;
    }
  } else if (Mt == 2) { // maternal AA
    if      (y == 2) p = (g==1 ? 1.0 : (g==2 ? 0.5 : 0.0));
    else if (y == 1) p = (g==3 ? 1.0 : (g==2 ? 0.5 : 0.0));
    else              p = 0.0;
  } else if (Mt == 0) { // maternal aa
    if      (y == 0) p = (g==3 ? 1.0 : (g==2 ? 0.5 : 0.0));
    else if (y == 1) p = (g==1 ? 1.0 : (g==2 ? 0.5 : 0.0));
    else              p = 0.0;
  } else {
    return 1.0;
  }
  return (1.0 - epsilon) * p + (epsilon / 2.0) * (1.0 - p);
}

// Map class index 0..9 to left and right paternal single-locus genotypes gL,gR ∈ {2,1,0}
// Order: 0:AABB 1:AABb 2:AAbb 3:AaBB 4:AaBb_C 5:AaBb_R 6:Aabb 7:aaBB 8:aaBb 9:aabb
static inline void class_lr_genos(int s, int &gL, int &gR) {
  switch (s) {
  case 0: gL=2; gR=2; break; // AABB
  case 1: gL=2; gR=1; break; // AABb
  case 2: gL=2; gR=0; break; // AAbb
  case 3: gL=1; gR=2; break; // AaBB
  case 4: gL=1; gR=1; break; // AaBb_C
  case 5: gL=1; gR=1; break; // AaBb_R (same marginals as _C)
  case 6: gL=1; gR=0; break; // Aabb
  case 7: gL=0; gR=2; break; // aaBB
  case 8: gL=0; gR=1; break; // aaBb
  case 9: gL=0; gR=0; break; // aabb
  }
}

// Collapse interval weights Pi(:,k) to per-marker genotype frequencies on the LEFT marker k
static inline void collapse_left(const RcppParallel::RMatrix<double>& Pi, int k, double out_piL[3]) {
  out_piL[0]=out_piL[1]=out_piL[2]=0.0;
  for (int s=0; s<10; ++s) { int gL,gR; class_lr_genos(s,gL,gR); out_piL[gL] += Pi(s,k); }
  double z = out_piL[0]+out_piL[1]+out_piL[2]; if (z>0) { out_piL[0]/=z; out_piL[1]/=z; out_piL[2]/=z; }
}

// Collapse interval weights Pi(:,k) to per-marker genotype frequencies on the RIGHT marker k+1
static inline void collapse_right(const RcppParallel::RMatrix<double>& Pi, int k, double out_piR[3]) {
  out_piR[0]=out_piR[1]=out_piR[2]=0.0;
  for (int s=0; s<10; ++s) { int gL,gR; class_lr_genos(s,gL,gR); out_piR[gR] += Pi(s,k); }
  double z = out_piR[0]+out_piR[1]+out_piR[2]; if (z>0) { out_piR[0]/=z; out_piR[1]/=z; out_piR[2]/=z; }
}

// Two-locus kernel used for Model B responsibilities.
// Return Pr( (Yk, Yk1)=(y, y1) | maternal sends (h,h1), sire class s, dam genotypes Mk,Mk1 ).
// Use epsilon=0 here, since genotyping error is already handled in the HMM single-locus emission.
static inline double K2_prod(int y, int y1, int h, int h1, int s, int Mk, int Mk1) {
  int gL, gR; class_lr_genos(s, gL, gR);
  double p1 = py_hgM_obs(y,  h,  gL+1, Mk, 0.0);
  double p2 = py_hgM_obs(y1, h1, gR+1, Mk1, 0.0);
  return p1 * p2;
}

// emission vector e[h] = sum_g pi[g,t] * py_hgM_obs(y, h, g+1, M_t, epsilon)
// Requires: if y == NA_INTEGER, py_hgM_obs returns 1 for all (h,g) so that e[h] == 1.
static inline void emission_vec_wrapped(int y,
                                        const RMatrix<double>& pi, // 3 x T
                                        int t, int Mt, double epsilon,
                                        double out_e[2]) {
  const double pg0 = pi(0, t);
  const double pg1 = pi(1, t);
  const double pg2 = pi(2, t);
  for (int h = 0; h < 2; ++h) {
    double acc = 0.0;
    acc += pg0 * py_hgM_obs(y, h, 1, Mt, epsilon); // g=AA -> code 1
    acc += pg1 * py_hgM_obs(y, h, 2, Mt, epsilon); // g=Aa -> code 2
    acc += pg2 * py_hgM_obs(y, h, 3, Mt, epsilon); // g=aa -> code 3
    out_e[h] = (acc > 1e-15) ? acc : 1e-15;
  }
}

// ------------------------------------------------------------------
// Parallel E-step worker
// ------------------------------------------------------------------
struct EStepWorker : public Worker {
  // inputs
  RMatrix<int>     G;        // n x T child genotypes in {0,1,2,NA}
  RVector<int>     M;        // T dam genotypes in {0,1,2,NA}
  RVector<int>     phase;    // T-1, 0=repulsion, 1=coupling
  RVector<double>  r;        // T-1 recombination fractions
  RMatrix<double>  pi_emis;  // 3 x T, per-marker paternal genotype freqs used in emissions
  RMatrix<double>  Pi;       // 10 x (T-1), interval mixture for Model B
  std::string      paternal_mode;
  double           epsilon;
  int              T;

  // partial reductions
  std::vector<double> same, diff;    // length T-1, expected nonrecombinant vs recombinant under local phase
  std::vector<double> Ng0, Ng1, Ng2; // length T, expected paternal genotype counts per marker (per_marker)
  std::vector<double> N_A, N_a;      // length T, expected transmitted paternal A-/a-gamete counts (HWE/gametic)
  std::vector< std::array<double,10> > Ns; // length T-1, expected class counts per interval (Model B)
  double ll_sum;

  // constructor
  EStepWorker(IntegerMatrix G_,
              IntegerVector M_,
              IntegerVector phase_,
              NumericVector r_,
              NumericMatrix pi_emis_,
              NumericMatrix Pi_,
              std::string paternal_mode_,
              double epsilon_)
    : G(G_), M(M_), phase(phase_), r(r_),
      pi_emis(pi_emis_), Pi(Pi_), paternal_mode(paternal_mode_), epsilon(epsilon_),
      T(G_.ncol()),
      same(T-1, 0.0), diff(T-1, 0.0),
      Ng0(T, 0.0), Ng1(T, 0.0), Ng2(T, 0.0),
      N_A(T, 0.0), N_a(T, 0.0),
      Ns(T-1), ll_sum(0.0)
  {
    for (int k=0; k<T-1; ++k) for (int s=0; s<10; ++s) Ns[k][s] = 0.0;
  }

  // split constructor
  EStepWorker(EStepWorker& rhs, Split)
    : G(rhs.G), M(rhs.M), phase(rhs.phase), r(rhs.r),
      pi_emis(rhs.pi_emis), Pi(rhs.Pi), paternal_mode(rhs.paternal_mode), epsilon(rhs.epsilon),
      T(rhs.T),
      same(T-1, 0.0), diff(T-1, 0.0),
      Ng0(T, 0.0), Ng1(T, 0.0), Ng2(T, 0.0),
      N_A(T, 0.0), N_a(T, 0.0),
      Ns(T-1), ll_sum(0.0)
  {
    for (int k=0; k<T-1; ++k) for (int s=0; s<10; ++s) Ns[k][s] = 0.0;
  }

  inline double fb_child_and_accumulate(std::size_t i) {
    std::vector<double> alpha(2*T), beta(2*T), scale(T);
    std::vector<double> E(2*T);

    // emissions at all loci
    for (int t=0; t<T; ++t) {
      const int y = G(i,t);
      double e[2];
      emission_vec_wrapped(y, pi_emis, t, M[t], epsilon, e);
      E[2*t+0] = e[0];
      E[2*t+1] = e[1];
    }

    // forward init with g_h = 1/2
    alpha[0] = 0.5 * E[0];
    alpha[1] = 0.5 * E[1];
    double s0 = alpha[0] + alpha[1];
    if (s0 <= 0.0) s0 = 1e-15;
    scale[0] = s0;
    alpha[0] /= s0;
    alpha[1] /= s0;

    // forward
    for (int t=0; t<T-1; ++t) {
      const double p_same = (phase[t] == 1 ? (1.0 - r[t]) : r[t]);
      const double p_diff = 1.0 - p_same;
      const double a0 = alpha[2*t+0]*p_same + alpha[2*t+1]*p_diff;
      const double a1 = alpha[2*t+1]*p_same + alpha[2*t+0]*p_diff;
      const double e0 = E[2*(t+1)+0];
      const double e1 = E[2*(t+1)+1];
      alpha[2*(t+1)+0] = a0 * e0;
      alpha[2*(t+1)+1] = a1 * e1;
      double sc = alpha[2*(t+1)+0] + alpha[2*(t+1)+1];
      if (sc <= 0.0) sc = 1e-15;
      scale[t+1] = sc;
      alpha[2*(t+1)+0] /= sc;
      alpha[2*(t+1)+1] /= sc;
    }

    // backward init
    beta[2*(T-1)+0] = 1.0;
    beta[2*(T-1)+1] = 1.0;

    // backward
    for (int t=T-2; t>=0; --t) {
      const double p_same = (phase[t] == 1 ? (1.0 - r[t]) : r[t]);
      const double p_diff = 1.0 - p_same;
      const double e0 = E[2*(t+1)+0], e1 = E[2*(t+1)+1];
      const double b0 = p_same*e0*beta[2*(t+1)+0] + p_diff*e1*beta[2*(t+1)+1];
      const double b1 = p_same*e1*beta[2*(t+1)+1] + p_diff*e0*beta[2*(t+1)+0];
      const double sc = scale[t+1];
      beta[2*t+0] = b0 / sc;
      beta[2*t+1] = b1 / sc;
    }

    // interval expectations
    for (int t=0; t<T-1; ++t) {
      const double p_same = (phase[t] == 1 ? (1.0 - r[t]) : r[t]);
      const double p_diff = 1.0 - p_same;
      const double e0 = E[2*(t+1)+0], e1 = E[2*(t+1)+1];

      const double n00 = alpha[2*t+0]*p_same*e0*beta[2*(t+1)+0];
      const double n01 = alpha[2*t+0]*p_diff*e1*beta[2*(t+1)+1];
      const double n11 = alpha[2*t+1]*p_same*e1*beta[2*(t+1)+1];
      const double n10 = alpha[2*t+1]*p_diff*e0*beta[2*(t+1)+0];
      double s = n00 + n01 + n10 + n11;
      if (s <= 0.0) s = 1e-15;

      // same vs diff under local phase
      same[t] += (n00 + n11) / s;
      diff[t] += (n01 + n10) / s;

      // Model B responsibilities and expected class counts
      if (paternal_mode == "two_locus") {
        const double xi00 = n00 / s, xi01 = n01 / s, xi10 = n10 / s, xi11 = n11 / s;

        const int yk  = G(i,t);
        const int yk1 = G(i,t+1);
        const int Mk  = M[t];
        const int Mk1 = M[t+1];

        struct Pair { int j, jp; double w; };
        Pair pairs[4] = { {0,0,xi00}, {0,1,xi01}, {1,0,xi10}, {1,1,xi11} };

        for (auto pr : pairs) {
          if (pr.w <= 0.0) continue;

          double num[10];
          double Z = 0.0;

          for (int sidx=0; sidx<10; ++sidx) {
            double acc = 0.0;
            if (yk == NA_INTEGER || yk1 == NA_INTEGER) {
              // integrate over missing calls
              for (int ya=0; ya<3; ++ya) {
                if (yk  != NA_INTEGER && ya != yk ) continue;
                for (int yb=0; yb<3; ++yb) {
                  if (yk1 != NA_INTEGER && yb != yk1) continue;
                  acc += K2_prod(ya, yb, pr.j, pr.jp, sidx, Mk, Mk1);
                }
              }
            } else {
              acc = K2_prod(yk, yk1, pr.j, pr.jp, sidx, Mk, Mk1);
            }
            const double val = Pi(sidx,t) * acc + 1e-15;
            num[sidx] = val;
            Z += val;
          }
          if (Z <= 0.0) Z = 1e-15;
          for (int sidx=0; sidx<10; ++sidx) {
            const double R_s = num[sidx] / Z;
            Ns[t][sidx] += pr.w * R_s;
          }
        }
      }
    }

    // Paternal responsibilities per marker (Model A modes).
    //
    // HWE / gametic mode (public, default): accumulate DIRECT transmitted-gamete
    // counts. The offspring's paternal gamete is A with probability q_t and a with
    // 1 - q_t; the per-offspring, per-marker responsibility that the transmitted
    // gamete was A is
    //   rho_A = sum_h gamma(h) * q_t P(y|h,A) / b(h),
    //   rho_a = sum_h gamma(h) * (1-q_t) P(y|h,a) / b(h),
    //   b(h)  = q_t P(y|h,A) + (1-q_t) P(y|h,a),
    // where P(y|h,A) = P(y | maternal state h, paternal gamete A) is the AA
    // single-locus emission (AA transmits A with probability 1) and P(y|h,a) the
    // aa emission. This is the correct E-step of the one-parameter gametic model
    // whose M-step is q = (N_A + alpha)/(N_A + N_a + alpha + beta); it does NOT
    // route q through latent diploid paternal-genotype counts, which would make
    // the effective prior strength ~2x the documented pseudocount (a genotype
    // observation carries two alleles). The per-state b(h) normalizes INSIDE the
    // sum over h; a state contributes only when gamma(h) > 0 AND b(h) > 0 (an
    // impossible state at epsilon = 0 is skipped rather than dividing by zero).
    //
    // per_marker mode (legacy, for reproducing historical fits): unchanged
    // expected diploid genotype responsibilities R_g = sum_h gamma(h) pi_g
    // P(y|h,g) / b(h), b(h) = sum_g pi_g P(y|h,g).
    if (paternal_mode == "HWE") {
      for (int t=0; t<T; ++t) {
        const double gg0 = alpha[2*t+0]*beta[2*t+0];
        const double gg1 = alpha[2*t+1]*beta[2*t+1];
        const double z = gg0 + gg1;
        if (z <= 0.0) continue;                 // no maternal-state info here
        const double ga0 = gg0 / z;
        const double ga1 = gg1 / z;

        const int y  = G(i,t);
        if (y == NA_INTEGER) continue;
        const int Mt = M[t];
        // pi_emis is HWE(q_t) in this mode, so q_t = P(AA) + 0.5 P(Aa).
        const double q = clamp01(pi_emis(0,t) + 0.5*pi_emis(1,t));

        double nA = 0.0, na = 0.0;              // contributions to N_A / N_a
        if (ga0 > 0.0) {
          const double eA = py_hgM_obs(y,0,1,Mt,epsilon); // P(y | maternal 0, gamete A)
          const double ea = py_hgM_obs(y,0,3,Mt,epsilon); // P(y | maternal 0, gamete a)
          const double b0 = q*eA + (1.0-q)*ea;
          if (b0 > 0.0) { nA += ga0 * q * eA / b0; na += ga0 * (1.0-q) * ea / b0; }
        }
        if (ga1 > 0.0) {
          const double eA = py_hgM_obs(y,1,1,Mt,epsilon);
          const double ea = py_hgM_obs(y,1,3,Mt,epsilon);
          const double b1 = q*eA + (1.0-q)*ea;
          if (b1 > 0.0) { nA += ga1 * q * eA / b1; na += ga1 * (1.0-q) * ea / b1; }
        }
        N_A[t] += nA;
        N_a[t] += na;
      }
    } else if (paternal_mode == "per_marker") {
      for (int t=0; t<T; ++t) {
        const double gg0 = alpha[2*t+0]*beta[2*t+0];
        const double gg1 = alpha[2*t+1]*beta[2*t+1];
        const double z = gg0 + gg1;
        if (z <= 0.0) continue;                 // no maternal-state info here
        const double ga0 = gg0 / z;
        const double ga1 = gg1 / z;

        const int y  = G(i,t);
        if (y == NA_INTEGER) continue;
        const int Mt = M[t];

        double c0 = 0.0, c1 = 0.0, c2 = 0.0;    // contributions to Ng{AA,Aa,aa}
        if (ga0 > 0.0) {
          const double pAA = py_hgM_obs(y,0,1,Mt,epsilon);
          const double pAa = py_hgM_obs(y,0,2,Mt,epsilon);
          const double paa = py_hgM_obs(y,0,3,Mt,epsilon);
          const double b0 = pi_emis(0,t)*pAA + pi_emis(1,t)*pAa + pi_emis(2,t)*paa;
          if (b0 > 0.0) {
            c0 += ga0 * pi_emis(0,t) * pAA / b0;
            c1 += ga0 * pi_emis(1,t) * pAa / b0;
            c2 += ga0 * pi_emis(2,t) * paa / b0;
          }
        }
        if (ga1 > 0.0) {
          const double pAA = py_hgM_obs(y,1,1,Mt,epsilon);
          const double pAa = py_hgM_obs(y,1,2,Mt,epsilon);
          const double paa = py_hgM_obs(y,1,3,Mt,epsilon);
          const double b1 = pi_emis(0,t)*pAA + pi_emis(1,t)*pAa + pi_emis(2,t)*paa;
          if (b1 > 0.0) {
            c0 += ga1 * pi_emis(0,t) * pAA / b1;
            c1 += ga1 * pi_emis(1,t) * pAa / b1;
            c2 += ga1 * pi_emis(2,t) * paa / b1;
          }
        }
        Ng0[t] += c0;
        Ng1[t] += c1;
        Ng2[t] += c2;
      }
    }

    // log-likelihood from forward scales
    double ll = 0.0;
    for (int t=0; t<T; ++t) ll += std::log(scale[t]);
    return ll;
  }

  // parallel loop
  void operator()(std::size_t begin, std::size_t end) {
    double loc = 0.0;
    for (std::size_t i=begin; i<end; ++i) loc += fb_child_and_accumulate(i);
    ll_sum += loc;
  }

  // join
  void join(const EStepWorker& rhs) {
    for (int t=0; t<T-1; ++t) {
      same[t] += rhs.same[t];
      diff[t] += rhs.diff[t];
    }
    for (int t=0; t<T; ++t) {
      Ng0[t]  += rhs.Ng0[t];
      Ng1[t]  += rhs.Ng1[t];
      Ng2[t]  += rhs.Ng2[t];
      N_A[t]  += rhs.N_A[t];
      N_a[t]  += rhs.N_a[t];
    }
    for (int k=0; k<T-1; ++k)
      for (int s=0; s<10; ++s)
        Ns[k][s] += rhs.Ns[k][s];
    ll_sum += rhs.ll_sum;
  }
};

// ------------------------------------------------------------------
// Main EM driver
// ------------------------------------------------------------------

// [[Rcpp::export]]
Rcpp::List hmm_hs_cpp_parallel(Rcpp::IntegerMatrix G,           // n x T
                               Rcpp::IntegerVector M,           // T
                               Rcpp::IntegerVector phase_vec,   // T-1, 0=repulsion, 1=coupling
                               double r_start = 0.05,
                               std::string pi_mode = "per_marker", // for Model A
                               Rcpp::Nullable<Rcpp::NumericMatrix> pi_prior_in = R_NilValue, // 3 x T
                               double lambda = 20.0,
                               double epsilon = 1e-3,
                               double tol = 1e-6,
                               int    maxit = 200,
                               std::string paternal_mode = "per_marker", // "per_marker" | "HWE" | "two_locus"
                               Rcpp::Nullable<Rcpp::NumericMatrix> Pi_prior_in = R_NilValue, // 10 x (T-1)
                               Rcpp::Nullable<Rcpp::NumericVector> r_init = R_NilValue) // optional T-1 warm start for r
{
  const int n = G.nrow();
  const int T = G.ncol();
  if (M.size() != T) stop("M length must equal ncol(G)");
  if (phase_vec.size() != T-1) stop("phase_vec must have length T-1");
  for (int t=0; t<T-1; ++t) {
    const int v = phase_vec[t];
    if (!(v==0 || v==1)) stop("phase_vec must be 0 (repulsion) or 1 (coupling)");
  }
  if (!(paternal_mode == "per_marker" || paternal_mode == "HWE" || paternal_mode == "two_locus"))
    stop("paternal_mode must be 'per_marker', 'HWE', or 'two_locus'");

  // initialize r in [1e-6, 0.5]. A per-interval `r_init` (optional) overrides the
  // scalar r_start; used to warm-start the EM at an existing solution (diagnostic).
  NumericVector r(T-1, r_start);
  for (int t=0; t<T-1; ++t) r[t] = std::min(0.5, std::max(1e-6, r_start));
  if (r_init.isNotNull()) {
    NumericVector ri(r_init.get());
    if (ri.size() != T-1) stop("r_init must have length T-1");
    for (int t=0; t<T-1; ++t) r[t] = std::min(0.5, std::max(1e-6, ri[t]));
  }

  // Model A priors per marker
  NumericMatrix pi_prior(3, T);
  if (pi_prior_in.isNotNull()) {
    NumericMatrix tmp(pi_prior_in.get());
    if (tmp.nrow()!=3 || tmp.ncol()!=T) stop("pi_prior must be 3 x T");
    for (int j=0; j<T; ++j) {
      double s = tmp(0,j)+tmp(1,j)+tmp(2,j);
      if (s<=0) { pi_prior(0,j)=pi_prior(1,j)=pi_prior(2,j)=1.0/3.0; }
      else { pi_prior(0,j)=tmp(0,j)/s; pi_prior(1,j)=tmp(1,j)/s; pi_prior(2,j)=tmp(2,j)/s; }
    }
  } else {
    for (int j=0; j<T; ++j) pi_prior(0,j)=pi_prior(1,j)=pi_prior(2,j)=1.0/3.0;
  }

  if (paternal_mode != "two_locus") {
    // Enforce consistency between the chosen Model A flavors
    if (paternal_mode == "HWE" && pi_mode != "HWE")
      pi_mode = "HWE";
    if (paternal_mode == "per_marker" && pi_mode != "per_marker")
      pi_mode = "per_marker";
  }

  // initialize pi (3xT)
  // Rows: 0=AA, 1=Aa, 2=aa
  // HWE branch convention:
  // p denotes allele-A frequency, so HW genotypes are
  //   P(AA)=p^2, P(Aa)=2p(1-p), P(aa)=(1-p)^2.
  // The prior p is computed from prior genotype frequencies as p = P(AA) + 0.5*P(Aa).
  NumericMatrix pi(3, T);
  if (pi_mode == "per_marker") {
    for (int j=0; j<T; ++j){
      pi(0,j)=pi_prior(0,j);
      pi(1,j)=pi_prior(1,j);
      pi(2,j)=pi_prior(2,j);
    }
  } else if (pi_mode == "HWE") {
    // Interpret p as allele-A frequency
    for (int j=0; j<T; ++j){
      double p = pi_prior(0,j) + 0.5*pi_prior(1,j); // p = P(A)
      p = clamp01(p);
      pi(0,j) = p*p;            // AA
      pi(1,j) = 2.0*p*(1.0-p);  // Aa
      pi(2,j) = (1.0-p)*(1.0-p);// aa
    }
  } else {
    stop("pi_mode must be 'per_marker' or 'HWE'");
  }

  // Model B priors and parameters at intervals
  NumericMatrix Pi(10, std::max(T-1, 1));
  NumericMatrix Pi_prior(10, std::max(T-1, 1));
  if (paternal_mode == "two_locus") {
    if (Pi_prior_in.isNotNull()) {
      NumericMatrix tmp(Pi_prior_in.get());
      if (tmp.nrow()!=10 || tmp.ncol()!=T-1) stop("Pi_prior must be 10 x (T-1)");
      for (int k=0; k<T-1; ++k) {
        double s = 0.0; for (int sidx=0; sidx<10; ++sidx) s += tmp(sidx,k);
        if (s <= 0.0) { for (int sidx=0; sidx<10; ++sidx) Pi_prior(sidx,k) = 0.1; }
        else          { for (int sidx=0; sidx<10; ++sidx) Pi_prior(sidx,k) = tmp(sidx,k)/s; }
      }
    } else {
      for (int k=0; k<T-1; ++k) for (int sidx=0; sidx<10; ++sidx) Pi_prior(sidx,k) = 0.1;
    }
    // initialize Pi from prior
    for (int k=0; k<T-1; ++k) for (int sidx=0; sidx<10; ++sidx) Pi(sidx,k) = Pi_prior(sidx,k);
  }

  // working emission genotype frequencies used in the E-step
  NumericMatrix pi_emis(3, T);

  // helper to rebuild pi_emis from interval Pi (Model B)
  auto rebuild_pi_from_Pi = [&](NumericMatrix& out) {
    RMatrix<double> PiRM(Pi);
    for (int t=0; t<T; ++t) {
      double piL[3]={0,0,0}, piR_[3]={0,0,0};
      const bool hasL = (t <= T-2);  // interval k=t provides LEFT for marker t
      const bool hasR = (t >= 1);    // interval k=t-1 provides RIGHT for marker t
      if (hasL) collapse_left (PiRM, t,   piL);
      if (hasR) collapse_right(PiRM, t-1, piR_);
      if (hasL && hasR) {
        out(0,t) = 0.5*(piL[0] + piR_[0]);
        out(1,t) = 0.5*(piL[1] + piR_[1]);
        out(2,t) = 0.5*(piL[2] + piR_[2]);
      } else if (hasL) {
        out(0,t) = piL[0]; out(1,t) = piL[1]; out(2,t) = piL[2];
      } else { // hasR only
        out(0,t) = piR_[0]; out(1,t) = piR_[1]; out(2,t) = piR_[2];
      }
    }
  };

  double old_ll = -INFINITY;
  int it = 0, iters_done = 0;
  bool converged = false;
  std::string conv_reason = "maxit_reached";
  std::vector<double> ll_trace, dr_trace;

  for (it=1; it<=maxit; ++it) {
    iters_done = it;
    // choose emission table for this E-step
    if (paternal_mode == "two_locus") {
      rebuild_pi_from_Pi(pi_emis);
    } else {
      for (int j=0; j<T; ++j) {
        pi_emis(0,j)=pi(0,j);
        pi_emis(1,j)=pi(1,j);
        pi_emis(2,j)=pi(2,j);
      }
    }

    // parallel E-step
    EStepWorker worker(G, M, phase_vec, r, pi_emis, Pi, paternal_mode, epsilon);
    parallelReduce(0, static_cast<std::size_t>(n), worker);

    // expected transition counts per interval
    NumericVector same(T-1), diff(T-1);
    for (int t=0; t<T-1; ++t) { same[t] = worker.same[t]; diff[t] = worker.diff[t]; }

    // M-step for r with phase logic, clamp to [1e-6, 0.5]
    NumericVector r_new(T-1);
    for (int t=0; t<T-1; ++t) {
      const double tot = same[t] + diff[t];
      if (tot <= 1e-15) {
        r_new[t] = r[t];
      } else {
        r_new[t] = (phase_vec[t] == 1) ? (diff[t] / tot) : (same[t] / tot);
        r_new[t] = std::min(0.5, std::max(1e-6, r_new[t]));
      }
    }

    // M-step for paternal parameters
    double dpi = 0.0; // only used for Model A
    if (paternal_mode == "two_locus") {
      // collect expected class counts
      NumericMatrix Ns(10, std::max(T-1, 1));
      for (int k=0; k<T-1; ++k) for (int sidx=0; sidx<10; ++sidx) Ns(sidx,k) = worker.Ns[k][sidx];

      // Dirichlet shrinkage toward Pi_prior
      for (int k=0; k<T-1; ++k) {
        double s = 0.0;
        for (int sidx=0; sidx<10; ++sidx) { Pi(sidx,k) = Ns(sidx,k) + lambda * Pi_prior(sidx,k); s += Pi(sidx,k); }
        if (s <= 0.0) { for (int sidx=0; sidx<10; ++sidx) Pi(sidx,k) = 0.1; s = 1.0; }
        for (int sidx=0; sidx<10; ++sidx) Pi(sidx,k) /= s;
      }
      // rebuild pi_emis for next iteration is done at top of loop
    } else if (pi_mode == "per_marker") {
      NumericMatrix pi_new(3, T);
      for (int t=0; t<T; ++t) {
        pi_new(0,t) = worker.Ng0[t] + lambda * pi_prior(0,t);
        pi_new(1,t) = worker.Ng1[t] + lambda * pi_prior(1,t);
        pi_new(2,t) = worker.Ng2[t] + lambda * pi_prior(2,t);
        double s = pi_new(0,t)+pi_new(1,t)+pi_new(2,t);
        if (s <= 0.0) { pi_new(0,t)=pi_new(1,t)=pi_new(2,t)=1.0/3.0; }
        else { pi_new(0,t)/=s; pi_new(1,t)/=s; pi_new(2,t)/=s; }
      }
      // track max change
      for (int t=0; t<T; ++t)
        for (int g=0; g<3; ++g)
          dpi = std::max(dpi, std::fabs(pi_new(g,t) - pi(g,t)));
      // commit
      for (int j=0; j<T; ++j) { pi(0,j)=pi_new(0,j); pi(1,j)=pi_new(1,j); pi(2,j)=pi_new(2,j); }
    } else { // HWE / gametic: direct penalized (MAP) update on q
      NumericMatrix pi_new(3, T);
      for (int t = 0; t < T; ++t) {
        // Direct transmitted-gamete counts and pseudocount prior:
        //   q = (N_A + alpha)/(N_A + N_a + alpha + beta),
        // with pseudocount target q0 = alpha/(alpha+beta) and total lambda,
        // so alpha = lambda*q0, beta = lambda*(1-q0). This maximizes
        //   logLik + sum_t[alpha log q_t + beta log(1-q_t)]
        // (posterior mode under Beta(alpha+1, beta+1)); lambda = 0 is the MLE.
        const double NA = worker.N_A[t];
        const double Na = worker.N_a[t];
        const double q0 = pi_prior(0,t) + 0.5*pi_prior(1,t);   // pseudocount target q0 = P(A)
        const double a  = lambda * q0;
        const double b  = lambda * (1.0 - q0);
        double q_post = (NA + a) / std::max(NA + Na + a + b, 1e-12);
        q_post = clamp01(q_post);

        // Derived Hardy–Weinberg emission form with p = q_post = P(A):
        // rows are 0=AA, 1=Aa, 2=aa
        pi_new(0,t) = q_post * q_post;                 // AA
        pi_new(1,t) = 2.0 * q_post * (1.0 - q_post);  // Aa
        pi_new(2,t) = (1.0 - q_post) * (1.0 - q_post);// aa
      }

      // convergence metric and update
      for (int t = 0; t < T; ++t)
        for (int g = 0; g < 3; ++g)
          dpi = std::max(dpi, std::fabs(pi_new(g,t) - pi(g,t)));

      for (int j = 0; j < T; ++j) {
        pi(0,j) = pi_new(0,j);
        pi(1,j) = pi_new(1,j);
        pi(2,j) = pi_new(2,j);
      }
    }

    // total_ll = observed log-likelihood at the CURRENT (r, pi) used in this E-step.
    const double total_ll = worker.ll_sum;
    ll_trace.push_back(total_ll);
    (void) dpi;   // paternal decomposition change is NOT a convergence gate (it is
                  // a movement along a non-identifiable decomposition of q)

    // r-change (the identifiable parameter) for this M-step
    double dr = 0.0;
    for (int t=0; t<T-1; ++t) dr = std::max(dr, std::fabs(r_new[t] - r[t]));
    dr_trace.push_back(dr);

    // Convergence: RELATIVE change in observed log-likelihood AND stable r.
    const double rel = std::fabs(total_ll - old_ll) / (1.0 + std::fabs(old_ll));
    const bool conv = (it > 1) && (rel < tol) && (dr < tol);

    r = r_new; old_ll = total_ll;
    if (conv) { converged = true; conv_reason = "relative_loglik_and_r_stable"; break; }
  }

  // Recompute the observed log-likelihood at the FINAL returned parameters (r, pi),
  // i.e. AFTER the last M-step, so the reported logLik matches the returned params.
  if (paternal_mode == "two_locus") rebuild_pi_from_Pi(pi_emis);
  else for (int j=0; j<T; ++j) { pi_emis(0,j)=pi(0,j); pi_emis(1,j)=pi(1,j); pi_emis(2,j)=pi(2,j); }
  double final_ll;
  {
    EStepWorker fw(G, M, phase_vec, r, pi_emis, Pi, paternal_mode, epsilon);
    parallelReduce(0, static_cast<std::size_t>(n), fw);
    final_ll = fw.ll_sum;
  }

  // Penalized objective (gametic/HWE q-penalty only): logLik + sum_t[a log q + b log(1-q)].
  const bool has_pen = (paternal_mode == "HWE" && lambda > 0.0);
  double pen_obj = final_ll;
  if (has_pen) {
    for (int t=0; t<T; ++t) {
      const double q  = clamp01(pi(0,t) + 0.5*pi(1,t));
      const double q0 = pi_prior(0,t) + 0.5*pi_prior(1,t);
      const double a  = lambda * q0, b = lambda * (1.0 - q0);
      pen_obj += a * std::log(q) + b * std::log(1.0 - q);
    }
  }

  // annotate pi rows
  CharacterVector rn = CharacterVector::create("AA","Aa","aa");
  pi.attr("dimnames") = List::create(rn, R_NilValue);

  // outputs
  List out = List::create(
    _["r"]              = r,
    _["pi"]             = pi,
    _["pi_mode"]        = pi_mode,
    _["logLik"]         = final_ll,          // at the FINAL parameters
    _["penalized_obj"]  = has_pen ? Rcpp::wrap(pen_obj) : Rcpp::wrap(NA_REAL),
    _["converged"]      = converged,
    _["iters"]          = iters_done,        // never exceeds maxit
    _["conv_reason"]    = conv_reason,
    _["loglik_trace"]   = Rcpp::wrap(ll_trace),
    _["max_dr_trace"]   = Rcpp::wrap(dr_trace),
    _["epsilon"]        = epsilon,
    _["paternal_mode"]  = paternal_mode
  );
  if (paternal_mode == "two_locus") {
    out["Pi_interval"] = Pi;
    out["pi_emission"] = pi_emis;
  }
  return out;
}


// [[Rcpp::export]]
Rcpp::NumericVector gamma_cpp(Rcpp::IntegerMatrix G,           // n x T offspring genotypes 0/1/2/NA
                              Rcpp::IntegerVector M,           // length T maternal genotypes 0/1/2/NA
                              Rcpp::IntegerVector phase_vec,   // length T-1, 0 = repulsion, 1 = coupling
                              Rcpp::NumericVector r,           // length T-1 recombination fractions
                              Rcpp::NumericMatrix pi_emis,     // 3 x T paternal genotype freqs (AA,Aa,aa)
                              double epsilon = 1e-3)           // genotyping error in emissions
{
  using namespace Rcpp;
  using namespace RcppParallel;

  const int n = G.nrow();
  const int T = G.ncol();

  if (M.size() != T)                     stop("M length must equal ncol(G)");
  if (phase_vec.size() != T-1)           stop("phase_vec must have length T-1");
  if (r.size() != T-1)                   stop("r must have length T-1");
  if (pi_emis.nrow() != 3 || pi_emis.ncol() != T)
    stop("pi_emis must be 3 x T");

  RMatrix<int>     Gm(G);
  RVector<int>     Mm(M);
  RVector<int>     Ph(phase_vec);
  RVector<double>  Rv(r);
  RMatrix<double>  Piem(pi_emis);

  // Orientation by XOR over repulsion intervals:
  // hap0_is_A[t] = 1 means "hap0 corresponds to allele 'A' at marker t".
  IntegerVector hap0_is_A(T);
  hap0_is_A[0] = 1;
  for (int t = 1; t < T; ++t) {
    const int flip = (Ph[t-1] == 0) ? 1 : 0;  // repulsion => flip
    hap0_is_A[t] = hap0_is_A[t-1] ^ flip;
  }

  // Output array: [hap(2), marker(T), ind(n)], rows are hap0, hap1
  NumericVector out(Dimension(2, T, n));

  // Per-individual work buffers
  std::vector<double> alpha(2*T), beta(2*T), scale(T), E(2*T);

  for (int i = 0; i < n; ++i) {
    // Emissions for this individual
    for (int t = 0; t < T; ++t) {
      const int y = Gm(i, t);
      double e[2];
      emission_vec_wrapped(y, Piem, t, Mm[t], epsilon, e);
      E[2*t + 0] = e[0];   // h = 0  (maternal 'a')
      E[2*t + 1] = e[1];   // h = 1  (maternal 'A')
    }

    // Forward init (uniform prior over h)
    alpha[0] = 0.5 * E[0];
    alpha[1] = 0.5 * E[1];
    double s0 = alpha[0] + alpha[1];
    if (s0 <= 0.0) s0 = 1e-15;
    scale[0] = s0;
    alpha[0] /= s0; alpha[1] /= s0;

    // Forward
    for (int t = 0; t < T-1; ++t) {
      const double p_same = (Ph[t] == 1 ? (1.0 - Rv[t]) : Rv[t]);
      const double p_diff = 1.0 - p_same;

      const double a0 = alpha[2*t+0]*p_same + alpha[2*t+1]*p_diff;
      const double a1 = alpha[2*t+1]*p_same + alpha[2*t+0]*p_diff;

      const double e0 = E[2*(t+1)+0];
      const double e1 = E[2*(t+1)+1];

      alpha[2*(t+1)+0] = a0 * e0;
      alpha[2*(t+1)+1] = a1 * e1;

      double sc = alpha[2*(t+1)+0] + alpha[2*(t+1)+1];
      if (sc <= 0.0) sc = 1e-15;
      scale[t+1] = sc;
      alpha[2*(t+1)+0] /= sc;
      alpha[2*(t+1)+1] /= sc;
    }

    // Backward init
    beta[2*(T-1)+0] = 1.0;
    beta[2*(T-1)+1] = 1.0;

    // Backward
    for (int t = T-2; t >= 0; --t) {
      const double p_same = (Ph[t] == 1 ? (1.0 - Rv[t]) : Rv[t]);
      const double p_diff = 1.0 - p_same;

      const double e0 = E[2*(t+1)+0];
      const double e1 = E[2*(t+1)+1];

      const double b0 = p_same*e0*beta[2*(t+1)+0] + p_diff*e1*beta[2*(t+1)+1];
      const double b1 = p_same*e1*beta[2*(t+1)+1] + p_diff*e0*beta[2*(t+1)+0];

      const double sc = scale[t+1];
      beta[2*t+0] = b0 / sc;
      beta[2*t+1] = b1 / sc;
    }

    // Allele-state gammas and map to haplotypes via XOR orientation
    for (int t = 0; t < T; ++t) {
      const double g_a = alpha[2*t + 0] * beta[2*t + 0]; // state 0 (maternal 'a')
      const double g_A = alpha[2*t + 1] * beta[2*t + 1]; // state 1 (maternal 'A')
      const double z   = std::max(g_a + g_A, 1e-15);
      const double p_a = g_a / z;
      const double p_A = g_A / z;

      // hap0 = homolog 1, hap1 = homolog 2
      double hap0, hap1;
      if (hap0_is_A[t]) { hap0 = p_A; hap1 = p_a; }
      else              { hap0 = p_a; hap1 = p_A; }

      // write out[ hap, t, i ]
      out[0 + 2*(t + T*i)] = hap0;
      out[1 + 2*(t + T*i)] = hap1;
    }
  }

  // Attach orientation vector for reference (1 => hap0 is 'A' at that marker)
  out.attr("hap0_is_A") = hap0_is_A;
  return out;
}


// ------------------------------------------------------------------
// Joint multi-dam EM
//
// One SHARED recombination map r across dams; dam-specific phase and
// dam-specific paternal parameters. Reuses the single-dam EStepWorker per dam
// and pools expected recombination counts (converted from allele-state same/diff
// via each dam's phase) into one shared M-step. For D = 1 this reproduces the
// single-dam estimator. Convergence is on (log-likelihood, r) only (the paternal
// dpi criterion is tracked as a diagnostic but does not gate stopping).
// ------------------------------------------------------------------

// free version of the Model-B emission rebuild (3 x T from 10 x (T-1))
static inline void rebuild_pi_from_Pi_free(Rcpp::NumericMatrix& Pi,
                                           Rcpp::NumericMatrix& out, int T) {
  RcppParallel::RMatrix<double> PiRM(Pi);
  for (int t=0; t<T; ++t) {
    double piL[3]={0,0,0}, piR_[3]={0,0,0};
    const bool hasL = (t <= T-2);
    const bool hasR = (t >= 1);
    if (hasL) collapse_left (PiRM, t,   piL);
    if (hasR) collapse_right(PiRM, t-1, piR_);
    if (hasL && hasR) {
      out(0,t)=0.5*(piL[0]+piR_[0]); out(1,t)=0.5*(piL[1]+piR_[1]); out(2,t)=0.5*(piL[2]+piR_[2]);
    } else if (hasL) {
      out(0,t)=piL[0]; out(1,t)=piL[1]; out(2,t)=piL[2];
    } else {
      out(0,t)=piR_[0]; out(1,t)=piR_[1]; out(2,t)=piR_[2];
    }
  }
}

// [[Rcpp::export]]
Rcpp::List hmm_hs_joint_cpp(Rcpp::List G_list,                                  // list of n_d x T int matrices
                            Rcpp::List M_list,                                  // list of length-T int vectors
                            Rcpp::List phase_list,                              // list of length-(T-1) int (0/1)
                            double r_start = 0.05,
                            std::string pi_mode = "per_marker",
                            Rcpp::Nullable<Rcpp::List> pi_prior_list_in = R_NilValue, // list of 3 x T
                            double lambda = 20.0,
                            double epsilon = 1e-3,
                            double tol = 1e-6,
                            int    maxit = 200,
                            std::string paternal_mode = "per_marker",
                            Rcpp::Nullable<Rcpp::List> Pi_prior_list_in = R_NilValue) // list of 10 x (T-1)
{
  const int D = G_list.size();
  if (D < 1) stop("G_list must contain at least one dam");
  if (M_list.size()     != D) stop("M_list length must equal G_list length");
  if (phase_list.size() != D) stop("phase_list length must equal G_list length");
  if (!(paternal_mode == "per_marker" || paternal_mode == "HWE" || paternal_mode == "two_locus"))
    stop("paternal_mode must be 'per_marker', 'HWE', or 'two_locus'");

  // pull per-dam data; T fixed by dam 1, all dams must match (shared marker order)
  std::vector<IntegerMatrix> G(D);
  std::vector<IntegerVector> M(D);
  std::vector<IntegerVector> Ph(D);
  int T = -1;
  for (int d=0; d<D; ++d) {
    IntegerMatrix Gd = G_list[d];
    IntegerVector Md = M_list[d];
    IntegerVector Pd = phase_list[d];
    if (d==0) T = Gd.ncol();
    if (Gd.ncol() != T)     stop("all dams must share the same number of markers T (same order)");
    if ((int)Md.size() != T)   stop("M_list[[d]] length must equal T");
    if ((int)Pd.size() != T-1) stop("phase_list[[d]] length must equal T-1");
    for (int t=0; t<T-1; ++t) { int v=Pd[t]; if (!(v==0 || v==1)) stop("phase entries must be 0 or 1"); }
    G[d]=Gd; M[d]=Md; Ph[d]=Pd;
  }

  // shared r
  NumericVector r(T-1);
  for (int t=0; t<T-1; ++t) r[t] = std::min(0.5, std::max(1e-6, r_start));

  // keep Model A flavor consistent
  if (paternal_mode == "HWE")        pi_mode = "HWE";
  if (paternal_mode == "per_marker") pi_mode = "per_marker";

  // per-dam priors and parameters
  std::vector<NumericMatrix> pi_prior(D), pi(D);          // Model A: 3 x T
  std::vector<NumericMatrix> Pi_prior(D), Pi(D);          // Model B: 10 x (T-1)

  List pi_prior_list, Pi_prior_list;
  const bool has_pi_prior = pi_prior_list_in.isNotNull();
  if (has_pi_prior) pi_prior_list = pi_prior_list_in.get();
  const bool has_Pi_prior = Pi_prior_list_in.isNotNull();
  if (has_Pi_prior) Pi_prior_list = Pi_prior_list_in.get();

  for (int d=0; d<D; ++d) {
    // Model A prior
    NumericMatrix pp(3, T);
    if (has_pi_prior && d < pi_prior_list.size() && !Rf_isNull(pi_prior_list[d])) {
      NumericMatrix tmp = pi_prior_list[d];
      if (tmp.nrow()!=3 || tmp.ncol()!=T) stop("pi_prior_list[[d]] must be 3 x T");
      for (int j=0; j<T; ++j) {
        double s = tmp(0,j)+tmp(1,j)+tmp(2,j);
        if (s<=0) { pp(0,j)=pp(1,j)=pp(2,j)=1.0/3.0; }
        else      { pp(0,j)=tmp(0,j)/s; pp(1,j)=tmp(1,j)/s; pp(2,j)=tmp(2,j)/s; }
      }
    } else {
      for (int j=0; j<T; ++j) pp(0,j)=pp(1,j)=pp(2,j)=1.0/3.0;
    }
    pi_prior[d]=pp;

    // Model A init
    NumericMatrix p0(3, T);
    if (pi_mode == "per_marker") {
      for (int j=0; j<T; ++j) { p0(0,j)=pp(0,j); p0(1,j)=pp(1,j); p0(2,j)=pp(2,j); }
    } else { // HWE
      for (int j=0; j<T; ++j) { double p=clamp01(pp(0,j)+0.5*pp(1,j));
        p0(0,j)=p*p; p0(1,j)=2.0*p*(1.0-p); p0(2,j)=(1.0-p)*(1.0-p); }
    }
    pi[d]=p0;

    // Model B prior + init (kept valid even for Model A; unused there)
    NumericMatrix PP(10, std::max(T-1,1)), P0(10, std::max(T-1,1));
    if (paternal_mode == "two_locus") {
      if (has_Pi_prior && d < Pi_prior_list.size() && !Rf_isNull(Pi_prior_list[d])) {
        NumericMatrix tmp = Pi_prior_list[d];
        if (tmp.nrow()!=10 || tmp.ncol()!=T-1) stop("Pi_prior_list[[d]] must be 10 x (T-1)");
        for (int k=0; k<T-1; ++k) { double s=0; for (int si=0; si<10; ++si) s+=tmp(si,k);
          if (s<=0) { for (int si=0; si<10; ++si) PP(si,k)=0.1; }
          else      { for (int si=0; si<10; ++si) PP(si,k)=tmp(si,k)/s; } }
      } else {
        for (int k=0; k<T-1; ++k) for (int si=0; si<10; ++si) PP(si,k)=0.1;
      }
      for (int k=0; k<T-1; ++k) for (int si=0; si<10; ++si) P0(si,k)=PP(si,k);
    }
    Pi_prior[d]=PP; Pi[d]=P0;
  }

  double old_ll = -INFINITY;
  int it = 0, iters_done = 0;
  bool converged = false;
  std::string conv_reason = "maxit_reached";
  std::vector<double> ll_trace, dr_trace;

  for (it=1; it<=maxit; ++it) {
    iters_done = it;
    std::vector<double> Nrec(T-1, 0.0), Nnon(T-1, 0.0);
    double total_ll = 0.0;
    double dpi = 0.0;

    for (int d=0; d<D; ++d) {
      // emission table for this dam
      NumericMatrix pi_emis(3, T);
      if (paternal_mode == "two_locus") {
        rebuild_pi_from_Pi_free(Pi[d], pi_emis, T);
      } else {
        for (int j=0; j<T; ++j) { pi_emis(0,j)=pi[d](0,j); pi_emis(1,j)=pi[d](1,j); pi_emis(2,j)=pi[d](2,j); }
      }

      EStepWorker worker(G[d], M[d], Ph[d], r, pi_emis, Pi[d], paternal_mode, epsilon);
      parallelReduce(0, static_cast<std::size_t>(G[d].nrow()), worker);
      total_ll += worker.ll_sum;

      // pool recombination counts, converting allele-state same/diff via THIS dam's phase
      for (int t=0; t<T-1; ++t) {
        if (Ph[d][t] == 1) { Nrec[t] += worker.diff[t]; Nnon[t] += worker.same[t]; }
        else               { Nrec[t] += worker.same[t]; Nnon[t] += worker.diff[t]; }
      }

      // per-dam paternal M-step (dam-specific)
      if (paternal_mode == "two_locus") {
        for (int k=0; k<T-1; ++k) {
          double s=0;
          for (int si=0; si<10; ++si) { Pi[d](si,k)=worker.Ns[k][si]+lambda*Pi_prior[d](si,k); s+=Pi[d](si,k); }
          if (s<=0) { for (int si=0; si<10; ++si) Pi[d](si,k)=0.1; s=1.0; }
          for (int si=0; si<10; ++si) Pi[d](si,k)/=s;
        }
      } else if (pi_mode == "per_marker") {
        for (int t=0; t<T; ++t) {
          double a=worker.Ng0[t]+lambda*pi_prior[d](0,t);
          double b=worker.Ng1[t]+lambda*pi_prior[d](1,t);
          double c=worker.Ng2[t]+lambda*pi_prior[d](2,t);
          double s=a+b+c; if (s<=0) { a=b=c=1.0/3.0; s=1.0; }
          double na=a/s, nb=b/s, nc=c/s;
          dpi=std::max(dpi, std::fabs(na-pi[d](0,t)));
          dpi=std::max(dpi, std::fabs(nb-pi[d](1,t)));
          dpi=std::max(dpi, std::fabs(nc-pi[d](2,t)));
          pi[d](0,t)=na; pi[d](1,t)=nb; pi[d](2,t)=nc;
        }
      } else { // HWE / gametic: direct per-dam MAP update on q (gamete counts NOT pooled)
        for (int t=0; t<T; ++t) {
          double NA=worker.N_A[t], Na=worker.N_a[t];
          double q0=pi_prior[d](0,t)+0.5*pi_prior[d](1,t);      // pseudocount target q0 = P(A)
          double a=lambda*q0, b=lambda*(1.0-q0);
          double q_post=clamp01((NA+a)/std::max(NA+Na+a+b, 1e-12));
          double na=q_post*q_post, nb=2.0*q_post*(1.0-q_post), nc=(1.0-q_post)*(1.0-q_post);
          dpi=std::max(dpi, std::fabs(na-pi[d](0,t)));
          dpi=std::max(dpi, std::fabs(nb-pi[d](1,t)));
          dpi=std::max(dpi, std::fabs(nc-pi[d](2,t)));
          pi[d](0,t)=na; pi[d](1,t)=nb; pi[d](2,t)=nc;
        }
      }
    } // dams

    // shared M-step for r
    NumericVector r_new(T-1);
    double dr=0.0;
    for (int t=0; t<T-1; ++t) {
      double tot=Nrec[t]+Nnon[t];
      r_new[t] = (tot<=1e-15) ? r[t] : std::min(0.5, std::max(1e-6, Nrec[t]/tot));
      dr=std::max(dr, std::fabs(r_new[t]-r[t]));
    }

    ll_trace.push_back(total_ll);
    dr_trace.push_back(dr);
    (void) dpi;  // paternal decomposition change does NOT gate convergence

    // Convergence: RELATIVE observed-log-likelihood change AND stable shared r.
    const double rel = std::fabs(total_ll - old_ll) / (1.0 + std::fabs(old_ll));
    const bool conv = (it > 1) && (rel < tol) && (dr < tol);
    r = r_new; old_ll = total_ll;
    if (conv) { converged = true; conv_reason = "relative_loglik_and_r_stable"; break; }
  }

  const int iters = iters_done;   // never exceeds maxit

  // Recompute observed log-likelihood at the FINAL (shared r, per-dam pi).
  double final_ll = 0.0;
  for (int d=0; d<D; ++d) {
    NumericMatrix pe(3, T);
    if (paternal_mode == "two_locus") rebuild_pi_from_Pi_free(Pi[d], pe, T);
    else for (int j=0; j<T; ++j) { pe(0,j)=pi[d](0,j); pe(1,j)=pi[d](1,j); pe(2,j)=pi[d](2,j); }
    EStepWorker fw(G[d], M[d], Ph[d], r, pe, Pi[d], paternal_mode, epsilon);
    parallelReduce(0, static_cast<std::size_t>(G[d].nrow()), fw);
    final_ll += fw.ll_sum;
  }
  // Penalized objective (gametic/HWE q-penalty), summed over dams and markers.
  const bool has_pen = (paternal_mode == "HWE" && lambda > 0.0);
  double pen_obj = final_ll;
  if (has_pen) {
    for (int d=0; d<D; ++d) for (int t=0; t<T; ++t) {
      const double q  = clamp01(pi[d](0,t) + 0.5*pi[d](1,t));
      const double q0 = pi_prior[d](0,t) + 0.5*pi_prior[d](1,t);
      const double a  = lambda * q0, b = lambda * (1.0 - q0);
      pen_obj += a * std::log(q) + b * std::log(1.0 - q);
    }
  }

  // outputs
  CharacterVector rn = CharacterVector::create("AA","Aa","aa");
  List pi_list(D);
  for (int d=0; d<D; ++d) {
    NumericMatrix pid = pi[d];
    pid.attr("dimnames") = List::create(rn, R_NilValue);
    pi_list[d] = pid;
  }
  SEXP dn = G_list.attr("names");
  if (!Rf_isNull(dn)) pi_list.attr("names") = dn;

  List out = List::create(
    _["r"]              = r,
    _["pi_list"]        = pi_list,
    _["pi_mode"]        = pi_mode,
    _["logLik"]         = final_ll,          // at the FINAL parameters
    _["penalized_obj"]  = has_pen ? Rcpp::wrap(pen_obj) : Rcpp::wrap(NA_REAL),
    _["iters"]          = iters,             // never exceeds maxit
    _["converged"]      = converged,
    _["conv_reason"]    = conv_reason,
    _["loglik_trace"]   = Rcpp::wrap(ll_trace),
    _["max_dr_trace"]   = Rcpp::wrap(dr_trace),
    _["epsilon"]        = epsilon,
    _["paternal_mode"]  = paternal_mode,
    _["n_dams"]         = D
  );
  if (paternal_mode == "two_locus") {
    List Pi_list(D), emis_list(D);
    for (int d=0; d<D; ++d) {
      Pi_list[d] = Pi[d];
      NumericMatrix em(3, T); rebuild_pi_from_Pi_free(Pi[d], em, T); emis_list[d] = em;
    }
    if (!Rf_isNull(dn)) { Pi_list.attr("names") = dn; emis_list.attr("names") = dn; }
    out["Pi_interval_list"] = Pi_list;
    out["pi_emission_list"] = emis_list;
  }
  return out;
}


// ------------------------------------------------------------------
// Total HMM log-likelihood for one dam at a FIXED recombination vector r.
// Forward pass with per-position scaling, summed over offspring. Used by the
// per-dam map-heterogeneity (eta) test, which evaluates the likelihood at
// scaled maps without re-running EM.
// ------------------------------------------------------------------

// [[Rcpp::export]]
double loglik_hs_cpp(Rcpp::IntegerMatrix G,           // n x T offspring genotypes 0/1/2/NA
                     Rcpp::IntegerVector M,           // length T maternal genotypes
                     Rcpp::IntegerVector phase_vec,   // length T-1, 0 = repulsion, 1 = coupling
                     Rcpp::NumericVector r,           // length T-1 recombination fractions
                     Rcpp::NumericMatrix pi_emis,     // 3 x T paternal genotype freqs (AA,Aa,aa)
                     double epsilon = 1e-3)
{
  const int n = G.nrow();
  const int T = G.ncol();
  if (M.size() != T)            stop("M length must equal ncol(G)");
  if (phase_vec.size() != T-1)  stop("phase_vec must have length T-1");
  if (r.size() != T-1)          stop("r must have length T-1");
  if (pi_emis.nrow() != 3 || pi_emis.ncol() != T) stop("pi_emis must be 3 x T");

  RMatrix<double> Piem(pi_emis);

  double total_ll = 0.0;
  for (int i=0; i<n; ++i) {
    double e[2];
    emission_vec_wrapped(G(i,0), Piem, 0, M[0], epsilon, e);
    double a0 = 0.5*e[0], a1 = 0.5*e[1];
    double sc = a0 + a1; if (sc <= 0.0) sc = 1e-15;
    double ll = std::log(sc);
    a0 /= sc; a1 /= sc;

    for (int t=0; t<T-1; ++t) {
      const double p_same = (phase_vec[t] == 1 ? (1.0 - r[t]) : r[t]);
      const double p_diff = 1.0 - p_same;
      const double na0 = a0*p_same + a1*p_diff;
      const double na1 = a1*p_same + a0*p_diff;
      emission_vec_wrapped(G(i,t+1), Piem, t+1, M[t+1], epsilon, e);
      double b0 = na0 * e[0];
      double b1 = na1 * e[1];
      double s = b0 + b1; if (s <= 0.0) s = 1e-15;
      ll += std::log(s);
      a0 = b0/s; a1 = b1/s;
    }
    total_ll += ll;
  }
  return total_ll;
}
