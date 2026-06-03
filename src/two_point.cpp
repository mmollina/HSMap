// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppParallel)]]

#include <Rcpp.h>
#include <RcppParallel.h>
#include <map>
#include <cmath>

using namespace Rcpp;
using namespace RcppParallel;

// -----------------------------------------------------------------------------
// Small utilities (used by the parallel worker)
// -----------------------------------------------------------------------------

// clamp x to [lo, hi]
static inline double clamp(double x, double lo, double hi) {
  return (x < lo) ? lo : ((x > hi) ? hi : x);
}

// Child single-locus genotype probs (0/1/2) given paternal allele A freq p
// and maternal transmitted allele m in {0,1} (0=a, 1=A).
// Output: out[0]=P(Y=0), out[1]=P(Y=1), out[2]=P(Y=2).
static inline void child_single_probs(double p, int m, double out[3]) {
  if (m == 0) { out[0] = 1.0 - p; out[1] = p;       out[2] = 0.0; }
  else        { out[0] = 0.0;     out[1] = 1.0 - p; out[2] = p;   }
}

// Joint child 3x3 probs P(Y_i=a, Y_j=b) for a maternal Aa×Aa double het,
// phaseType = 0 (coupling) or 1 (repulsion), recomb r, paternal A-freqs p_i, p_j.
// The maternal transmissions are mixed: (1-r)/2 over the two non-recombinant
// paths and r/2 over the two recombinant paths.
static inline void joint_child_probs_3x3(int phaseType, double r,
                                         double p_i, double p_j,
                                         double P[3][3]) {
  for (int a=0; a<3; ++a) for (int b=0; b<3; ++b) P[a][b] = 0.0;

  // maternal allele pairs (m_i, m_j) over the two non-recombinant and two recombinant paths
  const int pat_no_C[2][2] = {{0,0},{1,1}};
  const int pat_re_C[2][2] = {{0,1},{1,0}};
  const int pat_no_R[2][2] = {{1,0},{0,1}};
  const int pat_re_R[2][2] = {{0,0},{1,1}};

  const int (*pat_no)[2] = (phaseType==0 ? pat_no_C : pat_no_R);
  const int (*pat_re)[2] = (phaseType==0 ? pat_re_C : pat_re_R);

  const double w_no = 0.5 * (1.0 - r);
  const double w_re = 0.5 * r;

  // accumulate over the two non-recombinant paths
  for (int k=0; k<2; ++k) {
    const int mi = pat_no[k][0], mj = pat_no[k][1];
    double Fi[3], Fj[3];
    child_single_probs(p_i, mi, Fi);
    child_single_probs(p_j, mj, Fj);
    for (int a=0; a<3; ++a) for (int b=0; b<3; ++b) P[a][b] += w_no * (Fi[a]*Fj[b]);
  }

  // accumulate over the two recombinant paths
  for (int k=0; k<2; ++k) {
    const int mi = pat_re[k][0], mj = pat_re[k][1];
    double Fi[3], Fj[3];
    child_single_probs(p_i, mi, Fi);
    child_single_probs(p_j, mj, Fj);
    for (int a=0; a<3; ++a) for (int b=0; b<3; ++b) P[a][b] += w_re * (Fi[a]*Fj[b]);
  }

  // defensive normalization
  double s=0.0; for (int a=0; a<3; ++a) for (int b=0; b<3; ++b) s += P[a][b];
  if (s>0.0 && std::fabs(s-1.0) > 1e-12) {
    for (int a=0; a<3; ++a) for (int b=0; b<3; ++b) P[a][b] /= s;
  }
}

// Log-likelihood of 3x3 counts C against probs P; tiny guards zeros
static inline double ll_3x3(const int C[3][3], const double P[3][3], double tiny) {
  double ll = 0.0;
  for (int a=0; a<3; ++a) for (int b=0; b<3; ++b) {
    const int n = C[a][b];
    if (n <= 0) continue;
    double p = P[a][b];
    if (!(p > 0)) p = tiny;
    ll += n * std::log(p);
  }
  return ll;
}

// Golden-section maximization on [lo, hi] for a unimodal objective f
template <class F>
static inline double golden_max(F f, double lo, double hi,
                                double tol=1e-5, int maxit=200) {
  const double gr = (std::sqrt(5.0) + 1.0) / 2.0;
  double c = hi - (hi - lo) / gr;
  double d = lo + (hi - lo) / gr;
  double fc = f(c), fd = f(d);
  int it = 0;
  while ((hi - lo) > tol && it < maxit) {
    if (fc > fd) { hi = d; d = c; fd = fc; c = hi - (hi - lo) / gr; fc = f(c); }
    else         { lo = c; c = d; fc = fd; d = lo + (hi - lo) / gr; fd = f(d); }
    ++it;
  }
  return (fc > fd ? c : d);
}

// -----------------------------------------------------------------------------
// Parallel worker: computes pairwise two-point fits for all (i,j), i<j
// -----------------------------------------------------------------------------
struct PairwiseWorker : public RcppParallel::Worker {
  // inputs (read-only views)
  const std::vector<RMatrix<int>> kids;   // per-dam child genotypes [nKids x Tm]
  const std::vector<RVector<int>> moms;   // per-dam maternal genotypes [Tm]
  const int Gp;                           // number of dams (populations)
  const int Tm;                           // number of markers
  const double lambda, r_start, tol, tiny;
  const int maxit;
  const bool share_pi_across_dams;
  const bool verbose;

  // outputs (write through)
  RMatrix<double> R;      // [Tm x Tm] pairwise r-hat
  RMatrix<double> LODR;   // [Tm x Tm] LOD vs r=0.5
  RMatrix<double> LODPH;  // [Tm x Tm] phase LOD at r-hat (sum over dams)
  RMatrix<double> LL;     // [Tm x Tm] log-likelihood at r-hat
  std::vector<RMatrix<int>> momPhase; // per-dam phase calls (C=1, R=0, NA)

  PairwiseWorker(const std::vector<RMatrix<int>>& kids_,
                 const std::vector<RVector<int>>& moms_,
                 int Gp_, int Tm_,
                 double lambda_, double r_start_, double tol_,
                 int maxit_, double tiny_,
                 bool share_pi_, bool verbose_,
                 RMatrix<double> R_,
                 RMatrix<double> LODR_,
                 RMatrix<double> LODPH_,
                 RMatrix<double> LL_,
                 const std::vector<RMatrix<int>>& momPhase_)
    : kids(kids_), moms(moms_), Gp(Gp_), Tm(Tm_),
      lambda(lambda_), r_start(r_start_), tol(tol_), tiny(tiny_),
      maxit(maxit_),
      share_pi_across_dams(share_pi_), verbose(verbose_),
      R(R_), LODR(LODR_), LODPH(LODPH_), LL(LL_), momPhase(momPhase_) {}

  void operator()(std::size_t begin, std::size_t end) {
    for (int i = static_cast<int>(begin); i < static_cast<int>(end); ++i) {
      for (int j = i + 1; j < Tm; ++j) {

        // Build 3x3 counts for each dam and flag Aa×Aa double het moms
        struct DamRec { int C[3][3]; bool is_dhet; };
        std::vector<DamRec> dams; dams.reserve(Gp);

        long N_tot = 0, Yi2_sum = 0, Yj2_sum = 0;
        int any_used = 0;

        for (int g = 0; g < Gp; ++g) {
          DamRec d;
          for (int a=0; a<3; ++a) for (int b=0; b<3; ++b) d.C[a][b] = 0;

          const RVector<int>& Mi = moms[g];
          const int a = Mi[i];
          const int b = Mi[j];

          d.is_dhet = (a == 1 && b == 1);

          if (a != NA_INTEGER && b != NA_INTEGER) {
            const RMatrix<int>& Gg = kids[g];
            const int nKids = Gg.nrow();
            for (int r = 0; r < nKids; ++r) {
              const int yi = Gg(r, i);
              const int yj = Gg(r, j);
              if (yi == NA_INTEGER || yj == NA_INTEGER) continue;
              if (yi < 0 || yi > 2 || yj < 0 || yj > 2) continue;
              d.C[yi][yj] += 1;
              any_used = 1;
            }
          }

          if (d.is_dhet) {
            long n_all = 0;
            for (int a2=0; a2<3; ++a2) for (int b2=0; b2<3; ++b2) n_all += d.C[a2][b2];
            if (n_all > 0) {
              long yi2 = d.C[2][2] + d.C[2][1] + d.C[2][0];
              long yj2 = d.C[2][2] + d.C[1][2] + d.C[0][2];
              N_tot   += n_all;
              Yi2_sum += yi2;
              Yj2_sum += yj2;
            }
          }

          dams.push_back(d);
        }

        if (!any_used) {
          // leave NA at [i,j] and [j,i]
          continue;
        }

        // Choose paternal A-freqs from double-het marginals:
        // For mom Aa: P(Y=2) = p/2  => p = 2 * P(Y=2)
        if (N_tot == 0) {
          // nothing to estimate (no informative dams)
          continue;
        }
        const double p_i = clamp( 2.0 * ( (double)Yi2_sum / (double)N_tot ), 1e-6, 1.0 - 1e-6 );
        const double p_j = clamp( 2.0 * ( (double)Yj2_sum / (double)N_tot ), 1e-6, 1.0 - 1e-6 );

        // Objective: sum over dams of max(ll_C, ll_R)
        auto obj = [&](double r)->double {
          r = clamp(r, 1e-6, 0.49);
          double total = 0.0;
          for (int k=0; k<Gp; ++k) {
            if (!dams[k].is_dhet) continue;
            double PC[3][3], PR[3][3];
            joint_child_probs_3x3(0, r, p_i, p_j, PC);
            joint_child_probs_3x3(1, r, p_i, p_j, PR);
            const double llC = ll_3x3(dams[k].C, PC, tiny);
            const double llR = ll_3x3(dams[k].C, PR, tiny);
            total += (llC >= llR ? llC : llR);
          }
          return total;
        };

        const double r_hat  = golden_max(obj, 1e-6, 0.49, 1e-5, 200);
        const double ll_hat = obj(r_hat);
        const double ll_half = obj(0.5);
        const double lod_r  = (ll_hat - ll_half) / std::log(10.0);

        // Initialize outputs for this pair
        R(i,j)      = R(j,i)      = r_hat;
        LL(i,j)     = LL(j,i)     = ll_hat;
        LODR(i,j)   = LODR(j,i)   = lod_r;
        LODPH(i,j)  = LODPH(j,i)  = NA_REAL; // will fill below
        for (int g=0; g<Gp; ++g) { momPhase[g](i,j) = NA_INTEGER; momPhase[g](j,i) = NA_INTEGER; }

        // Phase LOD at r_hat and per-dam phase calls
        double lod_ph_sum = 0.0;
        for (int g=0; g<Gp; ++g) {
          if (!dams[g].is_dhet) continue;

          double PC[3][3], PR[3][3];
          joint_child_probs_3x3(0, r_hat, p_i, p_j, PC);
          joint_child_probs_3x3(1, r_hat, p_i, p_j, PR);

          const double llC = ll_3x3(dams[g].C, PC, tiny);
          const double llR = ll_3x3(dams[g].C, PR, tiny);
          const double lmax = std::max(llC, llR), lmin = std::min(llC, llR);
          lod_ph_sum += (lmax - lmin) / std::log(10.0);

          const int phase = (llC >= llR) ? 1 : 0; // C=1, R=0
          momPhase[g](i,j) = phase;
          momPhase[g](j,i) = phase;
        }
        LODPH(i,j) = LODPH(j,i) = lod_ph_sum;
      }

      // set diagonal for row i
      R(i,i) = 0.0;
      LL(i,i) = NA_REAL;
      LODR(i,i) = NA_REAL;
      LODPH(i,i) = NA_REAL;
      for (int g=0; g<Gp; ++g) momPhase[g](i,i) = NA_INTEGER;
    }
  }
};

// -----------------------------------------------------------------------------
// Exported parallel wrapper
//   - Aligns inputs by marker order
//   - Spawns PairwiseWorker across i in [0, Tm-2]
//   - Returns r, LOD_r, LOD_ph, logLik, and per-dam phase matrices
// -----------------------------------------------------------------------------

// [[Rcpp::export]]
Rcpp::List pairwise_rf_estimation_multi_parallel_cpp(Rcpp::List G_list,
                                                     Rcpp::List M_list,
                                                     double lambda = 20.0,
                                                     double r_start = 0.05,
                                                     double tol = 1e-6,
                                                     int    maxit = 200,
                                                     double tiny  = 1e-12,
                                                     bool   share_pi_across_dams = false,
                                                     bool   verbose = false) {
  using Rcpp::List; using Rcpp::IntegerVector; using Rcpp::IntegerMatrix;
  using Rcpp::NumericMatrix; using Rcpp::CharacterVector;

  const int Gp = G_list.size();
  if (Gp < 1) stop("G_list must have length >= 1");
  if (M_list.size() != Gp) stop("M_list length must equal G_list length");

  // Marker order = exactly names(M_list[[1]])
  IntegerVector M0 = M_list[0];
  SEXP nmsSEXP = M0.attr("names");
  if (Rf_isNull(nmsSEXP)) stop("M_list[[1]] must be a named integer vector (markers).");
  CharacterVector markers(nmsSEXP);
  const int Tm = markers.size();
  if (Tm < 2) stop("Need at least 2 markers.");

  // Dam names (for mom_phase_list naming)
  CharacterVector dam_names = G_list.attr("names");
  if (Rf_isNull(dam_names)) dam_names = M_list.attr("names");

  // Align moms
  std::vector<IntegerVector> momsR(Gp);
  for (int g=0; g<Gp; ++g) {
    IntegerVector Mg = M_list[g];
    if ((int)Mg.size() != Tm) stop("All M_list vectors must have same length.");
    momsR[g] = Mg;
  }

  // Align kids to 'markers'
  std::vector<IntegerMatrix> kidsR(Gp);
  for (int g=0; g<Gp; ++g) {
    IntegerMatrix Gg = G_list[g];
    SEXP dnSEXP = Gg.attr("dimnames");
    if (!Rf_isNull(dnSEXP)) {
      List dns(dnSEXP);
      CharacterVector cur = dns[1];
      if (cur.size() == Tm) {
        bool same = true;
        for (int c=0; c<Tm; ++c) if (cur[c] != markers[c]) { same=false; break; }
        if (!same) {
          std::map<std::string,int> pos;
          for (int c=0; c<Tm; ++c) pos[ as<std::string>(cur[c]) ] = c;
          IntegerMatrix Gnew(Gg.nrow(), Tm);
          for (int c=0; c<Tm; ++c) {
            const int from = pos.at( as<std::string>(markers[c]) );
            for (int r=0; r<Gg.nrow(); ++r) Gnew(r,c) = Gg(r,from);
          }
          Gnew.attr("dimnames") = List::create(dns[0], markers);
          kidsR[g] = Gnew;
        } else {
          kidsR[g] = Gg;
        }
      } else {
        kidsR[g] = Gg;
      }
    } else {
      kidsR[g] = Gg;
    }
  }

  // Wrap inputs for threads (RcppParallel views)
  std::vector<RVector<int>> moms; moms.reserve(Gp);
  std::vector<RMatrix<int>> kids; kids.reserve(Gp);
  for (int g=0; g<Gp; ++g) {
    moms.emplace_back(momsR[g]);
    kids.emplace_back(kidsR[g]);
  }

  // Outputs
  NumericMatrix R(Tm, Tm), LODR(Tm, Tm), LODPH(Tm, Tm), LL(Tm, Tm);
  std::fill(R.begin(),     R.end(),     NA_REAL);
  std::fill(LODR.begin(),  LODR.end(),  NA_REAL);
  std::fill(LODPH.begin(), LODPH.end(), NA_REAL);
  std::fill(LL.begin(),    LL.end(),    NA_REAL);
  R.attr("dimnames")     = List::create(markers, markers);
  LODR.attr("dimnames")  = List::create(markers, markers);
  LODPH.attr("dimnames") = List::create(markers, markers);
  LL.attr("dimnames")    = List::create(markers, markers);

  // Per-dam mom-phase outputs
  std::vector<IntegerMatrix> momPhaseR(Gp);
  for (int g=0; g<Gp; ++g) {
    IntegerMatrix M(Tm, Tm);
    std::fill(M.begin(), M.end(), NA_INTEGER);
    M.attr("dimnames") = List::create(markers, markers);
    momPhaseR[g] = M;
  }
  std::vector<RMatrix<int>> momPhase; momPhase.reserve(Gp);
  for (int g=0; g<Gp; ++g) momPhase.emplace_back(momPhaseR[g]);

  // Run parallel worker across i in [0, Tm-2]
  PairwiseWorker worker(kids, moms, Gp, Tm,
                        lambda, r_start, tol, maxit, tiny,
                        share_pi_across_dams, verbose,
                        RMatrix<double>(R), RMatrix<double>(LODR),
                        RMatrix<double>(LODPH), RMatrix<double>(LL),
                        momPhase);

  parallelFor(0, Tm - 1, worker);

  // Set last diagonal cell
  R(Tm-1,Tm-1)     = 0.0;
  LL(Tm-1,Tm-1)    = NA_REAL;
  LODR(Tm-1,Tm-1)  = NA_REAL;
  LODPH(Tm-1,Tm-1) = NA_REAL;
  for (int g=0; g<Gp; ++g) momPhaseR[g](Tm-1,Tm-1) = NA_INTEGER;

  // Wrap per-dam phase matrices
  List mom_phase_list(Gp);
  for (int g=0; g<Gp; ++g) mom_phase_list[g] = momPhaseR[g];
  if (!Rf_isNull(dam_names)) mom_phase_list.attr("names") = dam_names;

  return List::create(
    _["r"]              = R,
    _["lod_r"]          = LODR,
    _["lod_ph"]         = LODPH,
    _["logLik"]         = LL,
    _["mom_phase_list"] = mom_phase_list
  );
}
