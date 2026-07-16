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

// Coarse grid + bounded local refinement, for a possibly-MULTIMODAL objective on
// [lo, hi] (endpoints included). The profiled pairwise objective is a sum over dams
// of max(coupling, repulsion) log-likelihoods and need not be unimodal, so a single
// golden-section search is unreliable. Steps: (1) evaluate PAIRWISE_NGRID+1 grid
// points, (2) find every local maximum (endpoint-aware, including 0.5), (3) refine
// each promising neighborhood with a bounded golden section using the caller's tol
// and maxit, (4) keep the best evaluated objective. Diagnostics: best grid point,
// refined optimum, objective-evaluation count, at-boundary (0.5) and multimodality
// flags. The grid size is an internal constant, documented here.
static const int PAIRWISE_NGRID = 24;  // 25 points spanning [eps, 0.5], 0.5 included

struct PairOpt { double r_hat, f_hat, best_grid; int neval, n_local_max; bool at_half, multi_max; };

template <class F>
static inline PairOpt grid_refine_max(F f, double lo, double hi, int ngrid,
                                      double tol, int maxit) {
  PairOpt o; o.neval = 0;
  std::vector<double> xs(ngrid + 1), fs(ngrid + 1);
  for (int k = 0; k <= ngrid; ++k) {
    xs[k] = lo + (hi - lo) * ((double)k / (double)ngrid);
    fs[k] = f(xs[k]); ++o.neval;
  }
  int kbest = 0; for (int k = 1; k <= ngrid; ++k) if (fs[k] > fs[kbest]) kbest = k;
  o.best_grid = xs[kbest];

  std::vector<int> loc;                          // local maxima (endpoint-aware)
  for (int k = 0; k <= ngrid; ++k) {
    const bool okL = (k == 0)     || (fs[k] >= fs[k - 1]);
    const bool okR = (k == ngrid) || (fs[k] >= fs[k + 1]);
    if (okL && okR) loc.push_back(k);
  }
  o.n_local_max = (int)loc.size();

  double bx = xs[kbest], bf = fs[kbest];
  for (std::size_t t = 0; t < loc.size(); ++t) {
    const int m = loc[t];
    const double a = xs[std::max(0, m - 1)];
    const double b = xs[std::min(ngrid, m + 1)];
    int local_neval = 0;
    auto fc = [&](double r){ ++local_neval; return f(r); };
    const double gx = golden_max(fc, a, b, tol, maxit);
    const double gf = fc(gx);
    o.neval += local_neval;
    if (gf > bf) { bf = gf; bx = gx; }
  }
  o.r_hat = bx; o.f_hat = bf;
  // The no-linkage flag uses a FIXED boundary tolerance, independent of the optimizer
  // `tol`, so a coarse optimizer setting can never flag an r substantially below 0.5.
  static const double NO_LINKAGE_TOL = 1e-6;
  o.at_half = (bx >= hi - NO_LINKAGE_TOL);

  // conservative multimodality flag: >=2 local maxima whose GRID objective is within
  // a small margin of the best grid objective (a heads-up, not a guarantee).
  const double margin = 1e-6 * (1.0 + std::fabs(fs[kbest]));
  int ncomp = 0; for (std::size_t t = 0; t < loc.size(); ++t) if (fs[loc[t]] >= fs[kbest] - margin) ++ncomp;
  o.multi_max = (ncomp >= 2);
  return o;
}

// Single-locus marginal genotype probabilities P(Y=0/1/2) for a double-het dam
// (Aa x paternal gamete A-freq q), marginalized over the unobserved paternal
// gamete and the maternal transmission:
//   P(Y=2)=q/2, P(Y=1)=1/2, P(Y=0)=(1-q)/2.
// This is exactly sum_b PC[Y][b] = sum_b PR[Y][b]. CONDITIONAL ON FIXED q, a
// single-marker (partial) observation's likelihood contribution is therefore
// constant in r and in phase. Partial observations may still affect r and phase
// INDIRECTLY, because they contribute to the plug-in estimate of q (which is fixed
// before the r/phase optimization and enters both PC and PR).
static inline void marginal_single_dhet(double q, double out[3]) {
  out[0] = 0.5 * (1.0 - q);
  out[1] = 0.5;
  out[2] = 0.5 * q;
}

// -----------------------------------------------------------------------------
// Parallel worker: computes pairwise two-point fits for all (i,j), i<j
//
// Statistical core (this milestone):
//  * Dam-specific paternal gametic frequencies q_k^(d) are supplied precomputed
//    (per dam, per marker) and used to build each dam's joint genotype model;
//    they are NOT pooled across dams.
//  * Partial observations (one marker observed, the other missing) contribute to
//    the reported likelihood through the correct single-marker MARGINAL, and are
//    NOT discarded. Conditional on fixed q, that marginal is constant in r and
//    phase, so partial observations do not by themselves shift r-hat or the phase
//    LOD; they can still affect both INDIRECTLY through their contribution to the
//    plug-in estimate of q (computed before the r/phase optimization).
//  * Phase evidence is returned per dam (lod_ph_list, mom_phase_list); the pooled
//    lod_ph is the elementwise sum. An exact coupling-vs-repulsion tie yields
//    phase NA and phase-LOD 0 for that dam.
// The recombination-fraction optimizer (golden_max over sum_g max phase LL) and
// the two-marker likelihood kernel are unchanged.
// -----------------------------------------------------------------------------
struct PairwiseWorker : public RcppParallel::Worker {
  // inputs (read-only views)
  const std::vector<RMatrix<int>> kids;   // per-dam child genotypes [nKids x Tm]
  const std::vector<RVector<int>> moms;   // per-dam maternal genotypes [Tm]
  const std::vector<std::vector<double>>& qgk; // per-dam, per-marker q_k^(d) (NA where not het)
  const int Gp;                           // number of dams (populations)
  const int Tm;                           // number of markers
  const double lambda, tol, tiny;
  const int maxit;
  const bool want_diag;

  // outputs (write through)
  RMatrix<double> R;      // [Tm x Tm] pairwise r-hat
  RMatrix<double> LODR;   // [Tm x Tm] LOD vs r=0.5 (raw; may be tiny-negative)
  RMatrix<double> LODPH;  // [Tm x Tm] phase LOD at r-hat (sum over dams)
  RMatrix<double> LL;     // [Tm x Tm] log-likelihood at r-hat
  RMatrix<int>    NOLINK; // [Tm x Tm] 1 if r-hat is at the 0.5 boundary
  std::vector<RMatrix<int>>    momPhase;  // per-dam phase calls (C=1, R=0, NA)
  std::vector<RMatrix<double>> lodPhList; // per-dam phase LOD matrices
  // optional per-pair optimizer / count diagnostics (allocated only if want_diag)
  RMatrix<int>    NEVAL, MULTIM, NINFORM;
  RMatrix<double> BESTGR, NCOMP, NIONLY, NJONLY, NBOTH;

  PairwiseWorker(const std::vector<RMatrix<int>>& kids_,
                 const std::vector<RVector<int>>& moms_,
                 const std::vector<std::vector<double>>& qgk_,
                 int Gp_, int Tm_,
                 double lambda_, double tol_,
                 int maxit_, double tiny_, bool want_diag_,
                 RMatrix<double> R_, RMatrix<double> LODR_,
                 RMatrix<double> LODPH_, RMatrix<double> LL_, RMatrix<int> NOLINK_,
                 const std::vector<RMatrix<int>>& momPhase_,
                 const std::vector<RMatrix<double>>& lodPhList_,
                 RMatrix<int> NEVAL_, RMatrix<int> MULTIM_, RMatrix<int> NINFORM_,
                 RMatrix<double> BESTGR_, RMatrix<double> NCOMP_,
                 RMatrix<double> NIONLY_, RMatrix<double> NJONLY_, RMatrix<double> NBOTH_)
    : kids(kids_), moms(moms_), qgk(qgk_), Gp(Gp_), Tm(Tm_),
      lambda(lambda_), tol(tol_), tiny(tiny_), maxit(maxit_), want_diag(want_diag_),
      R(R_), LODR(LODR_), LODPH(LODPH_), LL(LL_), NOLINK(NOLINK_),
      momPhase(momPhase_), lodPhList(lodPhList_),
      NEVAL(NEVAL_), MULTIM(MULTIM_), NINFORM(NINFORM_),
      BESTGR(BESTGR_), NCOMP(NCOMP_), NIONLY(NIONLY_), NJONLY(NJONLY_), NBOTH(NBOTH_) {}

  void operator()(std::size_t begin, std::size_t end) {
    const double LOG10 = std::log(10.0);
    for (int i = static_cast<int>(begin); i < static_cast<int>(end); ++i) {
      for (int j = i + 1; j < Tm; ++j) {

        // Per-dam record: complete 3x3 counts, single-marker partial counts,
        // the constant (in r, phase) partial log-likelihood, and dam-specific q.
        struct DamRec {
          int  C[3][3];   // both markers observed
          int  niO[3];    // marker i observed, marker j missing
          int  njO[3];    // marker j observed, marker i missing
          long n_both;    // both markers missing (dhet dam)
          double ll_partial;  // sum niO[a] log marg_i[a] + njO[b] log marg_j[b]
          double q_i, q_j;
          bool is_dhet;
          bool has_complete;
        };
        std::vector<DamRec> dams(Gp);

        long n_complete_total = 0;   // complete observations across dhet dams
        long tot_complete = 0, tot_iOnly = 0, tot_jOnly = 0, tot_both = 0;
        int  n_informative = 0;      // dhet dams contributing >=1 complete observation

        for (int g = 0; g < Gp; ++g) {
          DamRec& d = dams[g];
          for (int a=0; a<3; ++a) for (int b=0; b<3; ++b) d.C[a][b] = 0;
          d.niO[0]=d.niO[1]=d.niO[2]=0;
          d.njO[0]=d.njO[1]=d.njO[2]=0;
          d.n_both = 0;
          d.ll_partial = 0.0;
          d.has_complete = false;

          const RVector<int>& Mi = moms[g];
          const int ma = Mi[i];
          const int mb = Mi[j];
          d.is_dhet = (ma == 1 && mb == 1);
          if (!d.is_dhet) continue;   // only double-het dams inform this pair

          d.q_i = clamp(qgk[g][i], 1e-6, 1.0 - 1e-6);
          d.q_j = clamp(qgk[g][j], 1e-6, 1.0 - 1e-6);

          const RMatrix<int>& Gg = kids[g];
          const int nKids = Gg.nrow();
          for (int r = 0; r < nKids; ++r) {
            const int yi = Gg(r, i);
            const int yj = Gg(r, j);
            const bool vi = (yi != NA_INTEGER && yi >= 0 && yi <= 2);
            const bool vj = (yj != NA_INTEGER && yj >= 0 && yj <= 2);
            if (vi && vj)      { d.C[yi][yj] += 1; d.has_complete = true; }
            else if (vi)       { d.niO[yi]   += 1; }   // i-only
            else if (vj)       { d.njO[yj]   += 1; }   // j-only
            else               { d.n_both    += 1; }   // both missing
          }

          // constant partial contribution via single-marker marginals
          double mi[3], mj[3];
          marginal_single_dhet(d.q_i, mi);
          marginal_single_dhet(d.q_j, mj);
          for (int a=0; a<3; ++a) if (d.niO[a] > 0) {
            double p = mi[a]; if (!(p > 0)) p = tiny;
            d.ll_partial += d.niO[a] * std::log(p);
          }
          for (int b=0; b<3; ++b) if (d.njO[b] > 0) {
            double p = mj[b]; if (!(p > 0)) p = tiny;
            d.ll_partial += d.njO[b] * std::log(p);
          }

          long nc = 0; for (int a=0;a<3;++a) for (int b=0;b<3;++b) nc += d.C[a][b];
          n_complete_total += nc;
          tot_complete += nc;
          tot_iOnly += d.niO[0] + d.niO[1] + d.niO[2];
          tot_jOnly += d.njO[0] + d.njO[1] + d.njO[2];
          tot_both  += d.n_both;
          if (d.has_complete) ++n_informative;
        }

        // A pair is fit only if some dhet dam has at least one complete (two-marker)
        // observation; otherwise r and phase are not identifiable -> leave NA.
        if (n_complete_total == 0) continue;

        // Objective: sum over dhet dams of max(ll_C, ll_R), using each dam's own
        // (already-estimated) q_i, q_j. Conditional on those fixed q, ll_partial is
        // constant in (r, phase) and does not affect r-hat; the partial data still
        // shaped r-hat indirectly, via its earlier contribution to q_i, q_j.
        // r is permitted up to EXACTLY 0.5 (true no-linkage null); it is not clamped
        // to 0.49. At r = 0.5 coupling and repulsion coincide, so the objective is
        // well defined and is the no-linkage likelihood.
        auto obj = [&](double r)->double {
          r = clamp(r, 1e-6, 0.5);
          double total = 0.0;
          for (int g=0; g<Gp; ++g) {
            if (!dams[g].is_dhet) continue;
            double PC[3][3], PR[3][3];
            joint_child_probs_3x3(0, r, dams[g].q_i, dams[g].q_j, PC);
            joint_child_probs_3x3(1, r, dams[g].q_i, dams[g].q_j, PR);
            const double llC = ll_3x3(dams[g].C, PC, tiny) + dams[g].ll_partial;
            const double llR = ll_3x3(dams[g].C, PR, tiny) + dams[g].ll_partial;
            total += (llC >= llR ? llC : llR);
          }
          return total;
        };

        // Grid + bounded local refinement over [eps, 0.5]; null evaluated at 0.5.
        const PairOpt po = grid_refine_max(obj, 1e-6, 0.5, PAIRWISE_NGRID, tol, maxit);
        const double r_hat  = po.r_hat;
        const double ll_hat = po.f_hat;
        const double ll_half = obj(0.5);                 // no-linkage null at EXACTLY 0.5
        const double lod_r  = (ll_hat - ll_half) / LOG10; // raw; may be tiny-negative from noise

        // Initialize outputs for this pair
        R(i,j)      = R(j,i)      = r_hat;
        LL(i,j)     = LL(j,i)     = ll_hat;
        LODR(i,j)   = LODR(j,i)   = lod_r;
        NOLINK(i,j) = NOLINK(j,i) = po.at_half ? 1 : 0;   // r-hat at the 0.5 boundary
        if (want_diag) {
          NEVAL(i,j)   = NEVAL(j,i)   = po.neval;
          BESTGR(i,j)  = BESTGR(j,i)  = po.best_grid;
          MULTIM(i,j)  = MULTIM(j,i)  = po.multi_max ? 1 : 0;
          NCOMP(i,j)   = NCOMP(j,i)   = (double)tot_complete;
          NIONLY(i,j)  = NIONLY(j,i)  = (double)tot_iOnly;
          NJONLY(i,j)  = NJONLY(j,i)  = (double)tot_jOnly;
          NBOTH(i,j)   = NBOTH(j,i)   = (double)tot_both;
          NINFORM(i,j) = NINFORM(j,i) = n_informative;
        }
        // At a fit pair, non-double-het dams (and ties) contribute phase-LOD 0, so
        // the pooled lod_ph equals the elementwise sum of lod_ph_list; their phase
        // call is NA. Unfit pairs and the diagonal stay NA in both.
        for (int g=0; g<Gp; ++g) {
          momPhase[g](i,j)  = momPhase[g](j,i)  = NA_INTEGER;
          lodPhList[g](i,j) = lodPhList[g](j,i) = 0.0;
        }

        // Per-dam phase LOD and phase call at r_hat, from THAT dam's likelihood
        // only. Partial terms cancel in (llC - llR), so use complete counts.
        double lod_ph_sum = 0.0;
        for (int g=0; g<Gp; ++g) {
          if (!dams[g].is_dhet) continue;

          double PC[3][3], PR[3][3];
          joint_child_probs_3x3(0, r_hat, dams[g].q_i, dams[g].q_j, PC);
          joint_child_probs_3x3(1, r_hat, dams[g].q_i, dams[g].q_j, PR);

          const double llC = ll_3x3(dams[g].C, PC, tiny);
          const double llR = ll_3x3(dams[g].C, PR, tiny);

          double lod; int phase;
          if (llC == llR) {            // exact tie (incl. no complete data)
            lod = 0.0; phase = NA_INTEGER;
          } else {
            lod   = std::fabs(llC - llR) / LOG10;
            phase = (llC > llR) ? 1 : 0;   // C=1, R=0
          }
          lodPhList[g](i,j) = lodPhList[g](j,i) = lod;
          momPhase[g](i,j)  = momPhase[g](j,i)  = phase;
          lod_ph_sum += lod;
        }
        LODPH(i,j) = LODPH(j,i) = lod_ph_sum;   // pooled = sum of per-dam LODs
      }

      // set diagonal for row i
      R(i,i) = 0.0;
      LL(i,i) = NA_REAL;
      LODR(i,i) = NA_REAL;
      LODPH(i,i) = NA_REAL;
      NOLINK(i,i) = NA_INTEGER;
      for (int g=0; g<Gp; ++g) { momPhase[g](i,i) = NA_INTEGER; lodPhList[g](i,i) = NA_REAL; }
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
                                                     double q0 = 0.5,
                                                     double tol = 1e-6,
                                                     int    maxit = 200,
                                                     double tiny  = 1e-12,
                                                     bool   share_q_across_dams = false,
                                                     bool   return_diagnostics = false,
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

  // Precompute dam-specific paternal gametic frequencies q_k^(d), per marker.
  // For a dam heterozygous (Aa) at marker k, a zero-error offspring genotype AA
  // observes one paternal A transmission and aa one paternal a transmission (Aa is
  // uninformative). Using EVERY offspring with an observed call at marker k
  // (including those missing at other markers):
  //   q_k^(d) = (n_AA + alpha) / (n_AA + n_aa + alpha + beta),
  //   alpha = lambda * q0,  beta = lambda * (1 - q0).
  // q is NA where the dam is not heterozygous at k (undefined there). These are
  // NOT pooled across dams. share_q_across_dams = true optionally pools the AA/aa
  // counts across dams (a single q per marker), retained for compatibility.
  const double a_pc = lambda * q0;
  const double b_pc = lambda * (1.0 - q0);
  std::vector<std::vector<double>> qgk(Gp, std::vector<double>(Tm, NA_REAL));
  {
    // per-marker pooled counts (only used when share_q_across_dams)
    std::vector<long> pooled_AA(Tm, 0), pooled_aa(Tm, 0);
    std::vector< std::vector<long> > nAA(Gp, std::vector<long>(Tm, 0));
    std::vector< std::vector<long> > naa(Gp, std::vector<long>(Tm, 0));
    for (int g=0; g<Gp; ++g) {
      IntegerVector Mg = momsR[g];
      IntegerMatrix Gg = kidsR[g];
      const int nKids = Gg.nrow();
      for (int k=0; k<Tm; ++k) {
        if (Mg[k] != 1) continue;              // dam not het at k -> q undefined
        long cAA=0, caa=0;
        for (int r=0; r<nKids; ++r) {
          const int y = Gg(r,k);
          if (y == 2) ++cAA; else if (y == 0) ++caa;   // Aa / NA: no info
        }
        nAA[g][k]=cAA; naa[g][k]=caa;
        pooled_AA[k]+=cAA; pooled_aa[k]+=caa;
      }
    }
    for (int g=0; g<Gp; ++g) {
      IntegerVector Mg = momsR[g];
      for (int k=0; k<Tm; ++k) {
        if (Mg[k] != 1) continue;
        long cAA = share_q_across_dams ? pooled_AA[k] : nAA[g][k];
        long caa = share_q_across_dams ? pooled_aa[k] : naa[g][k];
        const double denom = (double)cAA + (double)caa + a_pc + b_pc;
        double q = (denom > 0.0) ? ((double)cAA + a_pc) / denom : q0;
        qgk[g][k] = clamp(q, 1e-6, 1.0 - 1e-6);
      }
    }
  }

  // Outputs
  NumericMatrix R(Tm, Tm), LODR(Tm, Tm), LODPH(Tm, Tm), LL(Tm, Tm);
  IntegerMatrix NOLINK(Tm, Tm);
  std::fill(R.begin(),     R.end(),     NA_REAL);
  std::fill(LODR.begin(),  LODR.end(),  NA_REAL);
  std::fill(LODPH.begin(), LODPH.end(), NA_REAL);
  std::fill(LL.begin(),    LL.end(),    NA_REAL);
  std::fill(NOLINK.begin(),NOLINK.end(),NA_INTEGER);
  R.attr("dimnames")     = List::create(markers, markers);
  LODR.attr("dimnames")  = List::create(markers, markers);
  LODPH.attr("dimnames") = List::create(markers, markers);
  LL.attr("dimnames")    = List::create(markers, markers);
  NOLINK.attr("dimnames")= List::create(markers, markers);

  // Optional per-pair diagnostics: allocate full Tm x Tm only if requested, else a
  // 1x1 placeholder (kept out of the returned list unless return_diagnostics).
  const int dT = return_diagnostics ? Tm : 1;
  IntegerMatrix NEVAL(dT, dT), MULTIM(dT, dT), NINFORM(dT, dT);
  NumericMatrix BESTGR(dT, dT), NCOMP(dT, dT), NIONLY(dT, dT), NJONLY(dT, dT), NBOTH(dT, dT);
  if (return_diagnostics) {
    std::fill(NEVAL.begin(),  NEVAL.end(),  NA_INTEGER);
    std::fill(MULTIM.begin(), MULTIM.end(), NA_INTEGER);
    std::fill(NINFORM.begin(),NINFORM.end(),NA_INTEGER);
    std::fill(BESTGR.begin(), BESTGR.end(), NA_REAL);
    std::fill(NCOMP.begin(),  NCOMP.end(),  NA_REAL);
    std::fill(NIONLY.begin(), NIONLY.end(), NA_REAL);
    std::fill(NJONLY.begin(), NJONLY.end(), NA_REAL);
    std::fill(NBOTH.begin(),  NBOTH.end(),  NA_REAL);
    List dn = List::create(markers, markers);
    NEVAL.attr("dimnames")=dn; MULTIM.attr("dimnames")=dn; NINFORM.attr("dimnames")=dn;
    BESTGR.attr("dimnames")=dn; NCOMP.attr("dimnames")=dn; NIONLY.attr("dimnames")=dn;
    NJONLY.attr("dimnames")=dn; NBOTH.attr("dimnames")=dn;
  }

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

  // Per-dam phase-LOD outputs
  std::vector<NumericMatrix> lodPhR(Gp);
  for (int g=0; g<Gp; ++g) {
    NumericMatrix Lg(Tm, Tm);
    std::fill(Lg.begin(), Lg.end(), NA_REAL);
    Lg.attr("dimnames") = List::create(markers, markers);
    lodPhR[g] = Lg;
  }
  std::vector<RMatrix<double>> lodPhList; lodPhList.reserve(Gp);
  for (int g=0; g<Gp; ++g) lodPhList.emplace_back(lodPhR[g]);

  // Run parallel worker across i in [0, Tm-2]
  PairwiseWorker worker(kids, moms, qgk, Gp, Tm,
                        lambda, tol, maxit, tiny, return_diagnostics,
                        RMatrix<double>(R), RMatrix<double>(LODR),
                        RMatrix<double>(LODPH), RMatrix<double>(LL), RMatrix<int>(NOLINK),
                        momPhase, lodPhList,
                        RMatrix<int>(NEVAL), RMatrix<int>(MULTIM), RMatrix<int>(NINFORM),
                        RMatrix<double>(BESTGR), RMatrix<double>(NCOMP),
                        RMatrix<double>(NIONLY), RMatrix<double>(NJONLY), RMatrix<double>(NBOTH));

  parallelFor(0, Tm - 1, worker);

  // Set last diagonal cell
  R(Tm-1,Tm-1)     = 0.0;
  LL(Tm-1,Tm-1)    = NA_REAL;
  LODR(Tm-1,Tm-1)  = NA_REAL;
  LODPH(Tm-1,Tm-1) = NA_REAL;
  NOLINK(Tm-1,Tm-1)= NA_INTEGER;
  for (int g=0; g<Gp; ++g) { momPhaseR[g](Tm-1,Tm-1) = NA_INTEGER; lodPhR[g](Tm-1,Tm-1) = NA_REAL; }

  // Wrap per-dam phase matrices
  List mom_phase_list(Gp), lod_ph_list(Gp), q_list(Gp);
  for (int g=0; g<Gp; ++g) {
    mom_phase_list[g] = momPhaseR[g];
    lod_ph_list[g]    = lodPhR[g];
    NumericVector qv(Tm);
    for (int k=0; k<Tm; ++k) qv[k] = qgk[g][k];   // NA where dam not het at k
    qv.attr("names") = markers;
    q_list[g] = qv;
  }
  if (!Rf_isNull(dam_names)) {
    mom_phase_list.attr("names") = dam_names;
    lod_ph_list.attr("names")    = dam_names;
    q_list.attr("names")         = dam_names;
  }

  List out = List::create(
    _["r"]              = R,
    _["lod_r"]          = LODR,
    _["lod_ph"]         = LODPH,       // pooled = elementwise sum of lod_ph_list
    _["logLik"]         = LL,
    _["mom_phase_list"] = mom_phase_list,
    _["lod_ph_list"]    = lod_ph_list, // per-dam phase LOD matrices
    _["q_list"]         = q_list,      // per-dam, per-marker q_k^(d) (NA where not het)
    _["no_linkage"]     = NOLINK,      // 1 where r-hat is at the 0.5 boundary
    _["optimizer"]      = "grid+local-refine",
    _["n_grid"]         = PAIRWISE_NGRID
  );
  if (return_diagnostics) {
    out["diagnostics"] = List::create(
      _["n_eval"]        = NEVAL,      // objective evaluations per pair
      _["best_grid"]     = BESTGR,     // best grid point before refinement
      _["multi_maxima"]  = MULTIM,     // 1 if multiple comparable grid maxima
      _["n_complete"]    = NCOMP,      // two-marker observations (over dhet dams)
      _["n_i_only"]      = NIONLY,     // marker-i-only observations
      _["n_j_only"]      = NJONLY,     // marker-j-only observations
      _["n_both_missing"]= NBOTH,      // both-missing observations
      _["n_informative_dams"] = NINFORM
    );
  }
  return out;
}
