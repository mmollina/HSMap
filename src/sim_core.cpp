#include <Rcpp.h>
using namespace Rcpp;

static inline double clamp(double x, double lo, double hi){
  return (x<lo?lo:(x>hi?hi:x));
}

// [[Rcpp::export]]
Rcpp::List simulate_offspring_cpp(
    Rcpp::IntegerVector M,                 // 0/1/2/NA length T
    Rcpp::NumericMatrix pi_true,           // 3 x T (AA,Aa,aa)
    Rcpp::NumericVector r_true,            // length T-1
    int    n_offspring,
    double error_rate = 0.0,
    bool   keep_paths = false,
    Rcpp::Nullable<Rcpp::IntegerVector> z_phase_in = R_NilValue // <-- NEW (optional)
){
  RNGScope scope;
  const int T = M.size();
  if (T < 1) stop("M must have length >= 1");
  if (pi_true.nrow()!=3 || pi_true.ncol()!=T) stop("pi_true must be 3 x T");
  if (r_true.size() != std::max(0, T-1)) stop("r_true must have length T-1");
  if (error_rate < 0.0 || error_rate > 1.0) stop("error_rate must be in [0,1]");

  // orientation z_t in {0,1}: when mom is Aa, transmitted allele = H_t XOR z_t
  std::vector<int> z(T, 0);
  if (z_phase_in.isNotNull()){
    IntegerVector zin(z_phase_in.get());
    if ((int)zin.size()!=T) stop("z_phase must have length T");
    for (int t=0;t<T;++t) z[t] = (zin[t]==NA_INTEGER?0:(zin[t]&1));
  }


  IntegerMatrix G(n_offspring, T);
  std::fill(G.begin(), G.end(), NA_INTEGER);
  IntegerMatrix H;
  if (keep_paths){ H = IntegerMatrix(n_offspring, T); std::fill(H.begin(), H.end(), 0); }

  auto draw_paternal_allele = [&](int t)->int{
    double pAA = pi_true(0,t), pAa = pi_true(1,t), paa = pi_true(2,t);
    if (!R_finite(pAA) || !R_finite(pAa) || !R_finite(paa)) return NA_INTEGER;
    double s = pAA + pAa + paa;
    if (!R_finite(s) || s <= 0) return NA_INTEGER;
    pAA /= s; pAa /= s; // normalize
    double u = unif_rand();
    if (u < pAA) return 1;
    if (u < pAA + pAa) return (unif_rand() < 0.5 ? 1 : 0);
    return 0;
  };


  for (int n=0;n<n_offspring;++n){
    int h = (unif_rand()<0.5?0:1);
    if (keep_paths) H(n,0)=h;

    for (int t=0;t<T;++t){
      if (t>0){
        double r = clamp(r_true[t-1], 1e-12, 0.49);
        if (unif_rand() < r) h ^= 1;   // recombination flips H
        if (keep_paths) H(n,t)=h;
      }
      int Mt = M[t];
      if (Mt==NA_INTEGER) continue;

      // maternal allele (0/1)
      int a_mom = (Mt==0?0:(Mt==2?1:(h ^ z[t])));
      int a_pat = draw_paternal_allele(t);
      if (a_pat == NA_INTEGER) { G(n,t) = NA_INTEGER; continue; }
      int y = a_mom + a_pat; // 0/1/2


      // genotyping error -> replace with a *different* value
      if (error_rate>0.0 && unif_rand()<error_rate){
        if (y==0)      y = (unif_rand()<0.5?1:2);
        else if (y==1) y = (unif_rand()<0.5?0:2);
        else           y = (unif_rand()<0.5?0:1);
      }
      G(n,t)=y;
    }
  }
  return keep_paths ? List::create(_["G"]=G, _["H"]=H) : List::create(_["G"]=G);
}
