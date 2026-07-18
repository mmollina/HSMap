# Reproducible synthetic benchmark for the phase_from_pairwise() optimization.
#
# Measures, on the SAME ~500-marker signed phase graph and with NO private data, the
# three stages of the per-component phase solve, comparing the pre-optimization path
# (base eigen() + the R greedy .pf_greedy()) with the optimized path (RSpectra leading
# eigenvector + the C++ pf_greedy_cpp()):
#
#   1. signed-graph construction  (shared / unchanged)
#   2. spectral initialization    (eigen()            vs RSpectra::eigs_sym())
#   3. greedy phase refinement    (.pf_greedy() [R]   vs pf_greedy_cpp() [C++])
#
# Run from the package root, e.g.:
#   Rscript -e "devtools::load_all(quiet=TRUE); source('dev/benchmarks/phase_bench.R')"

if (!exists("phase_from_pairwise")) {
  if (requireNamespace("devtools", quietly = TRUE)) devtools::load_all(quiet = TRUE) else library(HSMap)
}
stopifnot(requireNamespace("RSpectra", quietly = TRUE))

# median elapsed seconds over `reps` runs (a light-touch timer; no extra deps)
.med_time <- function(expr, reps = 7L) {
  e <- substitute(expr); pe <- parent.frame()
  ts <- vapply(seq_len(reps), function(i) {
    t0 <- Sys.time(); eval(e, pe); as.numeric(difftime(Sys.time(), t0, units = "secs"))
  }, numeric(1))
  stats::median(ts)
}

# Build a reproducible n-marker signed, LOD-weighted phase graph shaped like a real
# linkage group: banded supported edges (marker i linked to markers within `band`),
# LOD weight decaying with distance, true homolog labels following a coupling/repulsion
# random walk along the chromosome, and a fraction `noise` of edge signs flipped. The
# banded structure gives a weaker global spectral signal than a complete graph, so the
# greedy refinement does real work.
make_phase_graph <- function(n = 500L, band = 25L, noise = 0.12, rep_rate = 0.3,
                             seed = 1L, w_lo = 0.5, w_hi = 5) {
  set.seed(seed)
  h <- cumsum(c(0L, stats::rbinom(n - 1L, 1L, rep_rate))) %% 2L   # true labels
  A <- matrix(NA_real_, n, n); W <- matrix(0, n, n)
  for (i in seq_len(n - 1L)) {
    js <- (i + 1L):min(n, i + band)
    for (j in js) {
      s <- as.integer(h[i] == h[j])
      if (stats::runif(1) < noise) s <- 1L - s
      w <- max(stats::runif(1, w_lo, w_hi) * (1 - (j - i) / (band + 1)), 0.01)
      A[i, j] <- A[j, i] <- s; W[i, j] <- W[j, i] <- w
    }
  }
  list(A = A, LOD = W, true_h = h)
}

# Stage 1: signed-graph construction (the matrix assembly done in phase_one()).
build_J <- function(A, LOD, thr = 1e-8) {
  Tn <- nrow(A)
  Sup <- (!is.na(A)) & is.finite(LOD) & (LOD > thr); diag(Sup) <- FALSE; Sup <- Sup | t(Sup)
  S <- matrix(0, Tn, Tn); S[!is.na(A) & A >= 0.5] <- 1; S[!is.na(A) & A < 0.5] <- -1
  W <- LOD; W[!is.finite(W)] <- 0; S[!Sup] <- 0; W[!Sup] <- 0; diag(S) <- 0; diag(W) <- 0
  J <- S * W; (J + t(J)) / 2
}

run_bench <- function(n = 500L, noise = 0.15, seed = 1L, reps = 7L) {
  g <- make_phase_graph(n, noise, seed)

  t_build <- .med_time(build_J(g$A, g$LOD), reps)
  J <- build_J(g$A, g$LOD)

  # spectral
  t_eigen    <- .med_time(eigen(J, symmetric = TRUE, only.values = FALSE)$vectors[, 1], reps)
  t_rspectra <- .med_time(HSMap:::.pf_leading_eigvec(J), reps)

  # greedy (same start from the leading-eigenvector sign)
  xc <- ifelse(HSMap:::.pf_leading_eigvec(J) >= 0, 1, -1)
  grR <- HSMap:::.pf_greedy(J, xc, 50L, 1e-9)
  grC <- HSMap:::pf_greedy_cpp(J, as.integer(xc), 50L, 1e-9)
  t_greedyR <- .med_time(HSMap:::.pf_greedy(J, xc, 50L, 1e-9), reps)
  t_greedyC <- .med_time(HSMap:::pf_greedy_cpp(J, as.integer(xc), 50L, 1e-9), reps)

  old_total <- t_build + t_eigen + t_greedyR
  new_total <- t_build + t_rspectra + t_greedyC
  sg_old <- t_eigen + t_greedyR                 # optimized stages only
  sg_new <- t_rspectra + t_greedyC

  same_partition <- all(outer(grR$x, grR$x) == outer(grC$x, grC$x))

  cat(sprintf("phase_from_pairwise optimization benchmark  (n=%d markers, noise=%.2f, seed=%d)\n",
              n, noise, seed))
  cat("-----------------------------------------------------------------------\n")
  cat(sprintf("  graph construction (shared) : %8.4f s\n", t_build))
  cat(sprintf("  spectral   eigen()          : %8.4f s\n", t_eigen))
  cat(sprintf("  spectral   RSpectra         : %8.4f s   (%.1fx)\n", t_rspectra, t_eigen / t_rspectra))
  cat(sprintf("  greedy     R  (.pf_greedy)  : %8.4f s   [flips=%d, passes=%d]\n", t_greedyR, grR$n_flips, grR$iters))
  cat(sprintf("  greedy     C++ (pf_greedy)  : %8.4f s   (%.1fx) [flips=%d, passes=%d]\n",
              t_greedyC, t_greedyR / t_greedyC, grC$n_flips, grC$iters))
  cat("-----------------------------------------------------------------------\n")
  cat(sprintf("  spectral+greedy  OLD -> NEW : %8.4f s -> %8.4f s   (%.1fx)\n", sg_old, sg_new, sg_old / sg_new))
  cat(sprintf("  per-component    OLD -> NEW : %8.4f s -> %8.4f s   (%.1fx)\n", old_total, new_total, old_total / new_total))
  cat(sprintf("  objective  OLD=%.4f  NEW=%.4f  (obj_new >= obj_old: %s)\n",
              grR$objective, grC$objective, grC$objective >= grR$objective - 1e-8))
  cat(sprintf("  same phase partition (up to component sign): %s\n", same_partition))
  invisible(list(t_build = t_build, t_eigen = t_eigen, t_rspectra = t_rspectra,
                 t_greedyR = t_greedyR, t_greedyC = t_greedyC,
                 spectral_greedy_speedup = sg_old / sg_new,
                 component_speedup = old_total / new_total,
                 same_partition = same_partition))
}

# default run
invisible(run_bench(n = 500L, noise = 0.15, seed = 1L))
