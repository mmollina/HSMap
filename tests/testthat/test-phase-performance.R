# Preservation tests for the phase_from_pairwise() performance optimization
# (RSpectra leading eigenvector + C++ greedy). The optimized workflow must give the
# SAME phase as the original eigen()+R-greedy implementation, retained here as the
# reference. See dev/benchmarks/phase_bench.R for timing.

# ---- helpers ---------------------------------------------------------------

# Build a one-dam HSMap.tpt from a phase-sign matrix A (0/1/NA) and a LOD matrix.
.mk_tpt <- function(A, LOD, dam = "D") {
  Tn <- nrow(A); mk <- sprintf("m%04d", seq_len(Tn))
  dimnames(A) <- dimnames(LOD) <- list(mk, mk)
  fit <- list(mom_phase_list = stats::setNames(list(A), dam),
              lod_ph_list = stats::setNames(list(LOD), dam),
              lod_ph = LOD, markers = mk, n_dams = 1L)
  structure(list(fit = fit, markers = mk), class = "HSMap.tpt")
}

# Build a signed weighted graph from a component assignment. `comp[k]` = component id of
# marker k (id 0 => isolated / no supported edge). Within each component (>1 marker) all
# pairs are connected with a positive random weight and a coupling/repulsion sign implied
# by random true homolog labels (so the optimum recovers those labels up to sign).
.build_graph <- function(comp, seed = 1L, w_lo = 0.5, w_hi = 5) {
  set.seed(seed)
  Tn <- length(comp)
  A <- matrix(NA_real_, Tn, Tn)
  LOD <- matrix(0, Tn, Tn)
  h <- sample(0:1, Tn, replace = TRUE)                 # true homolog labels
  for (cc in unique(comp[comp > 0])) {
    idx <- which(comp == cc)
    if (length(idx) < 2L) next
    m <- length(idx)
    # vectorized dense block: coupling(1)/repulsion(0) from labels, positive weights
    A[idx, idx]   <- outer(h[idx], h[idx], function(a, b) as.integer(a == b))
    Wb <- matrix(stats::runif(m * m, w_lo, w_hi), m, m); Wb <- (Wb + t(Wb)) / 2
    LOD[idx, idx] <- Wb
  }
  mk <- sprintf("m%04d", seq_len(Tn))
  dimnames(A) <- dimnames(LOD) <- list(mk, mk)
  list(A = A, LOD = LOD, mk = mk, true_h = h)
}

# Reference phaser: exactly the phase_one() core, but using base eigen() and the R
# reference greedy .pf_greedy() (the pre-optimization path).
.ref_phase <- function(A, LOD, order, thr = 1e-8, max_passes = 50L, tol = 1e-9) {
  o <- order; A <- A[o, o, drop = FALSE]; LOD <- LOD[o, o, drop = FALSE]; Tn <- length(o)
  Sup <- (!is.na(A)) & is.finite(LOD) & (LOD > thr); diag(Sup) <- FALSE; Sup <- Sup | t(Sup)
  S <- matrix(0, Tn, Tn); S[!is.na(A) & A >= 0.5] <- 1; S[!is.na(A) & A < 0.5] <- -1
  W <- LOD; W[!is.finite(W)] <- 0; S[!Sup] <- 0; W[!Sup] <- 0; diag(S) <- 0; diag(W) <- 0
  J <- S * W; J <- (J + t(J)) / 2
  comp <- HSMap:::.pf_components(Sup)
  x <- rep(1, Tn); cobj <- 0
  for (cc in sort(unique(comp))) {
    idx <- which(comp == cc); if (length(idx) == 1L) next
    Jc <- J[idx, idx, drop = FALSE]; if (all(Jc == 0)) next
    ev <- eigen(Jc, symmetric = TRUE, only.values = FALSE)
    xc <- ifelse(ev$vectors[, 1] >= 0, 1, -1)
    gr <- HSMap:::.pf_greedy(Jc, xc, max_passes, tol); xc <- gr$x
    if (xc[1] != 1) xc <- -xc
    x[idx] <- xc; cobj <- cobj + gr$objective
  }
  clusters <- ifelse(x == 1, 0L, 1L)
  resolved <- if (Tn >= 2) comp[-Tn] == comp[-1] else logical(0)
  phase_vec <- if (Tn >= 2)
    ifelse(resolved, as.integer(clusters[-Tn] == clusters[-1]), NA_integer_) else integer(0)
  list(component = comp, clusters = clusters, phase_vec = phase_vec,
       objective = cobj, resolved = resolved)
}

# clusters equal up to an independent sign flip within each component?
.clusters_equiv <- function(cn, cr, comp) {
  for (cc in unique(comp)) {
    idx <- which(comp == cc)
    if (!(all(cn[idx] == cr[idx]) || all(cn[idx] == 1L - cr[idx]))) return(FALSE)
  }
  TRUE
}

# ---- component-level greedy equivalence (incl. dense 500) ------------------
test_that("C++ greedy matches the R reference greedy (same partition and objective)", {
  mkJ <- function(n, seed) { set.seed(seed)
    S <- matrix(sample(c(-1, 0, 1), n * n, TRUE, c(.3, .4, .3)), n, n)
    W <- matrix(stats::runif(n * n, 0, 5), n, n); J <- S * W; J <- (J + t(J)) / 2; diag(J) <- 0; J }
  for (cfg in list(c(500, 11), c(80, 12), c(30, 13), c(2, 14))) {
    J <- mkJ(cfg[1], cfg[2])
    v <- HSMap:::.pf_leading_eigvec(J); xc <- ifelse(v >= 0, 1, -1)
    grR <- HSMap:::.pf_greedy(J, xc, 50L, 1e-9)
    grC <- HSMap:::pf_greedy_cpp(J, as.integer(xc), 50L, 1e-9)
    expect_true(all(outer(grR$x, grR$x) == outer(grC$x, grC$x)))   # phase up to sign
    expect_gte(grC$objective, grR$objective - 1e-8)                # objective >= reference
  }
})

test_that("leading-eigenvector helper gives the same sign pattern as base eigen()", {
  set.seed(7); n <- 120
  S <- matrix(sample(c(-1, 0, 1), n * n, TRUE), n, n); W <- matrix(stats::runif(n * n, 0, 4), n, n)
  J <- S * W; J <- (J + t(J)) / 2; diag(J) <- 0
  vv <- HSMap:::.pf_leading_eigvec(J)
  ve <- eigen(J, symmetric = TRUE)$vectors[, 1]; if (ve[which.max(abs(ve))] < 0) ve <- -ve
  expect_true(all((vv >= 0) == (ve >= 0)))
  # small components (< 8) use the eigen fallback
  Js <- J[1:5, 1:5]; vs <- HSMap:::.pf_leading_eigvec(Js)
  es <- eigen(Js, symmetric = TRUE)$vectors[, 1]; if (es[which.max(abs(es))] < 0) es <- -es
  expect_true(all((vs >= 0) == (es >= 0)))
})

# ---- full workflow: mixed component structure -----------------------------
test_that("optimized phase_from_pairwise() matches the reference (mixed components)", {
  # components: 1:40 dense, 41:55 dense(15), 56:57 pair, 58/59/60 singletons (no edges)
  comp <- c(rep(1L, 40), rep(2L, 15), rep(3L, 2), 0L, 0L, 0L)
  g <- .build_graph(comp, seed = 21)
  tpt <- .mk_tpt(g$A, g$LOD)
  ord <- g$mk
  ph  <- phase_from_pairwise(tpt, order = ord)
  ref <- .ref_phase(g$A, g$LOD, ord)

  # identical connected components (relabel-invariant): same partition of markers
  expect_true(all(outer(ph$component, ph$component, "==") == outer(ref$component, ref$component, "==")))
  expect_identical(ph$n_components, length(unique(ref$component)))
  # identical resolved / unresolved intervals
  expect_identical(ph$resolved_interval, ref$resolved)
  # phase vectors exactly equal (phase_vec is sign-invariant), clusters up to component flip
  expect_identical(as.integer(ph$phase_vec), as.integer(ref$phase_vec))
  expect_true(.clusters_equiv(ph$clusters, ref$clusters, ph$component))
  # optimized objective >= reference
  expect_gte(ph$objective, ref$objective - 1e-8)

  # deterministic
  ph2 <- phase_from_pairwise(tpt, order = ord)
  expect_identical(as.integer(ph$phase_vec), as.integer(ph2$phase_vec))
  expect_identical(ph$clusters, ph2$clusters)
  expect_equal(ph$objective, ph2$objective)
})

# Note: the dense 500-marker component is exercised at the component level (the greedy +
# leading-eigenvector equivalence test above). A full phase_from_pairwise() call on a
# 500-marker COMPLETE graph is dominated by the (unchanged) connected-component traversal
# .pf_components(), not by the optimized spectral/greedy stages, so it is measured in the
# benchmark under dev/benchmarks/ rather than run in the ordinary suite.

# ---- sizes 1 and 2 --------------------------------------------------------
test_that("size-1 and size-2 components behave correctly", {
  comp <- c(1L, 1L, 0L, 0L)       # one pair + two singletons
  g <- .build_graph(comp, seed = 3)
  tpt <- .mk_tpt(g$A, g$LOD); ord <- g$mk
  ph <- phase_from_pairwise(tpt, order = ord)
  expect_identical(ph$component_sizes[order(-ph$component_sizes)][1], 2L)  # the pair
  expect_identical(sum(ph$component_sizes == 1L), 2L)                       # two singletons
  expect_false(is.na(ph$phase_vec[1]))                                     # the pair resolves (m1-m2)
  expect_true(is.na(ph$phase_vec[2]))                                      # m2-m3 crosses components
})

# ---- no supported edges ---------------------------------------------------
test_that("a graph with no supported edges is entirely unresolved (with a warning)", {
  Tn <- 15L; mk <- sprintf("m%04d", seq_len(Tn))
  A <- matrix(NA_real_, Tn, Tn); LOD <- matrix(0, Tn, Tn)   # no edges above threshold
  tpt <- .mk_tpt(A, LOD)
  expect_warning(ph <- phase_from_pairwise(tpt, order = mk), "No supported phase edges")
  expect_identical(ph$n_components, Tn)                     # all singletons
  expect_true(all(is.na(ph$phase_vec)))                    # all unresolved
  expect_false(ph$fully_resolved)
})
