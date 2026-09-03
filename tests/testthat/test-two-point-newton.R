## Validation of the EXPERIMENTAL phase-decomposed Newton optimizer for the
## single-dam two-point likelihood (branch optimization/two-point-concave-optimizer)
## against (a) the untouched production grid optimizer run on the SAME pair via
## two_point_pair_grid_cpp(), and (b) an INDEPENDENT high-accuracy R reference.
##
## The grid method is itself approximate, so disagreements are adjudicated by
## the reference, not by assuming the grid is right.
##
## Both harnesses consume identical sufficient statistics: the 3x3 complete-pair
## counts C3, one-marker partial counts niO/njO, and plug-in q_i/q_j.

## ---- simulate counts for ONE double-het dam directly from the model ---------
.tp_sim_counts <- function(n, r, phase = c("C", "R"), q_i, q_j,
                           miss = 0, seed = 1) {
  phase <- match.arg(phase)
  set.seed(seed)
  nonrec <- if (phase == "C") list(c(0L, 0L), c(1L, 1L)) else list(c(1L, 0L), c(0L, 1L))
  rec    <- if (phase == "C") list(c(0L, 1L), c(1L, 0L)) else list(c(0L, 0L), c(1L, 1L))
  C3  <- matrix(0L, 3, 3)
  niO <- integer(3); njO <- integer(3)
  for (k in seq_len(n)) {
    m <- if (stats::runif(1) < 1 - r) nonrec[[sample.int(2, 1)]] else rec[[sample.int(2, 1)]]
    yi <- m[1] + stats::rbinom(1, 1, q_i)      # maternal allele + paternal draw
    yj <- m[2] + stats::rbinom(1, 1, q_j)
    oi <- stats::runif(1) >= miss
    oj <- stats::runif(1) >= miss
    if (oi && oj)      C3[yi + 1L, yj + 1L] <- C3[yi + 1L, yj + 1L] + 1L
    else if (oi)       niO[yi + 1L] <- niO[yi + 1L] + 1L
    else if (oj)       njO[yj + 1L] <- njO[yj + 1L] + 1L
  }
  list(C3 = C3, niO = niO, njO = njO)
}

## ---- independent high-accuracy reference (pure R, fine grid + optimize) -----
.tp_ref <- function(C3, niO, njO, q_i, q_j, tiny = 1e-12, ngrid = 20001L) {
  q_i <- min(max(q_i, 1e-6), 1 - 1e-6)
  q_j <- min(max(q_j, 1e-6), 1 - 1e-6)
  F0 <- function(q) c(1 - q, q, 0)
  F1 <- function(q) c(0, 1 - q, q)
  A <- outer(F0(q_i), F0(q_j)) + outer(F1(q_i), F1(q_j))
  B <- outer(F0(q_i), F1(q_j)) + outer(F1(q_i), F0(q_j))
  mi <- c((1 - q_i) / 2, 0.5, q_i / 2)
  mj <- c((1 - q_j) / 2, 0.5, q_j / 2)
  ll_partial <- sum(niO * log(pmax(mi, tiny))) + sum(njO * log(pmax(mj, tiny)))

  ll_phase <- function(r, coupling) {
    P <- if (coupling) (1 - r) / 2 * A + r / 2 * B else (1 - r) / 2 * B + r / 2 * A
    sum(C3 * log(pmax(P, tiny))) + ll_partial
  }
  solve_phase <- function(coupling) {
    xs <- seq(1e-6, 0.5, length.out = ngrid)
    fs <- vapply(xs, ll_phase, numeric(1), coupling = coupling)
    k  <- which.max(fs)
    lo <- xs[max(1L, k - 1L)]; hi <- xs[min(ngrid, k + 1L)]
    op <- stats::optimize(function(r) ll_phase(r, coupling),
                          interval = c(lo, hi), maximum = TRUE, tol = 1e-12)
    cand_r  <- c(xs[k], op$maximum, 1e-6, 0.5)
    cand_ll <- vapply(cand_r, ll_phase, numeric(1), coupling = coupling)
    b <- which.max(cand_ll)
    list(r = cand_r[b], ll = cand_ll[b])
  }
  sC <- solve_phase(TRUE); sR <- solve_phase(FALSE)
  if (sC$ll >= sR$ll) list(r_hat = sC$r, ll_hat = sC$ll, r_C = sC$r, r_R = sR$r,
                           ll_C = sC$ll, ll_R = sR$ll)
  else                list(r_hat = sR$r, ll_hat = sR$ll, r_C = sC$r, r_R = sR$r,
                           ll_C = sC$ll, ll_R = sR$ll)
}

## ---- run both cpp harnesses on one scenario and compare ---------------------
.tp_cmp <- function(C3, niO = integer(3), njO = integer(3), q_i, q_j,
                    check_ref = TRUE, label = "case") {
  g <- HSMap:::two_point_pair_grid_cpp(C3, niO, njO, q_i, q_j)
  n <- HSMap:::two_point_pair_newton_cpp(C3, niO, njO, q_i, q_j)

  ## the new optimizer must never be meaningfully worse than the grid
  expect_gte(n$logLik, g$logLik - 1e-9, label = paste(label, "ll not worse"))
  ## same boundary classification and (near-tie aside) same phase call
  expect_identical(n$no_linkage, g$no_linkage, label = paste(label, "no_linkage"))
  phase_same <- identical(n$phase, g$phase)
  near_tie   <- max(n$lod_ph, g$lod_ph, na.rm = TRUE) < 1e-3
  expect_true(phase_same || near_tie, label = paste(label, "phase (or near-tie)"))
  ## LOD_r shares ll_half exactly, so it differs only through logLik
  expect_lt(abs(n$lod_r - g$lod_r), 1e-5, label = paste(label, "lod_r"))
  if (phase_same) {
    expect_lt(abs(n$lod_ph - g$lod_ph), 1e-3, label = paste(label, "lod_ph"))
  }
  ## adjudicate against the independent reference
  if (check_ref) {
    ref <- .tp_ref(C3, niO, njO, q_i, q_j)
    expect_lt(abs(n$logLik - ref$ll_hat), 1e-6, label = paste(label, "vs reference ll"))
    expect_gte(n$logLik, ref$ll_hat - 1e-8, label = paste(label, "not below reference"))
  }
  invisible(data.frame(label = label,
                       d_r = abs(n$r_hat - g$r_hat),
                       d_ll = abs(n$logLik - g$logLik),
                       d_lod_r = abs(n$lod_r - g$lod_r),
                       d_lod_ph = abs(n$lod_ph - g$lod_ph),
                       phase_same = phase_same,
                       r_grid = g$r_hat, r_newton = n$r_hat,
                       ne_grid = g$n_eval, ne_newton = n$n_eval))
}

## =============================================================================
test_that("newton vs grid: engineered scenarios agree (adjudicated by reference)", {
  diffs <- list()

  ## 1. interior optimum
  s <- .tp_sim_counts(300, 0.15, "C", 0.4, 0.4, seed = 11)
  diffs$interior <- .tp_cmp(s$C3, s$niO, s$njO, 0.4, 0.4, label = "interior")

  ## 2. optimum at the lower bound: pure non-recombinant coupling counts
  C3 <- matrix(0L, 3, 3); diag(C3) <- c(50L, 100L, 50L)
  d <- .tp_cmp(C3, q_i = 0.5, q_j = 0.5, label = "at-lo")
  g <- HSMap:::two_point_pair_grid_cpp(C3, integer(3), integer(3), 0.5, 0.5)
  n <- HSMap:::two_point_pair_newton_cpp(C3, integer(3), integer(3), 0.5, 0.5)
  expect_identical(n$r_hat, 1e-6)          # exact boundary return
  expect_identical(g$r_hat, 1e-6)          # grid keeps its first grid point
  diffs$at_lo <- d

  ## 3. optimum at EXACTLY 0.5 with an exact phase tie (balanced classes, q=0.5)
  C3 <- matrix(0L, 3, 3); C3[1, 1] <- 25L; C3[3, 3] <- 25L; C3[1, 3] <- 25L; C3[3, 1] <- 25L
  d <- .tp_cmp(C3, q_i = 0.5, q_j = 0.5, label = "at-half-tie")
  n <- HSMap:::two_point_pair_newton_cpp(C3, integer(3), integer(3), 0.5, 0.5)
  expect_identical(n$r_hat, 0.5)
  expect_identical(n$no_linkage, 1L)
  expect_true(is.na(n$phase))              # phases coincide at 0.5 -> tie -> NA
  expect_identical(n$lod_ph, 0)
  diffs$at_half <- d

  ## 4. no recombinant observations, asymmetric boundary winner:
  ##    counts only in (0,2)/(2,0) -> repulsion non-recombinant -> r_hat at lo,
  ##    phase = repulsion (0)
  C3 <- matrix(0L, 3, 3); C3[1, 3] <- 40L; C3[3, 1] <- 40L
  d <- .tp_cmp(C3, q_i = 0.5, q_j = 0.5, label = "no-recomb-repulsion")
  n <- HSMap:::two_point_pair_newton_cpp(C3, integer(3), integer(3), 0.5, 0.5)
  expect_identical(n$r_hat, 1e-6)
  expect_identical(n$phase, 0L)
  diffs$norec <- d

  ## 5. flat likelihood: only double-het offspring, q = 0.5 exactly.
  ## The argmax is NOT unique: every r has identical likelihood. The grid picks
  ## an arbitrary FP-ulp-noise-selected grid point (observed: 0.0208, because
  ## P[1][1] computed there lands an ulp above its value at 1e-6), while the
  ## Newton path deterministically returns the lower bound. Comparable r values
  ## are therefore NOT required here -- equal likelihood and consistent calls are.
  C3 <- matrix(0L, 3, 3); C3[2, 2] <- 200L
  g <- HSMap:::two_point_pair_grid_cpp(C3, integer(3), integer(3), 0.5, 0.5)
  n <- HSMap:::two_point_pair_newton_cpp(C3, integer(3), integer(3), 0.5, 0.5)
  expect_identical(n$r_hat, 1e-6)          # deterministic convention (lower bound)
  expect_lt(abs(n$logLik - g$logLik), 1e-9)  # same (constant) likelihood anywhere
  expect_lt(abs(n$lod_r), 1e-12)           # FP-ulp noise only, exactly like the
  expect_lt(abs(g$lod_r), 1e-12)           #   production "tiny-negative" comment
  expect_true(is.na(n$phase) && is.na(g$phase))
  diffs$flat <- .tp_cmp(C3, q_i = 0.5, q_j = 0.5, label = "flat")

  ## 6/7. coupling strongly preferred / repulsion strongly preferred
  s <- .tp_sim_counts(400, 0.05, "C", 0.4, 0.6, seed = 12)
  diffs$coupC <- .tp_cmp(s$C3, s$niO, s$njO, 0.4, 0.6, label = "coupling-strong")
  n <- HSMap:::two_point_pair_newton_cpp(s$C3, s$niO, s$njO, 0.4, 0.6)
  expect_identical(n$phase, 1L)
  s <- .tp_sim_counts(400, 0.05, "R", 0.4, 0.6, seed = 13)
  diffs$coupR <- .tp_cmp(s$C3, s$niO, s$njO, 0.4, 0.6, label = "repulsion-strong")
  n <- HSMap:::two_point_pair_newton_cpp(s$C3, s$niO, s$njO, 0.4, 0.6)
  expect_identical(n$phase, 0L)

  ## 8. near-0.5 from data simulated at r = 0.5
  s <- .tp_sim_counts(400, 0.5, "C", 0.4, 0.4, seed = 14)
  diffs$half_sim <- .tp_cmp(s$C3, s$niO, s$njO, 0.4, 0.4, label = "sim-at-half")

  ## 9/10. zero cells & sparse genotype classes (tiny families)
  s <- .tp_sim_counts(15, 0.1, "C", 0.3, 0.7, seed = 15)
  diffs$sparse15 <- .tp_cmp(s$C3, s$niO, s$njO, 0.3, 0.7, label = "sparse-n15")
  s <- .tp_sim_counts(8, 0.2, "R", 0.5, 0.5, seed = 16)
  diffs$sparse8 <- .tp_cmp(s$C3, s$niO, s$njO, 0.5, 0.5, label = "sparse-n8")

  ## 11. extreme paternal allele frequencies
  s <- .tp_sim_counts(300, 0.1, "C", 0.02, 0.98, seed = 17)
  diffs$extremeq <- .tp_cmp(s$C3, s$niO, s$njO, 0.02, 0.98, label = "extreme-q")

  ## 12. missing observations: partials present; they must not move the optimum
  s <- .tp_sim_counts(300, 0.12, "C", 0.4, 0.4, miss = 0.3, seed = 18)
  diffs$missing <- .tp_cmp(s$C3, s$niO, s$njO, 0.4, 0.4, label = "missing-30pc")
  n_part <- HSMap:::two_point_pair_newton_cpp(s$C3, s$niO, s$njO, 0.4, 0.4)
  n_none <- HSMap:::two_point_pair_newton_cpp(s$C3, integer(3), integer(3), 0.4, 0.4)
  expect_identical(n_part$r_hat, n_none$r_hat)                 # r unmoved
  expect_lt(abs((n_part$logLik - n_none$logLik) - n_part$ll_partial), 1e-9)
  expect_lt(abs(n_part$lod_r - n_none$lod_r), 1e-12)           # partials cancel

  ## 13. lambda = 0 territory: q at its clamp, probabilities near the tiny floor
  C3 <- matrix(0L, 3, 3); C3[1, 1] <- 60L; C3[3, 1] <- 2L; C3[1, 3] <- 2L; C3[2, 2] <- 20L
  diffs$floor <- .tp_cmp(C3, q_i = 1e-6, q_j = 1e-6, label = "q-at-clamp")

  ## 14. default-lambda-style q values on a realistic family size
  s <- .tp_sim_counts(331, 0.12, "C", 0.35, 0.65, seed = 19)
  diffs$default <- .tp_cmp(s$C3, s$niO, s$njO, 0.35, 0.65, label = "default-lambda")

  ## 15. one-sided near-floor: q_i at clamp only
  C3 <- matrix(0L, 3, 3); C3[3, 2] <- 5L; C3[2, 2] <- 50L; C3[1, 2] <- 5L; C3[2, 1] <- 10L
  diffs$floor1 <- .tp_cmp(C3, q_i = 1e-6, q_j = 0.5, label = "one-sided-floor")

  dd <- do.call(rbind, diffs)
  testthat::skip_if(nrow(dd) == 0)
  message(sprintf(
    "engineered cases: max|dr|=%.3g max|dll|=%.3g max|dlod_r|=%.3g max|dlod_ph|=%.3g",
    max(dd$d_r), max(dd$d_ll), max(dd$d_lod_r), max(dd$d_lod_ph)))
})

test_that("newton vs grid: no-complete-observation guard matches production", {
  C3 <- matrix(0L, 3, 3)
  g <- HSMap:::two_point_pair_grid_cpp(C3, c(5L, 3L, 2L), integer(3), 0.4, 0.4)
  n <- HSMap:::two_point_pair_newton_cpp(C3, c(5L, 3L, 2L), integer(3), 0.4, 0.4)
  for (f in c("r_hat", "logLik", "lod_r", "lod_ph", "phase", "no_linkage")) {
    expect_true(is.na(g[[f]]), label = paste("grid NA", f))
    expect_true(is.na(n[[f]]), label = paste("newton NA", f))
  }
})

test_that("phase symmetry ll_R(r) = ll_C(1 - r) holds in the reference model", {
  set.seed(31)
  for (k in 1:5) {
    q_i <- stats::runif(1, 0.05, 0.95); q_j <- stats::runif(1, 0.05, 0.95)
    F0 <- function(q) c(1 - q, q, 0); F1 <- function(q) c(0, 1 - q, q)
    A <- outer(F0(q_i), F0(q_j)) + outer(F1(q_i), F1(q_j))
    B <- outer(F0(q_i), F1(q_j)) + outer(F1(q_i), F0(q_j))
    r <- stats::runif(1, 0.01, 0.49)
    PC_r  <- (1 - r) / 2 * A + r / 2 * B
    PR_1r <- (1 - (1 - r)) / 2 * B + (1 - r) / 2 * A
    expect_lt(max(abs(PC_r - PR_1r)), 1e-15)
  }
})

test_that("newton vs grid: randomized sweep (seeded) stays within tolerance", {
  set.seed(2024)
  grid_par <- expand.grid(
    n    = c(8L, 25L, 60L, 150L, 400L),
    r    = c(0.005, 0.03, 0.1, 0.2, 0.35, 0.49, 0.5),
    phase = c("C", "R"),
    q_i  = c(0.05, 0.3, 0.5, 0.95),
    miss = c(0, 0.25),
    stringsAsFactors = FALSE
  )
  pick <- grid_par[sample.int(nrow(grid_par), 250L), ]
  worst <- data.frame(d_r = 0, d_ll = 0, d_lod_r = 0, d_lod_ph = 0)
  n_phase_diff <- 0L; n_checked <- 0L

  for (k in seq_len(nrow(pick))) {
    p <- pick[k, ]
    q_j <- stats::runif(1, 0.05, 0.95)
    s <- .tp_sim_counts(p$n, p$r, p$phase, p$q_i, q_j,
                        miss = p$miss, seed = 5000 + k)
    if (sum(s$C3) == 0) next
    g <- HSMap:::two_point_pair_grid_cpp(s$C3, s$niO, s$njO, p$q_i, q_j)
    n <- HSMap:::two_point_pair_newton_cpp(s$C3, s$niO, s$njO, p$q_i, q_j)
    n_checked <- n_checked + 1L

    expect_gte(n$logLik, g$logLik - 1e-9)
    expect_identical(n$no_linkage, g$no_linkage)
    phase_same <- identical(n$phase, g$phase)
    if (!phase_same) {
      n_phase_diff <- n_phase_diff + 1L
      expect_lt(max(n$lod_ph, g$lod_ph, na.rm = TRUE), 1e-3)  # only near-ties may differ
    }
    expect_lt(abs(n$lod_r - g$lod_r), 1e-5)
    ## adjudicate any r disagreement beyond optimizer tolerance via the reference
    if (abs(n$r_hat - g$r_hat) > 1e-4) {
      ref <- .tp_ref(s$C3, s$niO, s$njO, p$q_i, q_j)
      expect_lt(abs(n$logLik - ref$ll_hat), 1e-6)
      expect_gte(n$logLik, ref$ll_hat - 1e-8)
    }
    worst$d_r      <- max(worst$d_r,      abs(n$r_hat - g$r_hat))
    worst$d_ll     <- max(worst$d_ll,     abs(n$logLik - g$logLik))
    worst$d_lod_r  <- max(worst$d_lod_r,  abs(n$lod_r - g$lod_r))
    if (phase_same && !is.na(n$lod_ph))
      worst$d_lod_ph <- max(worst$d_lod_ph, abs(n$lod_ph - g$lod_ph))
  }
  message(sprintf(
    "sweep: %d scenarios | max|dr|=%.3g max|dll|=%.3g max|dlod_r|=%.3g max|dlod_ph|=%.3g | phase diffs (near-tie only): %d",
    n_checked, worst$d_r, worst$d_ll, worst$d_lod_r, worst$d_lod_ph, n_phase_diff))
  expect_gt(n_checked, 200L)
})
