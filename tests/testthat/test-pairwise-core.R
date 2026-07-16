# Tests for the pairwise statistical core (src/two_point.cpp):
#   * dam-specific paternal gametic frequency q_k^(d) from AA/aa transmissions,
#   * partial-missing likelihood via the correct single-marker marginal,
#   * dam-specific phase evidence (lod_ph_list, mom_phase_list) with pooled lod_ph,
#   * exact phase tie -> NA / LOD 0.
#
# An independent R implementation of the two-marker genotype model mirrors the C++
# kernel so the engine can be cross-checked against explicit marginalization.

.pw <- function(G, M, ...) HSMap:::pairwise_rf_estimation_multi_parallel_cpp(G, M, ...)

# P(Y = 0,1,2) for one offspring: dam-transmitted maternal allele m (0=a,1=A) and
# paternal gamete A-frequency p.
.child_single <- function(p, m) if (m == 0) c(1 - p, p, 0) else c(0, 1 - p, p)

# Two-marker joint genotype probs; rows Y_i = 0,1,2, cols Y_j = 0,1,2.
# phase: 0 = coupling, 1 = repulsion. Mirrors joint_child_probs_3x3() in C++.
.joint33 <- function(phase, r, pi_, pj) {
  P <- matrix(0, 3, 3)
  if (phase == 0) { no <- list(c(0, 0), c(1, 1)); re <- list(c(0, 1), c(1, 0)) }
  else            { no <- list(c(1, 0), c(0, 1)); re <- list(c(0, 0), c(1, 1)) }
  wno <- 0.5 * (1 - r); wre <- 0.5 * r
  for (m in no) P <- P + wno * outer(.child_single(pi_, m[1]), .child_single(pj, m[2]))
  for (m in re) P <- P + wre * outer(.child_single(pi_, m[1]), .child_single(pj, m[2]))
  P
}

# build an integer offspring matrix with the two marker names
.gm <- function(mat, mk = c("m1", "m2")) {
  m <- matrix(as.integer(mat), ncol = length(mk), byrow = FALSE)
  colnames(m) <- mk
  m
}
.dhet2 <- function() list(P1 = c(m1 = 1L, m2 = 1L))


# 1 ---------------------------------------------------------------------------
test_that("partial observations equal explicit probability marginalization", {
  RcppParallel::setThreadOptions(numThreads = 1)
  mk <- c("m1", "m2"); M <- .dhet2()
  G <- rbind(
    matrix(rep(c(2L, 2L), 10), ncol = 2, byrow = TRUE),  # complete AA/AA
    matrix(rep(c(0L, 0L),  8), ncol = 2, byrow = TRUE),  # complete aa/aa
    matrix(rep(c(2L, 1L),  3), ncol = 2, byrow = TRUE),  # complete AA/Aa
    cbind(c(2L, 0L, 1L), NA_integer_),                   # i-only
    cbind(NA_integer_, c(2L, 0L))                        # j-only
  )
  colnames(G) <- mk; storage.mode(G) <- "integer"
  res <- .pw(list(P1 = G), M, lambda = 0, q0 = 0.5)

  r_hat <- res$r[1, 2]; qi <- res$q_list$P1[[1]]; qj <- res$q_list$P1[[2]]
  ll_phase <- function(phase) {
    PJ <- .joint33(phase, r_hat, qi, qj); ll <- 0
    for (rr in seq_len(nrow(G))) {
      yi <- G[rr, 1]; yj <- G[rr, 2]
      if (!is.na(yi) && !is.na(yj)) ll <- ll + log(PJ[yi + 1, yj + 1])
      else if (!is.na(yi))          ll <- ll + log(sum(PJ[yi + 1, ]))  # marginalize j
      else if (!is.na(yj))          ll <- ll + log(sum(PJ[, yj + 1]))  # marginalize i
    }
    ll
  }
  expect_equal(res$logLik[1, 2], max(ll_phase(0), ll_phase(1)), tolerance = 1e-6)
})


# 2 ---------------------------------------------------------------------------
test_that("a genotype observed at only one marker contributes to q at that marker", {
  RcppParallel::setThreadOptions(numThreads = 1)
  M <- .dhet2()
  base <- rbind(matrix(rep(c(2L, 2L), 5), ncol = 2, byrow = TRUE),
                matrix(rep(c(0L, 0L), 5), ncol = 2, byrow = TRUE))
  colnames(base) <- c("m1", "m2"); storage.mode(base) <- "integer"
  res0 <- .pw(list(P1 = base), M, lambda = 0, q0 = 0.5)

  add <- matrix(c(2L, 2L, 2L, NA_integer_, NA_integer_, NA_integer_), ncol = 2)  # 3 AA at m1 only
  colnames(add) <- c("m1", "m2")
  res1 <- .pw(list(P1 = rbind(base, add)), M, lambda = 0, q0 = 0.5)

  expect_equal(res0$q_list$P1[[1]], 0.5,   tolerance = 1e-9)  # 5/(5+5)
  expect_equal(res1$q_list$P1[[1]], 8 / 13, tolerance = 1e-9) # (5+3)/((5+3)+5)
  expect_gt(res1$q_list$P1[[1]], res0$q_list$P1[[1]])
})


# 3 ---------------------------------------------------------------------------
test_that("two dams with different paternal allele frequencies get different q", {
  RcppParallel::setThreadOptions(numThreads = 1)
  M <- list(A = c(m1 = 1L, m2 = 1L), B = c(m1 = 1L, m2 = 1L))
  gA <- rbind(matrix(rep(c(2L, 2L), 18), ncol = 2, byrow = TRUE),
              matrix(rep(c(0L, 0L),  2), ncol = 2, byrow = TRUE))   # paternal A-rich
  gB <- rbind(matrix(rep(c(2L, 2L),  2), ncol = 2, byrow = TRUE),
              matrix(rep(c(0L, 0L), 18), ncol = 2, byrow = TRUE))   # paternal a-rich
  colnames(gA) <- colnames(gB) <- c("m1", "m2")
  storage.mode(gA) <- storage.mode(gB) <- "integer"
  res <- .pw(list(A = gA, B = gB), M, lambda = 0, q0 = 0.5)
  expect_equal(res$q_list$A[[1]], 0.9, tolerance = 1e-9)  # 18/20
  expect_equal(res$q_list$B[[1]], 0.1, tolerance = 1e-9)  # 2/20
  expect_true(res$q_list$A[[1]] != res$q_list$B[[1]])
})


# 4 ---------------------------------------------------------------------------
test_that("each dam receives its own phase LOD matrix", {
  RcppParallel::setThreadOptions(numThreads = 1)
  M <- list(A = c(m1 = 1L, m2 = 1L), B = c(m1 = 1L, m2 = 1L))
  gA <- rbind(matrix(rep(c(2L, 2L), 30), ncol = 2, byrow = TRUE),   # coupling
              matrix(rep(c(0L, 0L), 30), ncol = 2, byrow = TRUE))
  gB <- rbind(matrix(rep(c(2L, 0L), 30), ncol = 2, byrow = TRUE),   # repulsion
              matrix(rep(c(0L, 2L), 30), ncol = 2, byrow = TRUE))
  colnames(gA) <- colnames(gB) <- c("m1", "m2")
  storage.mode(gA) <- storage.mode(gB) <- "integer"
  res <- .pw(list(A = gA, B = gB), M, lambda = 0, q0 = 0.5)
  expect_type(res$lod_ph_list, "list")
  expect_length(res$lod_ph_list, 2)
  expect_equal(names(res$lod_ph_list), c("A", "B"))
  expect_equal(dim(res$lod_ph_list$A), c(2L, 2L))
  # A infers coupling (1), B infers repulsion (0); both have strong, distinct LODs
  expect_identical(res$mom_phase_list$A[1, 2], 1L)
  expect_identical(res$mom_phase_list$B[1, 2], 0L)
  expect_gt(res$lod_ph_list$A[1, 2], 5)
  expect_gt(res$lod_ph_list$B[1, 2], 5)
})


# 5 ---------------------------------------------------------------------------
test_that("a dam with weak evidence does not inherit the strong LOD of another dam", {
  RcppParallel::setThreadOptions(numThreads = 1)
  M <- list(strong = c(m1 = 1L, m2 = 1L), weak = c(m1 = 1L, m2 = 1L))
  gS <- rbind(matrix(rep(c(2L, 2L), 40), ncol = 2, byrow = TRUE),
              matrix(rep(c(0L, 0L), 40), ncol = 2, byrow = TRUE))   # strong coupling
  # weak dam: near-balanced coupling/repulsion evidence -> small phase LOD
  gW <- .gm(c(2L, 2L, 0L, 2L, 0L,   2L, 2L, 0L, 0L, 2L))            # rows (2,2)(2,2)(0,0)(2,0)(0,2)
  colnames(gS) <- c("m1", "m2"); storage.mode(gS) <- "integer"
  res <- .pw(list(strong = gS, weak = gW), M, lambda = 0, q0 = 0.5)
  expect_gt(res$lod_ph_list$strong[1, 2], 20)                       # strong dam
  expect_lt(res$lod_ph_list$weak[1, 2], 0.5 * res$lod_ph_list$strong[1, 2])  # clearly weaker
  # the weak dam does NOT carry the pooled LOD
  expect_false(isTRUE(all.equal(res$lod_ph_list$weak[1, 2], res$lod_ph[1, 2])))
})


# 6 ---------------------------------------------------------------------------
test_that("pooled lod_ph equals the elementwise sum of lod_ph_list", {
  RcppParallel::setThreadOptions(numThreads = 1)
  mk <- c("m1", "m2", "m3")
  # dam A het at all 3; dam B not het at m2 (contributes 0 there)
  M <- list(A = c(m1 = 1L, m2 = 1L, m3 = 1L), B = c(m1 = 1L, m2 = 2L, m3 = 1L))
  mkG <- function(rows) { m <- matrix(as.integer(unlist(rows)), ncol = 3, byrow = TRUE); colnames(m) <- mk; m }
  gA <- mkG(list(c(2,2,2), c(0,0,0), c(2,2,0), c(0,0,2), c(2,2,2), c(0,0,0)))
  gB <- mkG(list(c(2,2,0), c(0,2,2), c(2,2,2), c(0,2,0)))
  res <- .pw(list(A = gA, B = gB), M, lambda = 20, q0 = 0.5)
  s <- Reduce(`+`, res$lod_ph_list)
  expect_equal(res$lod_ph, s)               # matches at every entry, incl. NA pattern
})


# 7 ---------------------------------------------------------------------------
test_that("an exact phase tie returns NA phase and LOD 0", {
  RcppParallel::setThreadOptions(numThreads = 1)
  M <- .dhet2()
  # symmetric evidence: one each of (AA,AA),(aa,aa),(AA,aa),(aa,AA) with q=0.5
  g <- .gm(c(2L, 0L, 2L, 0L,   2L, 0L, 0L, 2L))  # cols m1=(2,0,2,0), m2=(2,0,0,2)
  res <- .pw(list(P1 = g), M, lambda = 0, q0 = 0.5)
  expect_true(is.na(res$mom_phase_list$P1[1, 2]))
  expect_equal(res$lod_ph_list$P1[1, 2], 0)
  expect_equal(res$lod_ph[1, 2], 0)
})


# 8 ---------------------------------------------------------------------------
test_that("existing downstream pairwise matrices and names remain available", {
  RcppParallel::setThreadOptions(numThreads = 1)
  mk <- c("m1", "m2"); M <- list(A = c(m1 = 1L, m2 = 1L))
  g <- rbind(matrix(rep(c(2L, 2L), 10), ncol = 2, byrow = TRUE),
             matrix(rep(c(0L, 0L), 10), ncol = 2, byrow = TRUE))
  colnames(g) <- mk; storage.mode(g) <- "integer"
  res <- .pw(list(A = g), M, lambda = 20, q0 = 0.5)
  expect_true(all(c("r", "lod_r", "lod_ph", "logLik", "mom_phase_list") %in% names(res)))
  expect_identical(dimnames(res$r),      list(mk, mk))
  expect_identical(dimnames(res$lod_ph), list(mk, mk))
  expect_true(is.matrix(res$mom_phase_list$A))
  expect_identical(names(res$mom_phase_list), "A")
  # phase_from_pairwise (unmodified) still consumes lod_ph + mom_phase_list
  tpt <- structure(list(fit = res, markers = mk), class = "HSMap.tpt")
  ph <- phase_from_pairwise(tpt, order = mk)
  expect_s3_class(ph, "HSMap.phased")
})


# 9 ---------------------------------------------------------------------------
test_that("complete-data likelihood/LOD use the unchanged kernel given the same q", {
  RcppParallel::setThreadOptions(numThreads = 1)
  M <- .dhet2()
  g <- rbind(matrix(rep(c(2L, 2L), 12), ncol = 2, byrow = TRUE),
             matrix(rep(c(0L, 0L), 12), ncol = 2, byrow = TRUE),
             matrix(rep(c(2L, 0L),  3), ncol = 2, byrow = TRUE))
  colnames(g) <- c("m1", "m2"); storage.mode(g) <- "integer"
  res <- .pw(list(P1 = g), M, lambda = 0, q0 = 0.5)
  r_hat <- res$r[1, 2]; qi <- res$q_list$P1[[1]]; qj <- res$q_list$P1[[2]]
  C <- matrix(0, 3, 3)
  for (rr in seq_len(nrow(g))) C[g[rr, 1] + 1, g[rr, 2] + 1] <- C[g[rr, 1] + 1, g[rr, 2] + 1] + 1
  llp <- function(P) sum(ifelse(C > 0, C * log(pmax(P, 1e-12)), 0))
  llC <- llp(.joint33(0, r_hat, qi, qj)); llR <- llp(.joint33(1, r_hat, qi, qj))
  expect_equal(res$logLik[1, 2], max(llC, llR),           tolerance = 1e-8)  # no missing -> pure kernel
  expect_equal(res$lod_ph[1, 2], abs(llC - llR) / log(10), tolerance = 1e-8)
})


# 10 --------------------------------------------------------------------------
test_that("results are deterministic with one thread", {
  RcppParallel::setThreadOptions(numThreads = 1)
  M <- list(A = c(m1 = 1L, m2 = 1L), B = c(m1 = 1L, m2 = 1L))
  gA <- rbind(matrix(rep(c(2L, 2L), 15), ncol = 2, byrow = TRUE),
              matrix(rep(c(0L, 0L), 12), ncol = 2, byrow = TRUE),
              matrix(rep(c(2L, 0L),  4), ncol = 2, byrow = TRUE))
  gB <- rbind(matrix(rep(c(2L, 0L), 10), ncol = 2, byrow = TRUE),
              matrix(rep(c(0L, 2L),  9), ncol = 2, byrow = TRUE))
  colnames(gA) <- colnames(gB) <- c("m1", "m2")
  storage.mode(gA) <- storage.mode(gB) <- "integer"
  r1 <- .pw(list(A = gA, B = gB), M, lambda = 20, q0 = 0.5)
  r2 <- .pw(list(A = gA, B = gB), M, lambda = 20, q0 = 0.5)
  expect_identical(r1$r,       r2$r)
  expect_identical(r1$lod_ph,  r2$lod_ph)
  expect_identical(r1$logLik,  r2$logLik)
  expect_identical(r1$q_list,  r2$q_list)
  expect_identical(r1$lod_ph_list, r2$lod_ph_list)
})


# Three deterministic worked examples -----------------------------------------
test_that("deterministic example A: contrasting paternal q among dams", {
  RcppParallel::setThreadOptions(numThreads = 1)
  M <- list(highQ = c(m1 = 1L, m2 = 1L), lowQ = c(m1 = 1L, m2 = 1L))
  gH <- rbind(matrix(rep(c(2L, 2L), 8), ncol = 2, byrow = TRUE),
              matrix(rep(c(0L, 0L), 2), ncol = 2, byrow = TRUE))   # q = 8/10 = 0.8
  gL <- rbind(matrix(rep(c(2L, 2L), 2), ncol = 2, byrow = TRUE),
              matrix(rep(c(0L, 0L), 8), ncol = 2, byrow = TRUE))   # q = 2/10 = 0.2
  colnames(gH) <- colnames(gL) <- c("m1", "m2")
  storage.mode(gH) <- storage.mode(gL) <- "integer"
  res <- .pw(list(highQ = gH, lowQ = gL), M, lambda = 0, q0 = 0.5)
  expect_equal(unname(res$q_list$highQ), c(0.8, 0.8), tolerance = 1e-9)
  expect_equal(unname(res$q_list$lowQ),  c(0.2, 0.2), tolerance = 1e-9)
})

test_that("deterministic example B: a partial observation enters the reported likelihood", {
  RcppParallel::setThreadOptions(numThreads = 1)
  M <- .dhet2()
  base <- rbind(matrix(rep(c(2L, 2L), 20), ncol = 2, byrow = TRUE),
                matrix(rep(c(0L, 0L), 20), ncol = 2, byrow = TRUE))
  colnames(base) <- c("m1", "m2"); storage.mode(base) <- "integer"
  r0 <- .pw(list(P1 = base), M, lambda = 0, q0 = 0.5)
  # Add a balanced pair of marker-2-only offspring (one AA, one aa). Its single-
  # marker marginal is invariant to r and phase, and being balanced it leaves q_m2
  # at 0.5, so ONLY the reported logLik changes -- isolating the marginal term.
  add <- matrix(c(NA_integer_, NA_integer_, 2L, 0L), ncol = 2); colnames(add) <- c("m1", "m2")
  r1 <- .pw(list(P1 = rbind(base, add)), M, lambda = 0, q0 = 0.5)
  expect_equal(r1$r[1, 2],         r0$r[1, 2],         tolerance = 1e-9)  # r unchanged
  expect_equal(r1$lod_ph[1, 2],    r0$lod_ph[1, 2],    tolerance = 1e-9)  # phase LOD unchanged
  expect_equal(r1$q_list$P1[[2]],  r0$q_list$P1[[2]],  tolerance = 1e-9)  # balanced -> q unchanged
  expect_false(isTRUE(all.equal(r1$logLik[1, 2], r0$logLik[1, 2])))       # logLik changed
  # the exact logLik shift equals the two marginal terms  log(q/2) + log((1-q)/2)
  q2 <- r0$q_list$P1[[2]]
  expect_equal(r1$logLik[1, 2] - r0$logLik[1, 2],
               log(q2 / 2) + log((1 - q2) / 2), tolerance = 1e-8)
})

test_that("deterministic example C: one strong dam and one weak dam", {
  RcppParallel::setThreadOptions(numThreads = 1)
  M <- list(strong = c(m1 = 1L, m2 = 1L), weak = c(m1 = 1L, m2 = 1L))
  gS <- rbind(matrix(rep(c(2L, 2L), 50), ncol = 2, byrow = TRUE),
              matrix(rep(c(0L, 0L), 50), ncol = 2, byrow = TRUE))
  # weak dam: near-balanced phase evidence -> small LOD, well below the strong dam
  gW <- .gm(c(2L, 2L, 0L, 2L, 0L,   2L, 2L, 0L, 0L, 2L))   # rows (2,2)(2,2)(0,0)(2,0)(0,2)
  colnames(gS) <- c("m1", "m2"); storage.mode(gS) <- "integer"
  res <- .pw(list(strong = gS, weak = gW), M, lambda = 0, q0 = 0.5)
  expect_gt(res$lod_ph_list$strong[1, 2], 25)
  expect_lt(res$lod_ph_list$weak[1, 2], 0.25 * res$lod_ph_list$strong[1, 2])
  expect_equal(res$lod_ph[1, 2],
               res$lod_ph_list$strong[1, 2] + res$lod_ph_list$weak[1, 2], tolerance = 1e-9)
})
