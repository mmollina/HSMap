# Milestone 1: pairwise optimizer + API robustness.
#   * true no-linkage null at exactly r = 0.5 (not 0.49),
#   * grid + bounded local refinement (multimodal-safe),
#   * default threads = NULL, tol/maxit honored, deterministic, deprecations warn.

.pwo <- function(G, M, ...) HSMap:::pairwise_rf_estimation_multi_parallel_cpp(G, M, ...)
.cs  <- function(p, m) if (m == 0) c(1 - p, p, 0) else c(0, 1 - p, p)
.j33 <- function(phase, r, qi, qj) {
  P <- matrix(0, 3, 3)
  if (phase == 0) { no <- list(c(0, 0), c(1, 1)); re <- list(c(0, 1), c(1, 0)) }
  else            { no <- list(c(1, 0), c(0, 1)); re <- list(c(0, 0), c(1, 1)) }
  for (m in no) P <- P + 0.5 * (1 - r) * outer(.cs(qi, m[1]), .cs(qj, m[2]))
  for (m in re) P <- P + 0.5 * r       * outer(.cs(qi, m[1]), .cs(qj, m[2]))
  P
}
# reference single-dam objective at fixed q_i, q_j (no missing data)
.obj_ref <- function(C, qi, qj, r, tiny = 1e-12) {
  llp <- function(P) sum(ifelse(C > 0, C * log(pmax(P, tiny)), 0))
  max(llp(.j33(0, r, qi, qj)), llp(.j33(1, r, qi, qj)))
}
# expand a 3x3 count matrix to an offspring genotype matrix (cols m1, m2)
.expand <- function(C) {
  rows <- list()
  for (a in 0:2) for (b in 0:2) if (C[a + 1, b + 1] > 0)
    rows[[length(rows) + 1]] <- matrix(rep(c(a, b), C[a + 1, b + 1]), ncol = 2, byrow = TRUE)
  g <- do.call(rbind, rows); colnames(g) <- c("m1", "m2"); storage.mode(g) <- "integer"; g
}
.hsdata <- function(G_list, M_list) structure(list(G_list = G_list, M_list = M_list), class = "HSMap.data")


test_that("pairwise_rf works with default threads = NULL", {
  M <- list(P1 = c(m1 = 1L, m2 = 1L))
  g <- rbind(matrix(rep(c(2L, 2L), 15), ncol = 2, byrow = TRUE),
             matrix(rep(c(0L, 0L), 15), ncol = 2, byrow = TRUE))
  colnames(g) <- c("m1", "m2"); storage.mode(g) <- "integer"
  dat <- .hsdata(list(P1 = g), M)
  expect_silent(res <- pairwise_rf(dat))            # threads = NULL default
  expect_true(is.finite(res$fit$r[1, 2]))
  expect_identical(res$threads, 1L)
})

test_that("the no-linkage null is evaluated at exactly r = 0.5", {
  RcppParallel::setThreadOptions(numThreads = 1)
  M <- list(P1 = c(m1 = 1L, m2 = 1L))
  g <- rbind(matrix(rep(c(2L, 2L), 20), ncol = 2, byrow = TRUE),
             matrix(rep(c(0L, 0L), 20), ncol = 2, byrow = TRUE),
             matrix(rep(c(2L, 0L),  4), ncol = 2, byrow = TRUE))
  colnames(g) <- c("m1", "m2"); storage.mode(g) <- "integer"
  res <- .pwo(list(P1 = g), M, lambda = 0, q0 = 0.5)
  qi <- res$q_list$P1[[1]]; qj <- res$q_list$P1[[2]]
  C <- matrix(0, 3, 3); for (rr in seq_len(nrow(g))) C[g[rr, 1] + 1, g[rr, 2] + 1] <- C[g[rr, 1] + 1, g[rr, 2] + 1] + 1
  # LOD must use the null at EXACTLY 0.5, not 0.49
  lod_05  <- (res$logLik[1, 2] - .obj_ref(C, qi, qj, 0.50)) / log(10)
  lod_049 <- (res$logLik[1, 2] - .obj_ref(C, qi, qj, 0.49)) / log(10)
  expect_equal(res$lod_r[1, 2], lod_05, tolerance = 1e-6)
  expect_false(isTRUE(all.equal(res$lod_r[1, 2], lod_049)))
})

test_that("truly unlinked markers give r ~ 0.5 and LOD ~ 0", {
  RcppParallel::setThreadOptions(numThreads = 1)
  M <- list(P1 = c(m1 = 1L, m2 = 1L))
  # counts exactly proportional to the independent (r=0.5) joint at q=0.5:
  # marginal (0.25, 0.5, 0.25) outer product, scaled by 16 (then x3 for stability)
  base <- outer(c(1, 2, 1), c(1, 2, 1))            # = 16 * independent joint
  C <- base * 3
  g <- .expand(C)
  res <- .pwo(list(P1 = g), M, lambda = 0, q0 = 0.5)
  expect_gt(res$r[1, 2], 0.47)
  expect_equal(res$r[1, 2], 0.5, tolerance = 0.03)
  expect_lt(abs(res$lod_r[1, 2]), 0.05)            # near zero
  expect_identical(res$no_linkage[1, 2], 1L)       # flagged at the 0.5 boundary
})

test_that("grid + refinement matches a dense reference grid", {
  RcppParallel::setThreadOptions(numThreads = 1)
  M <- list(P1 = c(m1 = 1L, m2 = 1L))
  g <- rbind(matrix(rep(c(2L, 2L), 30), ncol = 2, byrow = TRUE),
             matrix(rep(c(0L, 0L), 25), ncol = 2, byrow = TRUE),
             matrix(rep(c(2L, 0L),  6), ncol = 2, byrow = TRUE),
             matrix(rep(c(0L, 2L),  5), ncol = 2, byrow = TRUE))
  colnames(g) <- c("m1", "m2"); storage.mode(g) <- "integer"
  res <- .pwo(list(P1 = g), M, lambda = 0, q0 = 0.5, tol = 1e-8)
  qi <- res$q_list$P1[[1]]; qj <- res$q_list$P1[[2]]
  C <- matrix(0, 3, 3); for (rr in seq_len(nrow(g))) C[g[rr, 1] + 1, g[rr, 2] + 1] <- C[g[rr, 1] + 1, g[rr, 2] + 1] + 1
  grid <- seq(1e-6, 0.5, length.out = 5001)
  fvals <- vapply(grid, function(r) .obj_ref(C, qi, qj, r), numeric(1))
  ref_r <- grid[which.max(fvals)]; ref_f <- max(fvals)
  expect_equal(res$r[1, 2], ref_r, tolerance = 2e-3)            # same optimum location
  expect_gte(res$logLik[1, 2], ref_f - 1e-6)                    # never worse than the dense grid
})

test_that("tol and maxit are actually used by the refinement", {
  RcppParallel::setThreadOptions(numThreads = 1)
  M <- list(P1 = c(m1 = 1L, m2 = 1L))
  g <- rbind(matrix(rep(c(2L, 2L), 40), ncol = 2, byrow = TRUE),
             matrix(rep(c(0L, 0L), 33), ncol = 2, byrow = TRUE),
             matrix(rep(c(2L, 0L),  7), ncol = 2, byrow = TRUE))
  colnames(g) <- c("m1", "m2"); storage.mode(g) <- "integer"
  coarse <- .pwo(list(P1 = g), M, lambda = 0, tol = 0.2,  maxit = 200, return_diagnostics = TRUE)
  fine   <- .pwo(list(P1 = g), M, lambda = 0, tol = 1e-9, maxit = 200, return_diagnostics = TRUE)
  # finer tol refines further: more objective evaluations and a >= objective
  expect_gt(fine$diagnostics$n_eval[1, 2], coarse$diagnostics$n_eval[1, 2])
  expect_gte(fine$logLik[1, 2], coarse$logLik[1, 2] - 1e-9)
  # maxit caps refinement iterations: maxit=1 evaluates fewer than maxit=200
  m1  <- .pwo(list(P1 = g), M, lambda = 0, tol = 1e-12, maxit = 1,   return_diagnostics = TRUE)
  m200<- .pwo(list(P1 = g), M, lambda = 0, tol = 1e-12, maxit = 200, return_diagnostics = TRUE)
  expect_lt(m1$diagnostics$n_eval[1, 2], m200$diagnostics$n_eval[1, 2])
})

test_that("pairwise optimizer is deterministic with one thread", {
  RcppParallel::setThreadOptions(numThreads = 1)
  M <- list(A = c(m1 = 1L, m2 = 1L), B = c(m1 = 1L, m2 = 1L))
  gA <- rbind(matrix(rep(c(2L, 2L), 15), ncol = 2, byrow = TRUE),
              matrix(rep(c(0L, 0L), 12), ncol = 2, byrow = TRUE),
              matrix(rep(c(2L, 0L),  4), ncol = 2, byrow = TRUE))
  gB <- rbind(matrix(rep(c(2L, 0L), 10), ncol = 2, byrow = TRUE),
              matrix(rep(c(0L, 2L),  9), ncol = 2, byrow = TRUE))
  colnames(gA) <- colnames(gB) <- c("m1", "m2")
  storage.mode(gA) <- storage.mode(gB) <- "integer"
  r1 <- .pwo(list(A = gA, B = gB), M, lambda = 20)
  r2 <- .pwo(list(A = gA, B = gB), M, lambda = 20)
  expect_identical(r1$r, r2$r); expect_identical(r1$lod_ph, r2$lod_ph)
  expect_identical(r1$no_linkage, r2$no_linkage)
})

test_that("compatibility output names remain present", {
  RcppParallel::setThreadOptions(numThreads = 1)
  M <- list(A = c(m1 = 1L, m2 = 1L))
  g <- rbind(matrix(rep(c(2L, 2L), 10), ncol = 2, byrow = TRUE),
             matrix(rep(c(0L, 0L), 10), ncol = 2, byrow = TRUE))
  colnames(g) <- c("m1", "m2"); storage.mode(g) <- "integer"
  res <- .pwo(list(A = g), M)
  expect_true(all(c("r", "lod_r", "lod_ph", "logLik", "mom_phase_list") %in% names(res)))
  expect_true(all(c("lod_ph_list", "q_list", "no_linkage", "optimizer", "n_grid") %in% names(res)))
})

test_that("deprecated pairwise arguments warn rather than silently do nothing", {
  M <- list(P1 = c(m1 = 1L, m2 = 1L))
  g <- rbind(matrix(rep(c(2L, 2L), 8), ncol = 2, byrow = TRUE),
             matrix(rep(c(0L, 0L), 8), ncol = 2, byrow = TRUE))
  colnames(g) <- c("m1", "m2"); storage.mode(g) <- "integer"
  dat <- .hsdata(list(P1 = g), M)
  expect_warning(pairwise_rf(dat, r_start = 0.1), "r_start")
  expect_warning(pairwise_rf(dat, share_pi_across_dams = TRUE), "share_pi_across_dams")
  expect_error(pairwise_rf(dat, totally_unknown_arg = 3), "unused argument")
})
