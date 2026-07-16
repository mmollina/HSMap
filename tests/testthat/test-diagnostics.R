# Commit 4: pairwise no_linkage flag independent of optimizer tol;
# phase support metadata distinguishing direct-edge vs path-based resolution.

test_that("the no_linkage flag is independent of the optimizer tol", {
  RcppParallel::setThreadOptions(numThreads = 1)
  M <- list(P = c(m1 = 1L, m2 = 1L))
  # moderate recombination (~0.3): coupling-leaning with ~30% recombinants
  g <- rbind(matrix(rep(c(2L, 2L), 35), ncol = 2, byrow = TRUE),
             matrix(rep(c(0L, 0L), 35), ncol = 2, byrow = TRUE),
             matrix(rep(c(2L, 0L), 15), ncol = 2, byrow = TRUE),
             matrix(rep(c(0L, 2L), 15), ncol = 2, byrow = TRUE))
  colnames(g) <- c("m1", "m2"); storage.mode(g) <- "integer"
  fine <- HSMap:::pairwise_rf_estimation_multi_parallel_cpp(list(P = g), M, lambda = 0, tol = 1e-8)
  big  <- HSMap:::pairwise_rf_estimation_multi_parallel_cpp(list(P = g), M, lambda = 0, tol = 0.3)
  # r is clearly well below 0.5 under both tolerances
  expect_gt(big$r[1, 2], 0.15); expect_lt(big$r[1, 2], 0.45)
  # a huge optimizer tol must NOT flag an r substantially below 0.5 as no linkage
  expect_identical(big$no_linkage[1, 2], 0L)
  expect_identical(fine$no_linkage[1, 2], 0L)
})

test_that("no_linkage still flags a genuine r = 0.5 regardless of tol", {
  RcppParallel::setThreadOptions(numThreads = 1)
  M <- list(P = c(m1 = 1L, m2 = 1L))
  C <- outer(c(1, 2, 1), c(1, 2, 1)) * 3          # independent (r=0.5) joint at q=0.5
  rows <- list()
  for (a in 0:2) for (b in 0:2) if (C[a + 1, b + 1] > 0)
    rows[[length(rows) + 1]] <- matrix(rep(c(a, b), C[a + 1, b + 1]), ncol = 2, byrow = TRUE)
  g <- do.call(rbind, rows); colnames(g) <- c("m1", "m2"); storage.mode(g) <- "integer"
  for (tl in c(1e-8, 0.3)) {
    res <- HSMap:::pairwise_rf_estimation_multi_parallel_cpp(list(P = g), M, lambda = 0, tol = tl)
    expect_identical(res$no_linkage[1, 2], 1L)
  }
})

# ---- direct vs path-based phase resolution -----------------------------------
.mk_tpt2 <- function(A_list, L_list, markers) {
  pooled <- Reduce(`+`, lapply(L_list, function(m) { m[is.na(m)] <- 0; m }))
  fit <- list(mom_phase_list = A_list, lod_ph_list = L_list, lod_ph = pooled, markers = markers)
  structure(list(fit = fit, markers = markers), class = "HSMap.tpt")
}

test_that("phase metadata distinguishes direct-edge from path-based resolution", {
  mk <- paste0("m", 1:4)
  A <- matrix(NA_real_, 4, 4, dimnames = list(mk, mk))
  L <- matrix(0,       4, 4, dimnames = list(mk, mk))
  # supported edges m1-m2 and m1-m3 (coupling, LOD 5); the adjacent edge m2-m3 has a
  # phase sign but ZERO LOD (unsupported). m4 is isolated.
  A[1, 2] <- A[2, 1] <- 1; L[1, 2] <- L[2, 1] <- 5
  A[1, 3] <- A[3, 1] <- 1; L[1, 3] <- L[3, 1] <- 5
  A[2, 3] <- A[3, 2] <- 1                                  # direct m2-m3 sign present, LOD 0
  tpt <- .mk_tpt2(list(D = A), list(D = L), mk)
  ph <- phase_from_pairwise(tpt, order = mk, dam = "D")

  # interval m1-m2: direct supported edge -> "direct"
  # interval m2-m3: resolved via the path through m1 (direct LOD is 0) -> "path"
  # interval m3-m4: m4 isolated -> "unresolved"
  expect_identical(ph$resolved_via, c("direct", "path", "unresolved"))
  expect_identical(ph$direct_edge,  c(TRUE, FALSE, FALSE))
  # the path-resolved interval has zero/NA direct adjacent support but IS resolved
  expect_true(is.na(ph$interval_support[2]) || ph$interval_support[2] == 0)
  expect_true(ph$resolved_interval[2])
  # a zero direct LOD is not reported as "resolved by a direct edge"
  expect_false(ph$direct_edge[2])
})
