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

test_that("phase-support fields separate raw direct LOD, direct support, and resolution route", {
  mk <- paste0("m", 1:5)
  A <- matrix(NA_real_, 5, 5, dimnames = list(mk, mk))
  L <- matrix(0,        5, 5, dimnames = list(mk, mk))
  # a coupling star centred on m1 supports m2, m3, m4 (LOD 5, threshold 3). The adjacent
  # edges m2-m3 (LOD 0) and m3-m4 (LOD 0.5, below threshold) carry a phase sign but are
  # NOT supported, so those intervals resolve only through the star. m5 is isolated.
  A[1, 2] <- A[2, 1] <- 1; L[1, 2] <- L[2, 1] <- 5          # supported
  A[1, 3] <- A[3, 1] <- 1; L[1, 3] <- L[3, 1] <- 5          # supported
  A[1, 4] <- A[4, 1] <- 1; L[1, 4] <- L[4, 1] <- 5          # supported
  A[2, 3] <- A[3, 2] <- 1                                    # adjacent sign, raw LOD 0
  A[3, 4] <- A[4, 3] <- 1; L[3, 4] <- L[4, 3] <- 0.5        # adjacent sign, raw LOD 0.5 (< 3)
  tpt <- .mk_tpt2(list(D = A), list(D = L), mk)
  ph <- phase_from_pairwise(tpt, order = mk, dam = "D", min_phase_lod = 3)

  # m1-m2 directly supported; m2-m3 path (zero direct LOD); m3-m4 path (raw LOD below
  # threshold); m4-m5 unresolved (m5 isolated).
  expect_identical(ph$resolved_via,     c("direct", "path", "path", "unresolved"))
  expect_identical(ph$direct_supported, c(TRUE, FALSE, FALSE, FALSE))
  expect_equal(ph$direct_lod, c(5, 0, 0.5, 0))              # RAW adjacent LOD, reported as-is
  expect_identical(ph$resolved_interval, c(TRUE, TRUE, TRUE, FALSE))

  # directly supported edge: raw LOD passes and phase is resolved by that edge
  expect_true(ph$direct_supported[1] && ph$resolved_via[1] == "direct")
  # zero direct LOD, resolved via a path -> not described as directly supported
  expect_true(ph$direct_lod[2] == 0 && !ph$direct_supported[2] && ph$resolved_interval[2])
  # raw direct LOD BELOW threshold -> reported (0.5) but not counted as support
  expect_true(ph$direct_lod[3] == 0.5 && !ph$direct_supported[3] && ph$resolved_via[3] == "path")
  # unresolved interval
  expect_true(!ph$resolved_interval[4] && ph$resolved_via[4] == "unresolved")
  # interval_support is a strict alias of the raw direct LOD
  expect_identical(ph$interval_support, ph$direct_lod)
})
