# Milestone 2: safe dam-specific phase assembly.
#   * dam-specific lod_ph_list edge weights (not pooled),
#   * connected components; unresolved relative phase across components,
#   * no-information graphs are NOT forced to all-coupling,
#   * whole-component sign-reversal invariance; greedy ascent never decreases obj.

# Build a minimal HSMap.tpt from per-dam phase-sign (A) and LOD (L) matrices.
.mk_tpt <- function(A_list, L_list, markers) {
  pooled <- Reduce(`+`, lapply(L_list, function(m) { m[is.na(m)] <- 0; m }))
  fit <- list(mom_phase_list = A_list, lod_ph_list = L_list,
              lod_ph = pooled, markers = markers)
  structure(list(fit = fit, markers = markers), class = "HSMap.tpt")
}
.sqm <- function(fun, markers) {
  n <- length(markers); m <- matrix(NA_real_, n, n, dimnames = list(markers, markers))
  for (i in 1:n) for (j in 1:n) if (i != j) m[i, j] <- fun(i, j)
  m
}

test_that("each dam uses its own lod_ph_list; a weak dam does not inherit another's edges", {
  mk <- paste0("m", 1:4)
  A_coup <- .sqm(function(i, j) 1, mk)                     # all coupling
  L_strong <- .sqm(function(i, j) 10, mk)                  # strong support (dam A)
  L_zero   <- .sqm(function(i, j) 0,  mk)                  # NO support (dam B)
  tpt <- .mk_tpt(list(A = A_coup, B = A_coup),
                 list(A = L_strong, B = L_zero), mk)
  # dam B has no supported edges -> a warning is expected
  expect_warning(res <- phase_from_pairwise(tpt, order = mk, dam = "all"), "No supported phase")
  # dam A: one component, fully resolved coupling
  expect_identical(res$A$n_components, 1L)
  expect_true(res$A$fully_resolved)
  expect_equal(res$A$phase_vec, c(1L, 1L, 1L))
  # dam B: zero support -> every marker isolated, NOTHING resolved (did NOT inherit A)
  expect_identical(res$B$n_components, 4L)
  expect_true(all(is.na(res$B$phase_vec)))
  expect_false(isTRUE(res$B$fully_resolved))
})

test_that("disconnected components are identified and cross-component phase is unresolved", {
  mk <- paste0("m", 1:4)
  A <- .sqm(function(i, j) 1, mk)                          # coupling everywhere it is called
  # support only within {m1,m2} and {m3,m4}; nothing across the m2|m3 gap
  L <- .sqm(function(i, j) if ((i <= 2 && j <= 2) || (i >= 3 && j >= 3)) 8 else 0, mk)
  tpt <- .mk_tpt(list(D = A), list(D = L), mk)
  res <- phase_from_pairwise(tpt, order = mk, dam = "D")
  expect_identical(res$n_components, 2L)
  expect_equal(sort(res$component_sizes), c(2L, 2L))
  # interval m1-m2 and m3-m4 resolved; m2-m3 (across components) unresolved
  expect_equal(res$resolved_interval, c(TRUE, FALSE, TRUE))
  expect_true(is.na(res$phase_vec[2]))
  expect_true(!is.na(res$phase_vec[1]) && !is.na(res$phase_vec[3]))
  expect_identical(res$unresolved_intervals, 2L)
})

test_that("an all-zero-evidence graph is NOT returned as resolved all-coupling", {
  mk <- paste0("m", 1:3)
  A <- .sqm(function(i, j) 1, mk)
  L <- .sqm(function(i, j) 0, mk)                          # no evidence anywhere
  tpt <- .mk_tpt(list(D = A), list(D = L), mk)
  expect_warning(res <- phase_from_pairwise(tpt, order = mk, dam = "D"), "No supported phase")
  expect_identical(res$n_components, 3L)                   # all isolated
  expect_true(all(is.na(res$phase_vec)))                   # NOT all-coupling
  expect_length(res$unresolved_markers, 3L)
})

test_that("resolved phase is invariant to a whole-component sign reversal (anchor flip)", {
  mk <- paste0("m", 1:4)
  A <- .sqm(function(i, j) if ((i + j) %% 2 == 0) 1 else 0, mk)   # mixed but connected
  L <- .sqm(function(i, j) 5, mk)
  tpt <- .mk_tpt(list(D = A), list(D = L), mk)
  r0 <- phase_from_pairwise(tpt, order = mk, dam = "D", anchor_label = 0L)
  r1 <- phase_from_pairwise(tpt, order = mk, dam = "D", anchor_label = 1L)
  # relative phase (phase_vec) invariant; absolute clusters flip within the component
  expect_equal(r0$phase_vec, r1$phase_vec)
  expect_equal(r0$objective, r1$objective, tolerance = 1e-9)
  expect_true(all(r0$clusters != r1$clusters) || all(r0$clusters == r1$clusters))
})

test_that("greedy coordinate ascent never decreases the within-component objective", {
  set.seed(1)  # deterministic construction (no engine RNG involved)
  n <- 6
  J <- matrix(0, n, n)
  vals <- c(2, -1, 3, -2, 1, 4, -3, 2, 1, -1, 2, 3, -2, 1, 2)
  k <- 1
  for (i in 1:(n - 1)) for (j in (i + 1):n) { J[i, j] <- vals[((k - 1) %% length(vals)) + 1]; J[j, i] <- J[i, j]; k <- k + 1 }
  x0 <- rep(1, n); x0[c(2, 4, 6)] <- -1
  obj0 <- sum(J[upper.tri(J)] * outer(x0, x0)[upper.tri(J)])
  gr <- HSMap:::.pf_greedy(J, x0, max_passes = 50L, tol = 1e-9)
  expect_gte(gr$objective, obj0)                           # never worse than the start
})

test_that("phase assembly is deterministic", {
  mk <- paste0("m", 1:4)
  A <- .sqm(function(i, j) if ((i + j) %% 2 == 0) 1 else 0, mk)
  L <- .sqm(function(i, j) 5, mk)
  tpt <- .mk_tpt(list(D = A), list(D = L), mk)
  r1 <- phase_from_pairwise(tpt, order = mk, dam = "D")
  r2 <- phase_from_pairwise(tpt, order = mk, dam = "D")
  expect_identical(r1$phase_vec, r2$phase_vec)
  expect_identical(r1$clusters, r2$clusters)
  expect_identical(r1$component, r2$component)
})
