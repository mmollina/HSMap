# Local order diagnostics: try_marker() and ripple_map().
#
# These validate the primitives on simulated data where the true order is
# known. Likelihood ties are expected by design, so tests assert membership of
# near-optimal sets and rank behaviour rather than exact floating-point values.

LAM <- 2      # test setting only; nothing here depends on a production default

## A well-spaced chromosome: order is identifiable at this information content.
sim_chr <- function(seed, Tm = 24L, n = 300L, r = 0.08, err = 0.01, miss = 0) {
  set.seed(seed)
  r_true <- rep(r, Tm - 1L)
  sim <- sim_multi_pop(T_markers = Tm, n_pops = 1, n_ind_per_pop = n,
                       marker_intersection = 1, r_vec = r_true,
                       phase_mode = "random", repulsion_rate = 0.3,
                       maternal_geno_mode = "all_het", maternal_pA = 0.5,
                       paternal_pA_base = 0.4, error_rate = err,
                       seed = seed, miss_rate = miss)
  mk <- sim$truth$markers_union
  list(dat = make_dat(sim),
       oph = oracle_phased(mk, 1L - sim$truth$v_true[[1]], "P1"),
       mk = mk, r_true = r_true)
}

## Rebuild an HSMap.phased for a permuted order, carrying phase correctly.
reorder_phased <- function(oph, new_order) {
  h <- HSMap:::.hsmap_marker_phase(oph$order, oph$phase_vec)
  oracle_phased(new_order, HSMap:::.hsmap_phase_for(new_order, h), oph$dam)
}

# ---------------------------------------------------------------------------
# phase handling
# ---------------------------------------------------------------------------
test_that("marker-level phase reconstruction round-trips exactly", {
  s <- sim_chr(101L, Tm = 20L, n = 120L)
  h <- HSMap:::.hsmap_marker_phase(s$oph$order, s$oph$phase_vec)
  expect_identical(HSMap:::.hsmap_phase_for(s$oph$order, h),
                   as.integer(s$oph$phase_vec))
  ## a permutation gets its own phase, and restoring the order restores it
  perm <- s$oph$order[c(1:4, 6L, 5L, 7:20)]
  pv_perm <- HSMap:::.hsmap_phase_for(perm, h)
  expect_length(pv_perm, length(perm) - 1L)
  expect_identical(HSMap:::.hsmap_phase_for(s$oph$order, h),
                   as.integer(s$oph$phase_vec))
  ## Swapping two markers changes the interval phase only when their marker
  ## states differ: if h[i] == h[j] the swap is phase-invariant, which is
  ## correct behaviour, not a bug. Pick a pair that does differ.
  j <- which(h[-1L] != h[-length(h)])[1L]      # h[j] != h[j+1]
  perm2 <- s$oph$order
  perm2[c(j, j + 1L)] <- perm2[c(j + 1L, j)]
  expect_false(identical(HSMap:::.hsmap_phase_for(perm2, h),
                         as.integer(s$oph$phase_vec)))
})

test_that("unresolved phase is rejected with a clear message", {
  s <- sim_chr(102L, Tm = 12L, n = 80L)
  bad <- s$oph; bad$phase_vec[3L] <- NA_integer_
  expect_error(try_marker(s$dat, bad, marker = s$mk[6L], window = 7L),
               "fully resolved phase")
})

# ---------------------------------------------------------------------------
# try_marker: placement recovery
# ---------------------------------------------------------------------------
test_that("try_marker puts a correctly placed marker at or near its position", {
  s <- sim_chr(201L, Tm = 24L, n = 300L)
  tr <- try_marker(s$dat, s$oph, marker = s$mk[12L], window = 11L,
                   lambda = LAM, epsilon = 0.01)
  expect_s3_class(tr, "HSMap.try")
  expect_equal(nrow(tr$profile), length(tr$framework) + 1L)
  expect_true(tr$current_pos %in% tr$near_optimal)
  expect_lte(abs(tr$displacement), 1L)
})

test_that("try_marker recovers a deliberately displaced marker", {
  s <- sim_chr(202L, Tm = 24L, n = 300L)
  ## move marker 12 three positions away in the SUPPLIED order
  ord <- s$mk
  moved <- ord[12L]
  wrong <- append(ord[-12L], moved, after = 14L)
  ph <- reorder_phased(s$oph, wrong)
  tr <- try_marker(s$dat, ph, marker = moved, window = 15L,
                   lambda = LAM, epsilon = 0.01)
  ## true neighbours are the markers flanking it in the TRUE order
  true_left <- ord[11L]
  pos_true <- match(true_left, tr$framework) + 1L
  expect_lte(abs(tr$best_pos - pos_true), 1L)
  ## and the true neighbourhood beats the supplied wrong placement
  expect_gt(tr$profile$logLik[tr$best_pos], tr$profile$logLik[tr$current_pos])
})

test_that("try_marker reports ambiguity for co-segregating markers", {
  ## three markers at r = 0 with each other: position among them is not
  ## identifiable and the tool must say so rather than pick one
  set.seed(303)
  Tm <- 15L
  r_true <- rep(0.10, Tm - 1L)
  r_true[7:8] <- 0
  sim <- sim_multi_pop(T_markers = Tm, n_pops = 1, n_ind_per_pop = 300,
                       marker_intersection = 1, r_vec = r_true,
                       phase_mode = "random", repulsion_rate = 0.3,
                       maternal_geno_mode = "all_het", maternal_pA = 0.5,
                       paternal_pA_base = 0.4, error_rate = 0, seed = 303)
  mk <- sim$truth$markers_union
  dat <- make_dat(sim)
  oph <- oracle_phased(mk, 1L - sim$truth$v_true[[1]], "P1")
  tr <- try_marker(dat, oph, marker = mk[8L], window = 11L,
                   lambda = LAM, epsilon = 0.01)
  expect_false(tr$unique_optimum)
  expect_gte(length(tr$near_optimal), 2L)
  expect_true(tr$current_pos %in% tr$near_optimal)
})

test_that("try_marker reports no credible internal placement for an unlinked marker", {
  s <- sim_chr(404L, Tm = 20L, n = 300L)
  other <- sim_chr(999L, Tm = 4L, n = 300L)     # independent chromosome
  ## splice an unlinked marker into the data and the order
  G <- s$dat$G_list[[1]]
  ug <- other$dat$G_list[[1]][, 2L, drop = FALSE]
  colnames(ug) <- "UNLINKED"
  G2 <- cbind(G, ug)
  M2 <- c(s$dat$M_list[[1]], stats::setNames(1L, "UNLINKED"))
  dat2 <- structure(list(G_list = list(P1 = G2), M_list = list(P1 = M2)),
                    class = "HSMap.data")
  ord <- append(s$mk, "UNLINKED", after = 10L)
  h <- HSMap:::.hsmap_marker_phase(s$oph$order, s$oph$phase_vec)
  h["UNLINKED"] <- 1L
  ph <- oracle_phased(ord, HSMap:::.hsmap_phase_for(ord, h), "P1")
  tr <- try_marker(dat2, ph, marker = "UNLINKED", window = 13L,
                   lambda = LAM, epsilon = 0.01)
  ## Empirical signature (measured, not assumed): an unlinked marker is
  ## strongly EDGE-favouring with both edges tied, because an internal
  ## insertion severs the chain into two no-linkage boundaries while an edge
  ## placement leaves the chromosome intact.
  expect_identical(tr$placement, "no credible internal placement (unlinked)")
  expect_true(tr$best_is_edge)
  expect_false(tr$unique_optimum)
  expect_setequal(tr$near_optimal, c(1L, nrow(tr$profile)))
  ## the decisive evidence: at its best placement the marker is unlinked to
  ## every neighbour it has
  expect_gte(tr$min_flank_r_at_best, 0.4)
  ## every interior placement is clearly worse than the edges
  int <- !tr$profile$at_edge
  expect_lt(max(tr$profile$dLL_best[int]), -10)
  expect_true(all(is.finite(tr$profile$local_cM)))   # never a huge cM value
  ## a linked control on the same data resolves to a unique interior optimum
  tl <- try_marker(dat2, ph, marker = s$mk[6L], window = 13L,
                   lambda = LAM, epsilon = 0.01)
  expect_identical(tl$placement, "resolved")
  expect_false(tl$best_is_edge)
})

test_that("try_marker handles missing data and an absent marker", {
  s <- sim_chr(505L, Tm = 20L, n = 250L, miss = 0.10)
  tr <- try_marker(s$dat, s$oph, marker = s$mk[10L], window = 9L,
                   lambda = LAM, epsilon = 0.01)
  expect_true(all(tr$profile$converged))
  ## marker absent from the supplied order: both orientations scored
  drop <- s$mk[10L]
  ord <- setdiff(s$mk, drop)
  ph <- reorder_phased(s$oph, ord)
  tr2 <- try_marker(s$dat, ph, marker = drop, framework = ord[5:14],
                    lambda = LAM, epsilon = 0.01)
  expect_true(is.na(tr2$current_pos))
  expect_true(all(tr2$profile$orientation %in% c(0L, 1L)))
  expect_equal(nrow(tr2$profile), 11L)
})

# ---------------------------------------------------------------------------
# ripple_map: order support
# ---------------------------------------------------------------------------
test_that("ripple_map does not invent changes when the order is correct", {
  s <- sim_chr(601L, Tm = 22L, n = 300L)
  rp <- ripple_map(s$dat, s$oph, window = 4L, anchor = 6L, by = 4L,
                   lambda = LAM, epsilon = 0.01)
  expect_s3_class(rp, "HSMap.ripple")
  expect_true(all(rp$windows$converged))
  ## no window may claim a supported reordering of an already-correct order
  expect_equal(sum(rp$windows$supported), 0L)
})

test_that("ripple_map recovers a deliberate adjacent swap", {
  s <- sim_chr(602L, Tm = 22L, n = 300L)
  ord <- s$mk; ord[10:11] <- ord[11:10]
  ph <- reorder_phased(s$oph, ord)
  rp <- ripple_map(s$dat, ph, window = 4L, anchor = 8L, positions = 9L,
                   lambda = LAM, epsilon = 0.01)
  expect_true(rp$windows$supported)
  expect_identical(rp$best_orders[[1]], s$mk[9:12])
})

test_that("ripple_map recovers a short local reversal", {
  s <- sim_chr(603L, Tm = 22L, n = 300L)
  ord <- s$mk; ord[9:12] <- rev(ord[9:12])
  ph <- reorder_phased(s$oph, ord)
  rp <- ripple_map(s$dat, ph, window = 4L, anchor = 8L, positions = 9L,
                   lambda = LAM, epsilon = 0.01)
  ## the true order (or its reverse, which is likelihood-equivalent locally)
  expect_true(rp$windows$supported)
  best <- rp$best_orders[[1]]
  expect_true(identical(best, s$mk[9:12]) || identical(best, rev(s$mk[9:12])))
})

test_that("ripple_map reports ties instead of a false unique order", {
  set.seed(704)
  Tm <- 16L
  r_true <- rep(0.10, Tm - 1L)
  r_true[8:9] <- 0                        # three co-segregating markers
  sim <- sim_multi_pop(T_markers = Tm, n_pops = 1, n_ind_per_pop = 300,
                       marker_intersection = 1, r_vec = r_true,
                       phase_mode = "random", repulsion_rate = 0.3,
                       maternal_geno_mode = "all_het", maternal_pA = 0.5,
                       paternal_pA_base = 0.4, error_rate = 0, seed = 704)
  mk <- sim$truth$markers_union
  dat <- make_dat(sim)
  oph <- oracle_phased(mk, 1L - sim$truth$v_true[[1]], "P1")
  rp <- ripple_map(dat, oph, window = 4L, anchor = 6L, positions = 7L,
                   lambda = LAM, epsilon = 0.01)
  expect_gt(rp$windows$n_near_optimal, 1L)
})

test_that("no-linkage intervals never become large cM values", {
  set.seed(805)
  Tm <- 14L
  r_true <- rep(0.05, Tm - 1L)
  r_true[7L] <- 0.5                       # a true no-linkage boundary
  sim <- sim_multi_pop(T_markers = Tm, n_pops = 1, n_ind_per_pop = 250,
                       marker_intersection = 1, r_vec = r_true,
                       phase_mode = "random", repulsion_rate = 0.3,
                       maternal_geno_mode = "all_het", maternal_pA = 0.5,
                       paternal_pA_base = 0.4, error_rate = 0, seed = 805)
  mk <- sim$truth$markers_union
  dat <- make_dat(sim)
  oph <- oracle_phased(mk, 1L - sim$truth$v_true[[1]], "P1")
  tr <- try_marker(dat, oph, marker = mk[7L], window = 11L,
                   lambda = LAM, epsilon = 0.01, no_linkage_r = 0.4)
  expect_true(all(is.finite(tr$profile$local_cM)))
  expect_true(all(tr$profile$local_cM < 400))
  expect_true(any(tr$profile$n_no_linkage > 0L))
})

# ---------------------------------------------------------------------------
# contract
# ---------------------------------------------------------------------------
test_that("results are deterministic under fixed settings", {
  s <- sim_chr(901L, Tm = 16L, n = 200L)
  a <- try_marker(s$dat, s$oph, marker = s$mk[8L], window = 9L, lambda = LAM)
  b <- try_marker(s$dat, s$oph, marker = s$mk[8L], window = 9L, lambda = LAM)
  expect_equal(a$profile$logLik, b$profile$logLik)
  expect_identical(a$best_pos, b$best_pos)
  r1 <- ripple_map(s$dat, s$oph, window = 4L, anchor = 5L, positions = 6L,
                   lambda = LAM)
  r2 <- ripple_map(s$dat, s$oph, window = 4L, anchor = 5L, positions = 6L,
                   lambda = LAM)
  expect_equal(r1$windows$best_logLik, r2$windows$best_logLik)
})

test_that("ripple_map never selects on map length", {
  s <- sim_chr(902L, Tm = 20L, n = 300L)
  ord <- s$mk; ord[10:11] <- ord[11:10]
  ph <- reorder_phased(s$oph, ord)
  rp <- ripple_map(s$dat, ph, window = 4L, anchor = 8L, positions = 9L,
                   lambda = LAM)
  ## the chosen order is the likelihood optimum whatever its length
  expect_gt(rp$windows$best_logLik, rp$windows$baseline_logLik)
  expect_true(is.finite(rp$windows$dcM))
})

test_that("apply mode is off by default and explicit when used", {
  s <- sim_chr(903L, Tm = 18L, n = 250L)
  ord <- s$mk; ord[8:9] <- ord[9:8]
  ph <- reorder_phased(s$oph, ord)
  rp <- ripple_map(s$dat, ph, window = 4L, anchor = 6L, positions = 7L,
                   lambda = LAM)
  expect_null(rp$applied_order)
  rp2 <- ripple_map(s$dat, ph, window = 4L, anchor = 6L, positions = 7L,
                    lambda = LAM, apply = TRUE)
  expect_type(rp2$applied_order, "character")
  expect_setequal(rp2$applied_order, ord)          # a permutation, nothing lost
  expect_length(rp2$applied_order, length(ord))
})

test_that("input validation", {
  s <- sim_chr(904L, Tm = 12L, n = 100L)
  expect_error(try_marker(list(), s$oph, marker = s$mk[5L]), "HSMap.data")
  expect_error(try_marker(s$dat, s$oph, marker = "NOPE"), "not present")
  expect_error(try_marker(s$dat, s$oph, marker = s$mk[5L], tie_tol = -1),
               "non-negative")
  expect_error(try_marker(s$dat, s$oph, marker = s$mk[5L], r_starts = c(0.6)),
               "must be in")
  expect_error(ripple_map(s$dat, s$oph, window = 9L), "refusing")
  expect_error(ripple_map(s$dat, s$oph, window = 4L, positions = 999L),
               "out of range")
  expect_error(try_marker(s$dat, "nope", marker = s$mk[5L]), "HSMap.phased")
})

test_that("print and plot methods work", {
  s <- sim_chr(905L, Tm = 14L, n = 150L)
  tr <- try_marker(s$dat, s$oph, marker = s$mk[7L], window = 9L, lambda = LAM)
  expect_output(print(tr), "HSMap.try")
  expect_output(print(tr), "Placement")
  rp <- ripple_map(s$dat, s$oph, window = 4L, anchor = 5L, positions = 5L,
                   lambda = LAM)
  expect_output(print(rp), "HSMap.ripple")
  pf <- tempfile(fileext = ".pdf"); grDevices::pdf(pf)
  expect_silent(plot(tr)); grDevices::dev.off(); unlink(pf)
})
