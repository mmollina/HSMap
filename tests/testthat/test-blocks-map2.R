# Commit 2: blockwise interval-status classification, gap_r, get_block_map,
# block-aware plotting, and joint multi-dam blockwise fitting.

.na_at <- function(v, i) { v[i] <- NA_integer_; as.integer(v) }

# ---- interval-status classification + get_block_map (mock, deterministic) ----
.mock_blocks <- function(gap_r = 0.499) {
  mk <- paste0("m", 1:6)
  f1 <- structure(list(order = mk[1:3], fit = list(r = c(0.05, 0.49999))), class = "HSMap.map")
  f2 <- structure(list(order = mk[4:6], fit = list(r = c(0.08, 0.03))),    class = "HSMap.map")
  blocks <- list(list(block = 1, markers = mk[1:3], contributing_dams = "A", fit = f1),
                 list(block = 2, markers = mk[4:6], contributing_dams = "A", fit = f2))
  block_id <- c(1L, 1L, 1L, 2L, 2L, 2L)
  pv <- c(1L, 1L, NA_integer_, 1L, 1L)                 # interval 3 is the block boundary
  itab <- HSMap:::.block_interval_table(mk, pv, block_id, blocks, gap_r = gap_r,
                                        boundary_status = "unresolved_phase")
  structure(list(blocks = blocks, block_id = stats::setNames(block_id, mk),
                 unresolved_boundaries = 3L, n_blocks = 2L, order = mk,
                 interval_table = itab, gap_r = gap_r, dams = "A"),
            class = "HSMap.map.blocks")
}

test_that("interval_table classifies every interval, not linked-just-because-fitted", {
  mb <- .mock_blocks()
  st <- mb$interval_table$status
  expect_identical(st, c("linked", "no_linkage_boundary", "unresolved_phase", "linked", "linked"))
  # the fitted r ~ 0.5 within block 1 is a gap, NOT linked
  expect_identical(st[2], "no_linkage_boundary")
})

test_that("get_block_map gives safe positions, lengths, and gaps (no huge cM)", {
  mb <- .mock_blocks()
  bm <- get_block_map(mb, "haldane")
  # positions reset per block: block 1 starts at 0, block 2 starts at 0
  p <- bm$positions
  expect_equal(p$pos[p$marker == "m1"], 0)
  expect_equal(p$pos[p$marker == "m4"], 0)               # block 2 reset
  # the r~0.5 interval contributes NO distance; total is finite and small
  expect_true(is.na(bm$interval_table$dist_cM[2]))
  expect_lt(bm$total_linked_length, 30)                  # never ~1577 cM
  expect_identical(bm$n_gaps, 2L)                        # r~0.5 + the phase boundary
  expect_true(all(c(2L, 3L) %in% bm$gap_intervals))
})

test_that("gap_r controls the no-linkage classification", {
  mb_hi <- .mock_blocks(gap_r = 0.499)   # r=0.49999 >= 0.499 -> gap
  mb_lo <- .mock_blocks(gap_r = 0.5)     # r=0.49999 <  0.5   -> linked
  expect_identical(mb_hi$interval_table$status[2], "no_linkage_boundary")
  expect_identical(mb_lo$interval_table$status[2], "linked")
})

test_that("plot_block_map returns a ggplot with one panel per block", {
  skip_if_not_installed("ggplot2")
  mb <- .mock_blocks()
  p <- plot_block_map(mb, "haldane")
  expect_s3_class(p, "ggplot")
})

# ---- joint multi-dam blockwise (real fits) -----------------------------------
.two_dam <- function(seed = 31L, Tm = 8L, n = 120L) {
  set.seed(seed)
  sim <- sim_multi_pop(T_markers = Tm, n_pops = 2, n_ind_per_pop = c(n, n),
                       marker_intersection = 1, r_vec = rep(0.06, Tm - 1),
                       phase_mode = "random", repulsion_rate = 0.3,
                       maternal_geno_mode = "all_het", paternal_pA_base = 0.4,
                       error_rate = 0.01, seed = seed)
  dat <- make_dat(sim)
  names(dat$G_list) <- names(dat$M_list) <- c("A", "B")   # match the phased dam names
  list(dat = dat, mk = sim$truth$markers_union,
       pvA = as.integer(1L - sim$truth$v_true[[1]]),
       pvB = as.integer(1L - sim$truth$v_true[[2]]))
}
.phased <- function(dam, mk, pv, comp)
  structure(list(dam = dam, order = mk, phase_vec = as.integer(pv), component = as.integer(comp)),
            class = "HSMap.phased")

test_that("joint blockwise: disconnected components split at the common boundary", {
  RcppParallel::setThreadOptions(numThreads = 1)
  d <- .two_dam(); mk <- d$mk
  phA <- .phased("A", mk, .na_at(d$pvA, 4L), c(1,1,1,1,2,2,2,2))  # A unresolved at interval 4
  phB <- .phased("B", mk, d$pvB,             rep(1L, 8))          # B fully resolved
  pm <- structure(list(A = phA, B = phB), class = "HSMap.phased.multi")
  mb <- hmm_map_blocks(d$dat, pm, epsilon = 0.01, paternal_mode = "HWE")
  expect_identical(mb$n_blocks, 2L)
  expect_true(4L %in% mb$unresolved_boundaries)
  expect_setequal(mb$blocks[[1]]$contributing_dams, c("A", "B"))
  expect_setequal(mb$blocks[[2]]$contributing_dams, c("A", "B"))
})

test_that("joint blockwise: a dam is excluded from a block where its phase is unresolved", {
  RcppParallel::setThreadOptions(numThreads = 1)
  d <- .two_dam(); mk <- d$mk
  phA <- .phased("A", mk, d$pvA, rep(1L, 8))                       # A resolved everywhere
  # B: markers 4,5 are singletons; its only NA (interval 4) touches two singletons ->
  # NON-informative there, so it does NOT force a boundary, but B is unresolved in the block.
  phB <- .phased("B", mk, .na_at(d$pvB, 4L), c(1,1,1,2,3,4,4,4))
  pm <- structure(list(A = phA, B = phB), class = "HSMap.phased.multi")
  mb <- hmm_map_blocks(d$dat, pm, epsilon = 0.01, paternal_mode = "HWE")
  expect_identical(mb$n_blocks, 1L)                                # conservative rule: no split
  expect_setequal(mb$blocks[[1]]$contributing_dams, "A")           # B excluded (unresolved in block)
})

test_that("fully resolved multi-dam: ordinary joint == blockwise", {
  RcppParallel::setThreadOptions(numThreads = 1)
  d <- .two_dam(); mk <- d$mk
  phA <- .phased("A", mk, d$pvA, rep(1L, 8))
  phB <- .phased("B", mk, d$pvB, rep(1L, 8))
  pm <- structure(list(A = phA, B = phB), class = "HSMap.phased.multi")
  mb <- hmm_map_blocks(d$dat, pm, epsilon = 0.01, paternal_mode = "HWE", tol = 1e-7)
  mj <- hmm_map_joint(d$dat, phased = pm, dam = "all", epsilon = 0.01,
                      paternal_mode = "HWE", tol = 1e-7)
  expect_identical(mb$n_blocks, 1L)
  expect_equal(as.numeric(mb$blocks[[1]]$fit$fit$r), as.numeric(mj$fit$r), tolerance = 1e-7)
  # blockwise positions/lengths available and finite
  bm <- get_block_map(mb, "haldane")
  expect_true(is.finite(bm$total_linked_length))
  expect_identical(nrow(bm$positions), length(mk))
})


# Commit 1 (final round): joint NA-phase rejection + corrected block-boundary rule.
test_that("hmm_map_joint rejects unresolved (NA) phase before the C++ engine", {
  RcppParallel::setThreadOptions(numThreads = 1)
  d <- .two_dam(); mk <- d$mk
  phA <- .phased("A", mk, .na_at(d$pvA, 4L), rep(1L, 8))   # A has an NA interval
  phB <- .phased("B", mk, d$pvB,             rep(1L, 8))
  pm <- structure(list(A = phA, B = phB), class = "HSMap.phased.multi")
  # direct joint call is rejected, naming the affected dam + count
  expect_error(hmm_map_joint(d$dat, phased = pm, dam = "all", epsilon = 0.01,
                             paternal_mode = "HWE"),
               "unresolved.*'A'.*hmm_map_blocks")
  # the multi-dam path through hmm_map(method='joint') is rejected too
  expect_error(hmm_map(d$dat, phased = pm, dam = "all", epsilon = 0.01,
                       paternal_mode = "HWE"),
               "unresolved")
})

test_that("joint block rule: an interval unresolved by EVERY dam is always a boundary", {
  RcppParallel::setThreadOptions(numThreads = 1)
  d <- .two_dam(); mk <- d$mk
  # both dams unresolved AND non-informative at interval 4 (m4,m5 are singletons):
  # the old rule (any informative-and-unresolved) would MISS this; the corrected rule
  # (no dam resolves -> boundary) catches it.
  comp <- c(1, 1, 1, 2, 3, 4, 4, 4)
  phA <- .phased("A", mk, .na_at(d$pvA, 4L), comp)
  phB <- .phased("B", mk, .na_at(d$pvB, 4L), comp)
  pm <- structure(list(A = phA, B = phB), class = "HSMap.phased.multi")
  mb <- hmm_map_blocks(d$dat, pm, epsilon = 0.01, paternal_mode = "HWE")
  expect_identical(mb$n_blocks, 2L)
  expect_true(4L %in% mb$unresolved_boundaries)
  expect_false(is.null(mb$blocks[[1]]$fit))               # valid blocks BOTH sides
  expect_false(is.null(mb$blocks[[2]]$fit))
  bi <- Filter(function(z) z$interval == 4L, mb$boundary_info)[[1]]
  expect_identical(bi$reason, "no_dam_resolved")
})

test_that("joint block rule: boundary metadata records resolved / informative-unresolved / no-info dams", {
  RcppParallel::setThreadOptions(numThreads = 1)
  d <- .two_dam(); mk <- d$mk
  # A informative-but-unresolved at interval 4 (component split) -> boundary;
  # B fully resolved (contributes to both blocks).
  phA <- .phased("A", mk, .na_at(d$pvA, 4L), c(1, 1, 1, 1, 2, 2, 2, 2))
  phB <- .phased("B", mk, d$pvB,             rep(1L, 8))
  pm <- structure(list(A = phA, B = phB), class = "HSMap.phased.multi")
  mb <- hmm_map_blocks(d$dat, pm, epsilon = 0.01, paternal_mode = "HWE")
  bi <- Filter(function(z) z$interval == 4L, mb$boundary_info)[[1]]
  expect_identical(bi$reason, "informative_dam_unresolved")
  expect_true("A" %in% bi$dams_informative_unresolved)    # A forced the boundary
  expect_true("B" %in% bi$dams_resolved)                  # B resolves it
  # a single unsupported interval does not discard the chromosome: two valid blocks
  expect_identical(mb$n_blocks, 2L)
})
