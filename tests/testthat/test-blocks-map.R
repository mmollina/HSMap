# Milestone 3: blockwise multipoint fitting + safe map reporting.
#   * hmm_map() rejects unresolved (NA) phase; hmm_map_blocks() splits at gaps,
#   * fully-resolved chromosome: ordinary == blockwise,
#   * r = 0.5 -> gap (NA distance), never ~1577 cM; cumulative resets between blocks.
# Shared helpers (one_dam etc.) live in helper-sim.R.

test_that("hmm_map() rejects unresolved (NA) phase with a clear error", {
  RcppParallel::setThreadOptions(numThreads = 1)
  d <- one_dam(1L, Tm = 8L, n = 120L)
  ph <- d$oph; ph$phase_vec[3] <- NA_integer_
  expect_error(
    hmm_map(d$dat, phased = ph, dam = 1, epsilon = 0.01, paternal_mode = "HWE"),
    "unresolved"
  )
})

test_that("blockwise fitting splits at unresolved intervals", {
  RcppParallel::setThreadOptions(numThreads = 1)
  d <- one_dam(2L, Tm = 8L, n = 120L)
  ph <- d$oph; ph$phase_vec[4] <- NA_integer_          # one unresolved interval
  mb <- hmm_map_blocks(d$dat, ph, epsilon = 0.01, paternal_mode = "HWE", tol = 1e-6)
  expect_s3_class(mb, "HSMap.map.blocks")
  expect_identical(mb$n_blocks, 2L)
  expect_true(4L %in% mb$unresolved_boundaries)
  expect_identical(mb$blocks[[1]]$markers, d$mk[1:4])
  expect_identical(mb$blocks[[2]]$markers, d$mk[5:8])
  # the between-block interval is reported as such
  expect_identical(mb$interval_table$status[4], "unresolved_phase")
})

test_that("a fully resolved chromosome gives the same result via ordinary and blockwise", {
  RcppParallel::setThreadOptions(numThreads = 1)
  d <- one_dam(3L, Tm = 10L, n = 150L)
  m  <- hmm_map(d$dat, phased = d$oph, dam = 1, epsilon = 0.01, paternal_mode = "HWE", tol = 1e-7)
  mb <- hmm_map_blocks(d$dat, d$oph, epsilon = 0.01, paternal_mode = "HWE", tol = 1e-7)
  expect_identical(mb$n_blocks, 1L)
  expect_equal(as.numeric(mb$blocks[[1]]$fit$fit$r), as.numeric(m$fit$r), tolerance = 1e-8)
})

test_that("r = 0.5 produces a map gap, not a huge centimorgan distance", {
  mk <- paste0("m", 1:4)
  map <- structure(list(order = mk, fit = list(r = c(0.05, 0.50, 0.08))),
                   class = "HSMap.map")
  pos <- get_map(map, "haldane")
  st  <- attr(pos, "status")
  dd  <- attr(pos, "dist_cM")
  expect_identical(unname(st[2]), "no_linkage_boundary")
  expect_true(is.na(dd[2]))                            # gap: NA distance
  expect_lt(attr(pos, "total_linked_length"), 20)     # NOT ~1577 cM
  expect_identical(attr(pos, "n_gaps"), 1L)
})

test_that("cumulative positions reset between blocks (after a gap)", {
  mk <- paste0("m", 1:4)
  map <- structure(list(order = mk, fit = list(r = c(0.05, 0.50, 0.05))),
                   class = "HSMap.map")
  pos <- get_map(map, "morgan")                        # 100*r within block
  # block 1: m1=0, m2=5 ; gap ; block 2 resets: m3=0, m4=5
  expect_equal(as.numeric(pos), c(0, 5, 0, 5), tolerance = 1e-8)
  expect_identical(unname(attr(pos, "block")), c(1L, 1L, 2L, 2L))
  expect_equal(attr(pos, "within_block_length"), c(5, 5), tolerance = 1e-8)
})

test_that("get_map handles unresolved-phase gaps without crashing", {
  mk <- paste0("m", 1:4)
  map <- structure(list(order = mk, fit = list(r = c(0.05, 0.10, 0.05)),
                        resolved_interval = c(TRUE, FALSE, TRUE)),
                   class = "HSMap.map")
  expect_silent(pos <- get_map(map, "kosambi"))
  expect_identical(unname(attr(pos, "status")[2]), "unresolved_phase")
  expect_true(is.na(attr(pos, "dist_cM")[2]))
  expect_identical(unname(attr(pos, "block")), c(1L, 1L, 2L, 2L))
})

test_that("block metadata and interval status remain length-aligned", {
  RcppParallel::setThreadOptions(numThreads = 1)
  d <- one_dam(4L, Tm = 9L, n = 120L)
  ph <- d$oph; ph$phase_vec[c(3L, 6L)] <- NA_integer_  # two gaps -> three blocks
  mb <- hmm_map_blocks(d$dat, ph, epsilon = 0.01, paternal_mode = "HWE", tol = 1e-6)
  expect_identical(mb$n_blocks, 3L)
  expect_identical(length(mb$block_id), length(d$mk))
  expect_identical(nrow(mb$interval_table), length(d$mk) - 1L)
  expect_identical(sum(mb$interval_table$status == "unresolved_phase"), 2L)
  # every within-block interval has a finite fitted r; boundary intervals are NA
  expect_true(all(is.na(mb$interval_table$r[mb$interval_table$status == "unresolved_phase"])))
  expect_true(all(is.finite(mb$interval_table$r[mb$interval_table$status == "linked"])))
})
