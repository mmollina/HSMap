# Commit 2 (final round): consistent gap_r override, map segments, gap-safe plotting.

.mk_seg_blocks <- function(r_mid, gap_r = 0.499) {
  mk <- paste0("m", 1:4)
  f  <- structure(list(order = mk, fit = list(r = c(0.05, r_mid, 0.05))), class = "HSMap.map")
  blocks <- list(list(block = 1, markers = mk, contributing_dams = "A", fit = f))
  block_id <- c(1L, 1L, 1L, 1L); pv <- c(1L, 1L, 1L)   # one fully resolved phase block
  itab <- HSMap:::.block_interval_table(mk, pv, block_id, blocks, gap_r = gap_r,
                                        boundary_status = "unresolved_phase")
  structure(list(blocks = blocks, block_id = stats::setNames(block_id, mk),
                 unresolved_boundaries = integer(0), n_blocks = 1L, order = mk,
                 interval_table = itab, gap_r = gap_r, dams = "A"),
            class = "HSMap.map.blocks")
}

test_that("overriding gap_r recomputes statuses, distances, positions, and gaps consistently", {
  mb <- .mk_seg_blocks(r_mid = 0.48, gap_r = 0.499)      # stored: 0.48 < 0.499 -> linked
  bm_def <- get_block_map(mb, "haldane")                 # stored threshold
  bm_lo  <- get_block_map(mb, "haldane", gap_r = 0.45)   # override: 0.48 >= 0.45 -> gap
  # default: all linked, one segment, no gaps
  expect_identical(bm_def$interval_table$status[2], "linked")
  expect_true(is.finite(bm_def$interval_table$dist_cM[2]))
  expect_identical(bm_def$n_gaps, 0L); expect_identical(bm_def$n_segments, 1L)
  # override recomputes: interval 2 becomes a gap -> NA distance, 2 segments, reset pos
  expect_identical(bm_lo$interval_table$status[2], "no_linkage_boundary")
  expect_true(is.na(bm_lo$interval_table$dist_cM[2]))
  expect_identical(bm_lo$n_gaps, 1L); expect_identical(bm_lo$n_segments, 2L)
  expect_equal(bm_lo$positions$pos[3], 0)                # segment reset after the gap
  expect_identical(bm_lo$gap_r, 0.45)
  expect_lt(bm_lo$total_linked_length, bm_def$total_linked_length)  # the gap drops a distance
})

test_that("a no-linkage interval inside one phase block creates two segments (one block)", {
  mb <- .mk_seg_blocks(r_mid = 0.5)                       # internal no-linkage at default gap_r
  bm <- get_block_map(mb, "haldane")
  expect_identical(bm$n_segments, 2L)                     # two segments...
  expect_true(all(bm$positions$phase_block == 1L))        # ...within ONE phase block
  expect_true(is.na(bm$interval_table$dist_cM[2]))
  expect_lt(bm$total_linked_length, 30)                  # no large finite distance
  # block & segment metadata stay aligned
  expect_identical(nrow(bm$positions), length(mb$order))
  expect_identical(nrow(bm$interval_table), length(mb$order) - 1L)
  expect_identical(length(bm$segment_lengths), bm$n_segments)
  skip_if_not_installed("ggplot2")
  expect_s3_class(plot_block_map(mb), "ggplot")           # facets by segment
})

test_that("get_block_map validates gap_r even when no block was fitted", {
  mk <- c("m1", "m2")
  blocks <- list(list(block = 1, markers = "m1", contributing_dams = character(0), fit = NULL),
                 list(block = 2, markers = "m2", contributing_dams = character(0), fit = NULL))
  itab <- data.frame(from = "m1", to = "m2", block = NA_integer_,
                     status = "unresolved_phase", r = NA_real_, stringsAsFactors = FALSE)
  mb <- structure(list(blocks = blocks, block_id = stats::setNames(c(1L, 2L), mk),
                       unresolved_boundaries = 1L, n_blocks = 2L, order = mk,
                       interval_table = itab, gap_r = 0.499, dams = "A"),
                  class = "HSMap.map.blocks")
  expect_error(get_block_map(mb, "haldane", gap_r = 0.9), "gap_r")
  expect_error(get_block_map(mb, "haldane", gap_r = -1),  "gap_r")
  expect_silent(get_block_map(mb, "haldane"))             # no fitted block, valid gap_r
})

test_that("plot_map_list is gap-safe: it stops (not silently continuous) on a gapped map", {
  mk <- paste0("m", 1:4)
  map_gap <- structure(list(order = mk, fit = list(r = c(0.05, 0.5, 0.05))), class = "HSMap.map")
  expect_error(plot_map_list(map_gap), "gap")
})
