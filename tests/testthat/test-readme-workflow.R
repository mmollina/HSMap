# The README's open-pollinated example workflow, run end to end on the packaged
# simulated example data. If this test breaks, the README example is out of date.

test_that("the full open-pollinated README workflow runs on the example data", {
  RcppParallel::setThreadOptions(numThreads = 1L)

  ped  <- system.file("extdata", "example_pedigree.csv",  package = "HSMap")
  geno <- system.file("extdata", "example_genotypes.csv", package = "HSMap")
  expect_true(nzchar(ped) && nzchar(geno))

  # 1. read pedigree + genotype files
  dat <- read_HSMap_data(ped, geno)
  expect_s3_class(dat, "HSMap.data")
  expect_gte(length(dat$G_list), 1L)

  # 2. pairwise (two-point) analysis
  tpt <- pairwise_rf(dat, threads = 1)
  expect_s3_class(tpt, "HSMap.tpt")

  # 3. filter the two-point table
  tptf <- tpt_filter(tpt, diagnostic.plot = FALSE)
  expect_s3_class(tptf, "HSMap.tpt")

  # 4. group markers into linkage groups (the example is one simulated chromosome)
  grp <- group_markers(tptf, k = 1, inter = FALSE)
  expect_s3_class(grp, "hsmap_group")

  # 5. MDS ordering within each linkage group
  ord <- mds_order(grp, tptf, plot_each = FALSE)
  expect_s3_class(ord, "hsmap_mds")
  o1 <- ord[[1]]
  expect_true(is.character(o1) && length(o1) >= 3L)

  # 6. dam-specific phase estimation
  ph <- phase_from_pairwise(tptf, order = o1, dam = "all")
  expect_s3_class(ph, "HSMap.phased.multi")

  # 7. blockwise multipoint HMM mapping (safe with unresolved phase)
  blocks <- hmm_map_blocks(dat, ph)
  expect_s3_class(blocks, "HSMap.map.blocks")
  expect_gte(blocks$n_blocks, 1L)

  # 8. safe map positions and lengths (gaps are NA, never huge cM)
  bm <- get_block_map(blocks, "haldane")
  expect_true(is.finite(bm$total_linked_length))
  expect_true(all(bm$interval_table$status %in%
    c("linked", "no_linkage_boundary", "unresolved_phase",
      "between_blocks", "insufficient_information")))
  # no gap becomes a huge finite distance
  gaps <- bm$interval_table$status != "linked"
  expect_true(all(is.na(bm$interval_table$dist_cM[gaps])))

  # 9. plotting phase blocks and map segments
  skip_if_not_installed("ggplot2")
  p <- plot_block_map(blocks, map.function = "haldane")
  expect_s3_class(p, "ggplot")
})
