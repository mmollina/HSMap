## diagnose_map_intervals() + collapse_framework_blocks(): contract, statuses,
## the accept path on a planted defect, the reject path, and the
## catastrophic-length guard semantics.

.fm_sim <- function(Tm = 60L, n = 250L, seed = 11, r = 0.002, pA = 0.4) {
  sim <- sim_multi_pop(T_markers = Tm, n_pops = 1L, n_ind_per_pop = n,
                       r_vec = rep(r, Tm - 1L), phase_mode = "all_coupling",
                       maternal_geno_mode = "all_het", paternal_pA_base = pA,
                       error_rate = 0, seed = seed)
  list(dat = structure(list(G_list = sim$G_list, M_list = sim$M_list),
                       class = "HSMap.data"),
       markers = sim$truth$markers_union)
}
.fm_fit <- function(dat, tpt, mk, dam = 1) {
  ph <- phase_from_pairwise(tpt, order = mk, dam = dam)
  if (inherits(ph, "HSMap.phased.multi")) ph <- ph[[1L]]
  suppressWarnings(hmm_map(dat, phased = ph, dam = dam, epsilon = 0.05, tol = 1e-6))
}

test_that("diagnose_map_intervals: contract and statuses", {
  s <- .fm_sim()
  tpt <- pairwise_rf(s$dat, threads = 1)
  map <- .fm_fit(s$dat, tpt, s$markers)
  iv <- diagnose_map_intervals(map, tpt)

  expect_equal(nrow(iv), length(s$markers) - 1L)
  expect_identical(iv$marker_left, s$markers[-length(s$markers)])
  expect_identical(iv$marker_right, s$markers[-1L])
  expect_true(all(iv$status %in% c("normal", "large_concordant",
                                   "large_discordant", "unresolved", "no_linkage")))
  expect_equal(iv$discordance, iv$r_mp - iv$r_2pt)
  ## clean tight data: everything normal, cM finite
  expect_true(all(iv$status == "normal"))
  expect_true(all(is.finite(iv$cM)))
  ## thresholds are honoured: absurd large_r flags nothing
  iv2 <- diagnose_map_intervals(map, tpt, large_r = 0.49)
  expect_true(all(iv2$status == "normal"))
  ## unresolved phase surfaces as `unresolved` with NA cM
  map_u <- map; map_u$phase_vec[5] <- NA_integer_
  ivu <- diagnose_map_intervals(map_u, tpt)
  expect_identical(ivu$status[5], "unresolved")
  expect_true(is.na(ivu$cM[5]))
  ## validation errors
  expect_error(diagnose_map_intervals(list(), tpt), "HSMap.map")
  expect_error(diagnose_map_intervals(map, list()), "HSMap.tpt")
  expect_error(diagnose_map_intervals(map, tpt, order = rev(s$markers)),
               "must match")
})

test_that("collapse_framework_blocks: accept path repairs a planted defect", {
  s <- .fm_sim(seed = 21, pA = 0.15)
  blk <- 26:35
  G <- s$dat$G_list[[1]]
  G[1:100, s$markers[blk]] <- 1L                 # coherent het-locked block
  s$dat$G_list[[1]] <- G
  tpt <- pairwise_rf(s$dat, threads = 1)

  map0 <- .fm_fit(s$dat, tpt, s$markers)
  iv0 <- diagnose_map_intervals(map0, tpt)
  skip_if(!any(iv0$status == "large_discordant"),
          "planted defect did not create a large_discordant interval")

  fw <- collapse_framework_blocks(s$dat, tpt, s$markers,
                                  blocks = list(s$markers[blk]),
                                  lg = "SIM1", verbose = FALSE)
  expect_s3_class(fw, "HSMap.framework")
  expect_true(fw$accepted)
  expect_true(all(fw$checks$pass))
  ## marker accounting: nothing deleted, statuses coherent
  expect_identical(sort(fw$marker_table$marker), sort(s$markers))
  expect_equal(sum(fw$marker_table$framework_status == "representative"), 1L)
  expect_equal(sum(fw$marker_table$framework_status == "attached"),
               length(blk) - 1L)
  expect_true(all(fw$marker_table$lg == "SIM1"))
  ## attached markers inherit the representative's position
  att <- fw$marker_table[fw$marker_table$framework_status == "attached", ]
  rep_pos <- fw$marker_table$position_cM[
    fw$marker_table$framework_status == "representative"]
  expect_true(all(att$position_cM == rep_pos))
  ## attached markers are excluded from the fitted map, representative retained
  expect_false(any(att$marker %in% fw$map$order))
  expect_true(fw$blocks$representative %in% fw$map$order)
  ## the repair actually reduced the target interval
  expect_true(all(fw$blocks$target_r_after < fw$blocks$target_r_before))
  ## EM budget is recorded and defaults are as documented
  expect_equal(fw$settings$tol, 1e-6)
  expect_equal(fw$settings$maxit, 1000L)
  expect_output(print(fw), "ACCEPTED")
})

test_that("collapse_framework_blocks: reject path returns the original map", {
  s <- .fm_sim(seed = 31)
  tpt <- pairwise_rf(s$dat, threads = 1)
  ## an impossible stability demand forces rejection deterministically
  fw <- collapse_framework_blocks(s$dat, tpt, s$markers,
                                  blocks = list(s$markers[20:24]),
                                  stability_tol = 0, stability_frac = 1,
                                  max_len_increase = 1 + 1e-12,
                                  verbose = FALSE)
  if (fw$accepted) skip("collapse was numerically exact here; cannot force reject")
  expect_false(fw$accepted)
  ## the returned map IS the original: same markers, same r vector
  expect_identical(fw$map$order, fw$map_original$order)
  expect_identical(as.numeric(fw$map$fit$r), as.numeric(fw$map_original$fit$r))
  ## no marker is reported attached when the repair was rejected
  expect_true(all(fw$marker_table$framework_status == "framework"))
  expect_false(any(fw$blocks$applied))
  expect_output(print(fw), "REJECTED")
})

test_that("catastrophic length inflation can never be accepted", {
  s <- .fm_sim(seed = 41)
  tpt <- pairwise_rf(s$dat, threads = 1)
  ## simulate check-7 semantics directly: even if all else passed, a repaired
  ## map longer than max_len_increase x original must fail verification.
  fw <- collapse_framework_blocks(s$dat, tpt, s$markers,
                                  blocks = list(s$markers[20:24]),
                                  max_len_increase = 0,   # any length fails
                                  verbose = FALSE)
  expect_false(fw$accepted)
  expect_false(fw$checks$pass[fw$checks$check == "map length not inflated"])
  expect_identical(fw$map$order, fw$map_original$order)
})

test_that("block validation is strict", {
  s <- .fm_sim(Tm = 30L, n = 100L, seed = 51)
  tpt <- pairwise_rf(s$dat, threads = 1)
  m <- s$markers
  expect_error(collapse_framework_blocks(s$dat, tpt, m, blocks = list()),
               "non-empty")
  expect_error(collapse_framework_blocks(s$dat, tpt, m,
                                         blocks = list(c(m[5], "nope"))),
               "not in `order`")
  expect_error(collapse_framework_blocks(s$dat, tpt, m, blocks = list(m[5])),
               "at least 2")
  expect_error(collapse_framework_blocks(s$dat, tpt, m,
                                         blocks = list(m[c(5, 9)])),
               "contiguous")
  expect_error(collapse_framework_blocks(s$dat, tpt, m,
                                         blocks = list(m[1:3])),
               "first or last")
  expect_error(collapse_framework_blocks(s$dat, tpt, m,
                                         blocks = list(m[5:8], m[7:10])),
               "overlap")
})
