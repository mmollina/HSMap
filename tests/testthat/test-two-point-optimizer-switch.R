## Optimizer selection tests (production behavior).
##
## Default is "auto": a SINGLE family uses the exact concave phase-decomposed
## Newton solve; several families fall back to the grid + refinement search.
## grid_refine_max() stays in the package as the reference implementation,
## reachable through the internal option HSMap.two_point_optimizer = "grid"
## (regression/comparison only -- not public API).

.mk_dat1 <- function(Tm = 60L, n = 200L, seed = 71, r_vec = NULL) {
  if (is.null(r_vec)) r_vec <- rep(0.06, Tm - 1L)
  sim <- sim_multi_pop(T_markers = Tm, n_pops = 1L, n_ind_per_pop = n,
                       r_vec = r_vec, phase_mode = "random",
                       repulsion_rate = 0.3, maternal_geno_mode = "all_het",
                       paternal_pA_base = 0.4, error_rate = 0.01, seed = seed)
  structure(list(G_list = sim$G_list, M_list = sim$M_list), class = "HSMap.data")
}

test_that("single family defaults to the newton path; grid stays available as reference", {
  dat <- .mk_dat1()
  old <- options(HSMap.two_point_optimizer = NULL); on.exit(options(old), add = TRUE)

  tn <- pairwise_rf(dat, threads = 1)                     # production default
  expect_identical(tn$fit$optimizer, "newton-concave")

  options(HSMap.two_point_optimizer = "grid")             # internal reference
  tg <- pairwise_rf(dat, threads = 1)
  expect_identical(tg$fit$optimizer, "grid+local-refine")
  options(HSMap.two_point_optimizer = NULL)

  ut <- upper.tri(tg$fit$r)
  ## same estimates to optimizer tolerance; identical calls
  expect_lt(max(abs(tg$fit$r - tn$fit$r)[ut], na.rm = TRUE), 1e-4)
  expect_lt(max(abs(tg$fit$lod_r - tn$fit$lod_r)[ut], na.rm = TRUE), 1e-6)
  expect_lt(max(abs(tg$fit$lod_ph - tn$fit$lod_ph)[ut], na.rm = TRUE), 1e-3)
  pg <- tg$fit$mom_phase_list[[1]][ut]
  pn <- tn$fit$mom_phase_list[[1]][ut]
  expect_identical(is.na(pg), is.na(pn))
  expect_true(all(pg == pn, na.rm = TRUE))
  expect_identical(tg$fit$no_linkage[ut], tn$fit$no_linkage[ut])
})

test_that("multi-family data automatically falls back to the grid", {
  sim <- sim_multi_pop(T_markers = 20L, n_pops = 2L, n_ind_per_pop = c(40L, 40L),
                       r_vec = rep(0.05, 19L), phase_mode = "random",
                       repulsion_rate = 0.3, maternal_geno_mode = "all_het",
                       paternal_pA_base = 0.4, error_rate = 0.01, seed = 72)
  dat3 <- structure(list(G_list = sim$G_list, M_list = sim$M_list),
                    class = "HSMap.data")
  old <- options(HSMap.two_point_optimizer = NULL); on.exit(options(old), add = TRUE)

  tpt <- pairwise_rf(dat3, threads = 1)                   # no option, no error
  expect_s3_class(tpt, "HSMap.tpt")
  expect_identical(tpt$fit$optimizer, "grid+local-refine")

  ## explicit "newton" on multi-family data must still refuse
  options(HSMap.two_point_optimizer = "newton")
  expect_error(pairwise_rf(dat3, threads = 1), "single family")
})

test_that("invalid optimizer option is rejected", {
  dat <- .mk_dat1(Tm = 10L, n = 40L, seed = 73)
  old <- options(HSMap.two_point_optimizer = "banana"); on.exit(options(old), add = TRUE)
  expect_error(pairwise_rf(dat, threads = 1), "grid")
})

test_that("r_hat exactly 0.5 yields phase NA and phase LOD 0 (production path)", {
  ## Deterministic construction of an exact-null pair. With balanced counts
  ## C[0][0] = C[2][2] = C[0][2] = C[2][0] = 25 and balanced AA/aa marginals,
  ## q-hat = (50 + 10)/(100 + 20) = 0.5 exactly, and every probability at
  ## q = 0.5 is an exact dyadic rational (0.25, 0.125, 0.0625). The derivative
  ## terms at r = 0.5 are then the exact integers -50, -50, +50, +50, so
  ## l'(0.5) == 0 in floating point and the Newton boundary branch returns
  ## BIT-EXACT r = 0.5. The explicit null rule must then report phase NA and
  ## phase LOD 0 -- never a summation-order ulp call.
  G <- rbind(
    matrix(rep(c(0L, 0L, 1L), each = 25), 25, 3, byrow = FALSE),
    matrix(rep(c(2L, 2L, 1L), each = 25), 25, 3, byrow = FALSE),
    matrix(rep(c(0L, 2L, 1L), each = 25), 25, 3, byrow = FALSE),
    matrix(rep(c(2L, 0L, 1L), each = 25), 25, 3, byrow = FALSE)
  )
  colnames(G) <- paste0("m", 1:3)
  dat <- structure(list(
    G_list = list(P1 = G),
    M_list = list(P1 = stats::setNames(c(1L, 1L, 1L), colnames(G)))
  ), class = "HSMap.data")
  old <- options(HSMap.two_point_optimizer = NULL); on.exit(options(old), add = TRUE)

  tpt <- pairwise_rf(dat, threads = 1)                    # production default
  expect_identical(tpt$fit$optimizer, "newton-concave")
  expect_identical(tpt$fit$r["m1", "m2"], 0.5)            # bit-exact null
  expect_identical(tpt$fit$no_linkage["m1", "m2"], 1L)
  expect_true(is.na(tpt$fit$mom_phase_list[[1]]["m1", "m2"]))
  expect_identical(tpt$fit$lod_ph["m1", "m2"], 0)
  expect_identical(tpt$fit$lod_ph_list[[1]]["m1", "m2"], 0)

  ## the grid reference obeys the same rule wherever it lands on bit-exact 0.5
  options(HSMap.two_point_optimizer = "grid")
  tg <- pairwise_rf(dat, threads = 1)
  at_half <- !is.na(tg$fit$r) & tg$fit$r == 0.5 & upper.tri(tg$fit$r)
  expect_true(all(is.na(tg$fit$mom_phase_list[[1]][at_half])))
  expect_true(all(tg$fit$lod_ph[at_half] == 0))
  options(HSMap.two_point_optimizer = NULL)
})
