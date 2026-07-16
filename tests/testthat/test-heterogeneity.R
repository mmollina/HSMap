# Tests for the per-dam recombination-map heterogeneity (eta) test.
# Shared helpers (make_dat, oracle_phased, oracle_multi) live in helper-sim.R.

test_that("heterogeneity test runs on a homogeneous joint map and finds eta ~ 1", {
  RcppParallel::setThreadOptions(numThreads = 1)
  set.seed(101); Tm <- 40
  r_true <- rep(0.04, Tm - 1)
  sim <- sim_multi_pop(T_markers = Tm, n_pops = 3, n_ind_per_pop = c(150, 150, 150),
                       marker_intersection = 1, r_vec = r_true,
                       phase_mode = "random", repulsion_rate = 0.3,
                       maternal_geno_mode = "all_het",
                       paternal_pA_base = 0.4, error_rate = 0.01, seed = 101)
  dat <- make_dat(sim); mk <- sim$truth$markers_union
  jm  <- hmm_map(dat, phased = oracle_multi(sim, mk), dam = "all",
                 epsilon = 0.01, paternal_mode = "gametic", maxit = 3000)

  het <- test_map_heterogeneity(dat, jm)
  expect_s3_class(het, "HSMap.hetero")
  expect_equal(length(het$eta), 3L)
  expect_true(all(het$eta > 0.6 & het$eta < 1.6))     # near 1 under homogeneity
  expect_true(het$LR >= 0)
  expect_gte(het$p_value, 0); expect_lte(het$p_value, 1)
  expect_equal(het$df, 2L)
})

test_that("heterogeneity test detects dams with different maps", {
  RcppParallel::setThreadOptions(numThreads = 1)
  Tm <- 40
  simA <- sim_multi_pop(T_markers = Tm, n_pops = 1, n_ind_per_pop = 200,
                        marker_intersection = 1, r_const = 0.02,
                        phase_mode = "random", repulsion_rate = 0.3,
                        maternal_geno_mode = "all_het",
                        paternal_pA_base = 0.4, error_rate = 0.01, seed = 11)
  simB <- sim_multi_pop(T_markers = Tm, n_pops = 1, n_ind_per_pop = 200,
                        marker_intersection = 1, r_const = 0.10,   # denser/longer map
                        phase_mode = "random", repulsion_rate = 0.3,
                        maternal_geno_mode = "all_het",
                        paternal_pA_base = 0.4, error_rate = 0.01, seed = 22)
  mk  <- simA$truth$markers_union
  dat <- structure(list(
    G_list = list(A = simA$G_list[[1]], B = simB$G_list[[1]]),
    M_list = list(A = simA$M_list[[1]], B = simB$M_list[[1]])), class = "HSMap.data")
  oph <- structure(list(
    A = oracle_phased(mk, 1L - simA$truth$v_true[[1]], "A"),
    B = oracle_phased(mk, 1L - simB$truth$v_true[[1]], "B")),
    class = "HSMap.phased.multi")

  jm  <- hmm_map(dat, phased = oph, dam = "all", epsilon = 0.01,
                 paternal_mode = "gametic", maxit = 3000)
  het <- test_map_heterogeneity(dat, jm)

  expect_gt(het$eta[["B"]], het$eta[["A"]])           # B has the longer map
  expect_gt(het$LR, stats::qchisq(0.99, df = 1))      # strong, non-flaky signal
  expect_lt(het$p_value, 0.01)
})


# Milestone 5: honest metadata for the conditional global-scale test.
test_that("heterogeneity test reports honest conditional global-scale metadata", {
  RcppParallel::setThreadOptions(numThreads = 1)
  set.seed(303); Tm <- 30
  sim <- sim_multi_pop(T_markers = Tm, n_pops = 3, n_ind_per_pop = c(120, 120, 120),
                       marker_intersection = 1, r_vec = rep(0.05, Tm - 1),
                       phase_mode = "random", repulsion_rate = 0.3,
                       maternal_geno_mode = "all_het",
                       paternal_pA_base = 0.4, error_rate = 0.01, seed = 303)
  dat <- make_dat(sim); mk <- sim$truth$markers_union
  jm  <- hmm_map(dat, phased = oracle_multi(sim, mk), dam = "all",
                 epsilon = 0.01, paternal_mode = "gametic", maxit = 3000)
  het <- test_map_heterogeneity(dat, jm)

  # test type + hypothesis structure
  expect_identical(het$test_type, "conditional global-scale LRT")
  expect_false(het$interval_specific)                       # NOT interval-specific
  expect_identical(het$n_params_null, 1L)                   # one common global scale
  expect_identical(het$n_params_alt, 3L)                    # one scale per dam (D=3)
  expect_identical(het$df, het$n_params_alt - het$n_params_null)  # df = D - 1
  # conditional-on and calibration are stated
  expect_identical(het$calibration, "asymptotic")
  expect_true(any(grepl("phase", het$conditional_on)))
  expect_true(any(grepl("paternal", het$conditional_on)))
  # per-dam boundary flag present and (here) not at a boundary
  expect_true("at_boundary" %in% names(het$per_dam))
  expect_false(het$any_boundary)
})

test_that("a boundary eta estimate is flagged", {
  RcppParallel::setThreadOptions(numThreads = 1)
  set.seed(404); Tm <- 25
  sim <- sim_multi_pop(T_markers = Tm, n_pops = 2, n_ind_per_pop = c(120, 120),
                       marker_intersection = 1, r_vec = rep(0.05, Tm - 1),
                       phase_mode = "random", repulsion_rate = 0.3,
                       maternal_geno_mode = "all_het",
                       paternal_pA_base = 0.4, error_rate = 0.01, seed = 404)
  dat <- make_dat(sim); mk <- sim$truth$markers_union
  jm  <- hmm_map(dat, phased = oracle_multi(sim, mk), dam = "all",
                 epsilon = 0.01, paternal_mode = "gametic", maxit = 3000)
  # a narrow eta_range forces the optimum to the boundary -> flagged
  het <- test_map_heterogeneity(dat, jm, eta_range = c(1.5, 3))
  expect_true(het$any_boundary)
  expect_true(any(het$per_dam$at_boundary))
})
