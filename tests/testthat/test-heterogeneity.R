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


# Commit 3: gap-safe heterogeneity test.
test_that("heterogeneity test accepts a valid block but rejects r=0.5 / NA / non-converged", {
  RcppParallel::setThreadOptions(numThreads = 1)
  set.seed(505); Tm <- 25
  sim <- sim_multi_pop(T_markers = Tm, n_pops = 3, n_ind_per_pop = c(120, 120, 120),
                       marker_intersection = 1, r_vec = rep(0.05, Tm - 1),
                       phase_mode = "random", repulsion_rate = 0.3,
                       maternal_geno_mode = "all_het",
                       paternal_pA_base = 0.4, error_rate = 0.01, seed = 505)
  dat <- make_dat(sim); mk <- sim$truth$markers_union
  jm  <- hmm_map(dat, phased = oracle_multi(sim, mk), dam = "all",
                 epsilon = 0.01, paternal_mode = "gametic", maxit = 3000)

  # a valid, linked, resolved, converged joint map is accepted
  expect_s3_class(test_map_heterogeneity(dat, jm), "HSMap.hetero")

  # r at/above the no-linkage threshold is rejected (r = 0.5 is NOT silently clamped)
  jm_gap <- jm; jm_gap$fit$r[3] <- 0.5
  expect_error(test_map_heterogeneity(dat, jm_gap), "no-linkage")

  # non-finite r is rejected
  jm_na <- jm; jm_na$fit$r[2] <- NA_real_
  expect_error(test_map_heterogeneity(dat, jm_na), "non-finite")

  # a non-converged fit is rejected
  jm_nc <- jm; jm_nc$fit$converged <- FALSE
  expect_error(test_map_heterogeneity(dat, jm_nc), "did not converge")

  # unresolved-phase metadata is rejected
  jm_up <- jm; jm_up$resolved_interval <- c(FALSE, rep(TRUE, Tm - 2))
  expect_error(test_map_heterogeneity(dat, jm_up), "unresolved")
})


# Commit 3 (final round): uncapped scaling; objective-decreased rejection.
test_that("heterogeneity scaled map uses the uncapped biological formula (not capped at gap_r)", {
  RcppParallel::setThreadOptions(numThreads = 1)
  Tm <- 20
  # a shared high-recombination interval (large Haldane m) + dam B a denser map,
  # so the scaled r can exceed gap_r for large eta.
  rv <- c(0.42, runif(Tm - 2, 0.02, 0.08))
  simA <- sim_multi_pop(T_markers = Tm, n_pops = 1, n_ind_per_pop = 200, marker_intersection = 1,
                        r_vec = rv, phase_mode = "random", repulsion_rate = 0.3,
                        maternal_geno_mode = "all_het", paternal_pA_base = 0.4, error_rate = 0.01, seed = 71)
  simB <- sim_multi_pop(T_markers = Tm, n_pops = 1, n_ind_per_pop = 200, marker_intersection = 1,
                        r_vec = pmin(rv * 2, 0.45), phase_mode = "random", repulsion_rate = 0.3,
                        maternal_geno_mode = "all_het", paternal_pA_base = 0.4, error_rate = 0.01, seed = 72)
  mk <- simA$truth$markers_union
  dat <- structure(list(G_list = list(A = simA$G_list[[1]], B = simB$G_list[[1]]),
                        M_list = list(A = simA$M_list[[1]], B = simB$M_list[[1]])), class = "HSMap.data")
  oph <- structure(list(A = oracle_phased(mk, 1L - simA$truth$v_true[[1]], "A"),
                        B = oracle_phased(mk, 1L - simB$truth$v_true[[1]], "B")),
                   class = "HSMap.phased.multi")
  jm <- hmm_map(dat, phased = oph, dam = "all", epsilon = 0.01, paternal_mode = "gametic", maxit = 3000)
  het <- test_map_heterogeneity(dat, jm, eta_range = c(0.1, 10))

  r <- pmin(as.numeric(jm$fit$r), 0.5 - 1e-12); m <- -0.5 * log(1 - 2 * r)
  r_scaled <- function(eta) pmin(0.5 * (1 - exp(-2 * eta * m)), 0.5 - 1e-12)
  eps <- jm$fit$epsilon %||% 0.01
  align <- function(dn) {
    G <- matrix(NA_integer_, nrow(dat$G_list[[dn]]), Tm, dimnames = list(NULL, mk))
    cg <- intersect(mk, colnames(dat$G_list[[dn]])); G[, cg] <- dat$G_list[[dn]][, cg]; storage.mode(G) <- "integer"
    M <- rep(NA_integer_, Tm); names(M) <- mk; cm <- intersect(mk, names(dat$M_list[[dn]])); M[cm] <- as.integer(dat$M_list[[dn]][cm])
    list(G = G, M = M, ph = as.integer(jm$phase_list[[dn]]), emis = as.matrix(jm$fit$pi_list[[dn]]))
  }
  for (i in seq_len(nrow(het$per_dam))) {
    dn <- het$per_dam$dam[i]; eta <- het$eta[[dn]]; a <- align(dn)
    ll <- HSMap:::loglik_hs_cpp(a$G, a$M, a$ph, r_scaled(eta), a$emis, eps)
    expect_equal(ll, het$per_dam$ll_alt[i], tolerance = 1e-5)   # function uses r_scaled(eta), uncapped
  }
  # the biological scaling is NOT capped at gap_r: at a large eta it reaches above gap_r
  expect_gt(max(r_scaled(9)), 0.499)
  expect_lt(max(r_scaled(9)), 0.5)
})

test_that("heterogeneity rejects a fit whose objective decreased", {
  RcppParallel::setThreadOptions(numThreads = 1)
  set.seed(606); Tm <- 20
  sim <- sim_multi_pop(T_markers = Tm, n_pops = 3, n_ind_per_pop = c(120, 120, 120),
                       marker_intersection = 1, r_vec = rep(0.05, Tm - 1),
                       phase_mode = "random", repulsion_rate = 0.3, maternal_geno_mode = "all_het",
                       paternal_pA_base = 0.4, error_rate = 0.01, seed = 606)
  dat <- make_dat(sim); mk <- sim$truth$markers_union
  jm  <- hmm_map(dat, phased = oracle_multi(sim, mk), dam = "all",
                 epsilon = 0.01, paternal_mode = "gametic", maxit = 3000)
  expect_s3_class(test_map_heterogeneity(dat, jm), "HSMap.hetero")   # normal fit accepted
  jm_bad <- jm; jm_bad$fit$objective_decreased <- TRUE
  expect_error(test_map_heterogeneity(dat, jm_bad), "objective_decreased")
})
