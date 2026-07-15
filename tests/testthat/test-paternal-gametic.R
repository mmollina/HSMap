# Tests for the identifiable paternal gametic allele-frequency model (M1).
#
# The paternal contribution enters the likelihood only through the sire gametic
# allele frequency q_k = P(paternal gamete transmits A) = pi_AA + 0.5*pi_Aa.
# `gametic` (default) and `HWE` are the same identifiable estimator; the free
# 3-genotype `per_marker` model and the 10-class `two_locus` model are
# non-identifiable and disabled at the public API.
#
# Shared helpers (make_dat, oracle_phased, oracle_multi) live in helper-sim.R.

test_that("emissions/likelihood depend on pi only through q = piAA + 0.5*piAa", {
  RcppParallel::setThreadOptions(numThreads = 1)
  set.seed(1); Tm <- 12
  sim <- sim_multi_pop(T_markers = Tm, n_pops = 1, n_ind_per_pop = 300,
                       marker_intersection = 1, r_const = 0.08,
                       phase_mode = "all_coupling", maternal_geno_mode = "all_het",
                       paternal_pA_base = 0.6, error_rate = 0, seed = 1)
  G <- sim$G_list[[1]]; storage.mode(G) <- "integer"
  M <- as.integer(sim$M_list[[1]])
  phase <- rep(1L, Tm - 1); r <- rep(0.08, Tm - 1)
  mk_pi <- function(a, b, c) matrix(c(a, b, c), nrow = 3, ncol = Tm)

  llA <- HSMap:::loglik_hs_cpp(G, M, phase, r, mk_pi(0.25, 0.50, 0.25), 1e-3)  # q = 0.50
  llB <- HSMap:::loglik_hs_cpp(G, M, phase, r, mk_pi(0.50, 0.00, 0.50), 1e-3)  # q = 0.50
  llC <- HSMap:::loglik_hs_cpp(G, M, phase, r, mk_pi(0.36, 0.48, 0.16), 1e-3)  # q = 0.60

  expect_equal(llA, llB, tolerance = 1e-9)   # same q  -> identical likelihood
  expect_gt(abs(llA - llC), 1e-6)            # different q -> different likelihood
})

test_that("gametic and HWE are numerically identical (n_threads = 1)", {
  RcppParallel::setThreadOptions(numThreads = 1)
  set.seed(7); Tm <- 40
  r_true <- runif(Tm - 1, 0.01, 0.15)
  sim <- sim_multi_pop(T_markers = Tm, n_pops = 1, n_ind_per_pop = 200,
                       marker_intersection = 1, r_vec = r_true,
                       phase_mode = "random", repulsion_rate = 0.3,
                       maternal_geno_mode = "all_het",
                       paternal_pA_base = 0.4, error_rate = 0.01, seed = 7)
  dat <- make_dat(sim); mk <- sim$truth$markers_union
  oph <- oracle_phased(mk, 1L - sim$truth$v_true[[1]], "P1")

  mg <- hmm_map(dat, phased = oph, dam = 1, epsilon = 0.01,
                paternal_mode = "gametic", tol = 1e-6, maxit = 500)
  mh <- hmm_map(dat, phased = oph, dam = 1, epsilon = 0.01,
                paternal_mode = "HWE", tol = 1e-6, maxit = 500)

  expect_identical(mg$fit$r, mh$fit$r)
  expect_identical(mg$fit$logLik, mh$fit$logLik)
})

test_that("changing only the nonidentifiable pi decomposition cannot change r", {
  RcppParallel::setThreadOptions(numThreads = 1)
  set.seed(7); Tm <- 30
  r_true <- runif(Tm - 1, 0.01, 0.15)
  sim <- sim_multi_pop(T_markers = Tm, n_pops = 1, n_ind_per_pop = 200,
                       marker_intersection = 1, r_vec = r_true,
                       phase_mode = "random", repulsion_rate = 0.3,
                       maternal_geno_mode = "all_het",
                       paternal_pA_base = 0.4, error_rate = 0.01, seed = 7)
  dat <- make_dat(sim); mk <- sim$truth$markers_union
  oph <- oracle_phased(mk, 1L - sim$truth$v_true[[1]], "P1")

  # Two priors with identical q = 0.5 but different genotype split.
  piA <- matrix(c(0.50, 0.00, 0.50), 3, Tm)
  piB <- matrix(c(0.25, 0.50, 0.25), 3, Tm)
  mA <- hmm_map(dat, phased = oph, dam = 1, epsilon = 0.01, paternal_mode = "gametic",
                pi_prior_in = piA, lambda = 50, tol = 1e-6, maxit = 500)
  mB <- hmm_map(dat, phased = oph, dam = 1, epsilon = 0.01, paternal_mode = "gametic",
                pi_prior_in = piB, lambda = 50, tol = 1e-6, maxit = 500)

  expect_identical(mA$fit$r, mB$fit$r)
})

test_that("D = 1 joint equals single-dam under gametic", {
  RcppParallel::setThreadOptions(numThreads = 1)
  set.seed(11); Tm <- 40
  r_true <- runif(Tm - 1, 0.01, 0.15)
  sim <- sim_multi_pop(T_markers = Tm, n_pops = 1, n_ind_per_pop = 150,
                       marker_intersection = 1, r_vec = r_true,
                       phase_mode = "random", repulsion_rate = 0.3,
                       maternal_geno_mode = "all_het",
                       paternal_pA_base = 0.4, error_rate = 0.01, seed = 11)
  dat <- make_dat(sim); mk <- sim$truth$markers_union
  oph <- oracle_phased(mk, 1L - sim$truth$v_true[[1]], "P1")

  ms <- hmm_map(dat, phased = oph, dam = 1, epsilon = 0.01,
                paternal_mode = "gametic", tol = 1e-6, maxit = 500)
  mj <- hmm_map_joint(dat, phased = oph, dam = "all", epsilon = 0.01,
                      paternal_mode = "gametic", tol = 1e-6, maxit = 500)

  expect_equal(as.numeric(mj$fit$r), as.numeric(ms$fit$r), tolerance = 1e-8)
})

test_that("per_marker is disabled with an informative error", {
  set.seed(1); Tm <- 10
  sim <- sim_multi_pop(T_markers = Tm, n_pops = 1, n_ind_per_pop = 50,
                       marker_intersection = 1, r_const = 0.1,
                       phase_mode = "all_coupling", maternal_geno_mode = "all_het",
                       paternal_pA_base = 0.5, seed = 1)
  dat <- make_dat(sim); mk <- sim$truth$markers_union
  oph <- oracle_phased(mk, 1L - sim$truth$v_true[[1]], "P1")
  expect_error(
    hmm_map(dat, phased = oph, dam = 1, paternal_mode = "per_marker"),
    "not identifiable|deprecated")
})

test_that("two_locus is disabled with an informative error", {
  set.seed(1); Tm <- 10
  sim <- sim_multi_pop(T_markers = Tm, n_pops = 1, n_ind_per_pop = 50,
                       marker_intersection = 1, r_const = 0.1,
                       phase_mode = "all_coupling", maternal_geno_mode = "all_het",
                       paternal_pA_base = 0.5, seed = 1)
  dat <- make_dat(sim); mk <- sim$truth$markers_union
  oph <- oracle_phased(mk, 1L - sim$truth$v_true[[1]], "P1")
  expect_error(
    hmm_map(dat, phased = oph, dam = 1, paternal_mode = "two_locus"),
    "disabled|identif")
})

test_that("fit$q / fit$q_list are the canonical paternal outputs", {
  RcppParallel::setThreadOptions(numThreads = 1)
  set.seed(7); Tm <- 30
  r_true <- runif(Tm - 1, 0.01, 0.15)
  sim <- sim_multi_pop(T_markers = Tm, n_pops = 2, n_ind_per_pop = c(120, 100),
                       marker_intersection = 1, r_vec = r_true,
                       phase_mode = "random", repulsion_rate = 0.3,
                       maternal_geno_mode = "all_het",
                       paternal_pA_base = 0.4, paternal_pA_sd = 0.05,
                       error_rate = 0.01, seed = 7)
  dat <- make_dat(sim); mk <- sim$truth$markers_union

  # single dam
  oph1 <- oracle_phased(mk, 1L - sim$truth$v_true[[1]], names(sim$G_list)[1])
  m1 <- hmm_map(dat, phased = oph1, dam = 1, epsilon = 0.01, paternal_mode = "gametic")
  expect_false(is.null(m1$fit$q))
  expect_length(m1$fit$q, Tm)
  expect_true(all(m1$fit$q >= 0 & m1$fit$q <= 1))
  expect_equal(unname(m1$fit$q),
               unname(as.numeric(m1$fit$pi["AA", ] + 0.5 * m1$fit$pi["Aa", ])),
               tolerance = 1e-12)

  # joint: per-dam q_list
  oph <- oracle_multi(sim, mk)
  mj <- hmm_map(dat, phased = oph, dam = "all", epsilon = 0.01, paternal_mode = "gametic")
  expect_false(is.null(mj$fit$q_list))
  expect_length(mj$fit$q_list, 2L)
  expect_length(mj$fit$q_list[[1]], Tm)
  expect_true(all(mj$fit$q_list[[1]] >= 0 & mj$fit$q_list[[1]] <= 1))
  expect_equal(unname(mj$fit$q_list[[1]]),
               unname(as.numeric(mj$fit$pi_list[[1]]["AA", ] + 0.5 * mj$fit$pi_list[[1]]["Aa", ])),
               tolerance = 1e-12)
})
