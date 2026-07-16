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

test_that("per_marker is deprecated: warns and routes to gametic", {
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

  mg <- hmm_map(dat, phased = oph, dam = 1, epsilon = 0.01,
                paternal_mode = "gametic", tol = 1e-6, maxit = 500)
  expect_warning(
    mp <- hmm_map(dat, phased = oph, dam = 1, epsilon = 0.01,
                  paternal_mode = "per_marker", tol = 1e-6, maxit = 500),
    "deprecated")
  expect_identical(mp$fit$r, mg$fit$r)   # routed to gametic -> byte-identical
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

# ---- M1.1: first-class gametic pseudocount priors -----------------------------

test_that("gametic recovers q near 0 and near 1 (no shrinkage)", {
  RcppParallel::setThreadOptions(numThreads = 1)
  Tm <- 30; r_true <- rep(0.05, Tm - 1)
  fit_mean_q <- function(pA, seed) {
    sim <- sim_multi_pop(T_markers = Tm, n_pops = 1, n_ind_per_pop = 400,
                         marker_intersection = 1, r_vec = r_true,
                         phase_mode = "all_coupling", maternal_geno_mode = "all_het",
                         paternal_pA_base = pA, error_rate = 0, seed = seed)
    dat <- make_dat(sim); mk <- sim$truth$markers_union
    oph <- oracle_phased(mk, 1L - sim$truth$v_true[[1]], "P1")
    m <- hmm_map(dat, phased = oph, dam = 1, epsilon = 1e-6, paternal_mode = "gametic",
                 q_prior_in = list(alpha = 0, beta = 0), tol = 1e-7, maxit = 800)
    expect_true(all(m$fit$q >= 0 & m$fit$q <= 1))
    mean(m$fit$q)
  }
  expect_lt(fit_mean_q(0.05, 3), 0.15)   # true q ~ 0.05
  expect_gt(fit_mean_q(0.95, 4), 0.85)   # true q ~ 0.95
})

test_that("uninformative markers fall back to the q prior mean", {
  RcppParallel::setThreadOptions(numThreads = 1)
  set.seed(5); Tm <- 20; r_true <- rep(0.05, Tm - 1)
  sim <- sim_multi_pop(T_markers = Tm, n_pops = 1, n_ind_per_pop = 200,
                       marker_intersection = 1, r_vec = r_true,
                       phase_mode = "all_coupling", maternal_geno_mode = "all_het",
                       paternal_pA_base = 0.4, error_rate = 0, seed = 5)
  dat <- make_dat(sim); mk <- sim$truth$markers_union
  dat$G_list[[1]][, 10] <- NA_integer_          # marker 10 uninformative for q
  oph <- oracle_phased(mk, 1L - sim$truth$v_true[[1]], "P1")
  m <- hmm_map(dat, phased = oph, dam = 1, epsilon = 1e-3, paternal_mode = "gametic",
               q_prior_in = 0.3, lambda = 40, tol = 1e-7, maxit = 800)
  expect_equal(m$fit$q[[10]], 0.3, tolerance = 1e-6)
})

test_that("no shrinkage (lambda = 0): the prior mean is ignored at convergence", {
  RcppParallel::setThreadOptions(numThreads = 1)
  set.seed(9); Tm <- 25; r_true <- runif(Tm - 1, 0.02, 0.12)
  sim <- sim_multi_pop(T_markers = Tm, n_pops = 1, n_ind_per_pop = 300,
                       marker_intersection = 1, r_vec = r_true,
                       phase_mode = "random", repulsion_rate = 0.3,
                       maternal_geno_mode = "all_het",
                       paternal_pA_base = 0.4, error_rate = 0.01, seed = 9)
  dat <- make_dat(sim); mk <- sim$truth$markers_union
  oph <- oracle_phased(mk, 1L - sim$truth$v_true[[1]], "P1")
  m_a <- hmm_map(dat, phased = oph, dam = 1, epsilon = 0.01, paternal_mode = "gametic",
                 q_prior_in = 0.2, lambda = 0, tol = 1e-8, maxit = 1500)
  m_b <- hmm_map(dat, phased = oph, dam = 1, epsilon = 0.01, paternal_mode = "gametic",
                 q_prior_in = 0.8, lambda = 0, tol = 1e-8, maxit = 1500)
  expect_equal(as.numeric(m_a$fit$r), as.numeric(m_b$fit$r), tolerance = 1e-4)
  expect_equal(as.numeric(m_a$fit$q), as.numeric(m_b$fit$q), tolerance = 1e-4)
})

test_that("marker-specific q priors set uninformative markers to their own means", {
  RcppParallel::setThreadOptions(numThreads = 1)
  set.seed(6); Tm <- 20; r_true <- rep(0.05, Tm - 1)
  sim <- sim_multi_pop(T_markers = Tm, n_pops = 1, n_ind_per_pop = 150,
                       marker_intersection = 1, r_vec = r_true,
                       phase_mode = "all_coupling", maternal_geno_mode = "all_het",
                       paternal_pA_base = 0.5, error_rate = 0, seed = 6)
  dat <- make_dat(sim); mk <- sim$truth$markers_union
  dat$G_list[[1]][, 5]  <- NA_integer_
  dat$G_list[[1]][, 15] <- NA_integer_
  oph <- oracle_phased(mk, 1L - sim$truth$v_true[[1]], "P1")
  means <- rep(0.5, Tm); means[5] <- 0.2; means[15] <- 0.8
  conc <- 50
  qp <- list(alpha = conc * means, beta = conc * (1 - means))   # constant concentration
  fit <- hmm_map(dat, phased = oph, dam = 1, epsilon = 1e-3, paternal_mode = "gametic",
                 q_prior_in = qp, tol = 1e-7, maxit = 800)
  expect_equal(fit$fit$q[[5]],  0.2, tolerance = 1e-6)
  expect_equal(fit$fit$q[[15]], 0.8, tolerance = 1e-6)
})

test_that("per-marker total pseudocount (alpha + beta) must be constant across markers", {
  RcppParallel::setThreadOptions(numThreads = 1)
  set.seed(1); Tm <- 10
  sim <- sim_multi_pop(T_markers = Tm, n_pops = 1, n_ind_per_pop = 50,
                       marker_intersection = 1, r_const = 0.1,
                       phase_mode = "all_coupling", maternal_geno_mode = "all_het",
                       paternal_pA_base = 0.5, seed = 1)
  dat <- make_dat(sim); mk <- sim$truth$markers_union
  oph <- oracle_phased(mk, 1L - sim$truth$v_true[[1]], "P1")
  bad <- list(alpha = c(5, rep(10, Tm - 1)), beta = rep(10, Tm))  # alpha+beta not constant
  expect_error(
    hmm_map(dat, phased = oph, dam = 1, paternal_mode = "gametic", q_prior_in = bad),
    "pseudocount")
})

test_that("legacy pi priors with the same induced q give identical results", {
  RcppParallel::setThreadOptions(numThreads = 1)
  set.seed(7); Tm <- 25; r_true <- runif(Tm - 1, 0.01, 0.15)
  sim <- sim_multi_pop(T_markers = Tm, n_pops = 1, n_ind_per_pop = 200,
                       marker_intersection = 1, r_vec = r_true,
                       phase_mode = "random", repulsion_rate = 0.3,
                       maternal_geno_mode = "all_het",
                       paternal_pA_base = 0.4, error_rate = 0.01, seed = 7)
  dat <- make_dat(sim); mk <- sim$truth$markers_union
  oph <- oracle_phased(mk, 1L - sim$truth$v_true[[1]], "P1")
  piA <- matrix(c(0.50, 0.00, 0.50), 3, Tm)   # induced q = 0.5
  piB <- matrix(c(0.25, 0.50, 0.25), 3, Tm)   # induced q = 0.5
  mA <- suppressWarnings(hmm_map(dat, phased = oph, dam = 1, epsilon = 0.01,
          paternal_mode = "per_marker", pi_prior_in = piA, tol = 1e-6, maxit = 500))
  mB <- suppressWarnings(hmm_map(dat, phased = oph, dam = 1, epsilon = 0.01,
          paternal_mode = "per_marker", pi_prior_in = piB, tol = 1e-6, maxit = 500))
  expect_identical(mA$fit$r, mB$fit$r)
})

# ---- M1.2: pseudocount (MAP) interpretation and the mild default --------------

test_that("the default gametic prior is alpha = beta = 10 (historical lambda = 20)", {
  RcppParallel::setThreadOptions(numThreads = 1)
  set.seed(7); Tm <- 30; r_true <- runif(Tm - 1, 0.01, 0.15)
  sim <- sim_multi_pop(T_markers = Tm, n_pops = 1, n_ind_per_pop = 200,
                       marker_intersection = 1, r_vec = r_true,
                       phase_mode = "random", repulsion_rate = 0.3,
                       maternal_geno_mode = "all_het",
                       paternal_pA_base = 0.4, error_rate = 0.01, seed = 7)
  dat <- make_dat(sim); mk <- sim$truth$markers_union
  oph <- oracle_phased(mk, 1L - sim$truth$v_true[[1]], "P1")
  m_def <- hmm_map(dat, phased = oph, dam = 1, epsilon = 0.01,
                   paternal_mode = "gametic", tol = 1e-6, maxit = 500)
  m_10  <- hmm_map(dat, phased = oph, dam = 1, epsilon = 0.01, paternal_mode = "gametic",
                   q_prior_in = list(alpha = 10, beta = 10), tol = 1e-6, maxit = 500)
  expect_identical(m_def$fit$r, m_10$fit$r)
  expect_identical(m_def$fit$q, m_10$fit$q)
})

test_that("alpha = beta = 0 reproduces the unpenalized engine update", {
  RcppParallel::setThreadOptions(numThreads = 1)
  set.seed(9); Tm <- 25; r_true <- runif(Tm - 1, 0.02, 0.12)
  sim <- sim_multi_pop(T_markers = Tm, n_pops = 1, n_ind_per_pop = 300,
                       marker_intersection = 1, r_vec = r_true,
                       phase_mode = "random", repulsion_rate = 0.3,
                       maternal_geno_mode = "all_het",
                       paternal_pA_base = 0.4, error_rate = 0.01, seed = 9)
  dat <- make_dat(sim); mk <- sim$truth$markers_union
  oph <- oracle_phased(mk, 1L - sim$truth$v_true[[1]], "P1")
  G <- dat$G_list[[1]][, mk, drop = FALSE]; storage.mode(G) <- "integer"
  M <- as.integer(dat$M_list[[1]][mk]); ph <- as.integer(1L - sim$truth$v_true[[1]])
  m0 <- hmm_map(dat, phased = oph, dam = 1, epsilon = 0.01, paternal_mode = "gametic",
                q_prior_in = list(alpha = 0, beta = 0), r_start = 0.05, tol = 1e-7, maxit = 1000)
  eng <- HSMap:::hmm_hs_cpp_parallel(G, M, ph, r_start = 0.05, pi_mode = "HWE",
             pi_prior_in = NULL, lambda = 0, epsilon = 0.01, tol = 1e-7, maxit = 1000,
             paternal_mode = "HWE")
  expect_equal(as.numeric(m0$fit$r), as.numeric(eng$r), tolerance = 1e-8)
})

test_that("the equivalent probability prior is Beta(alpha + 1, beta + 1)", {
  # The M-step maximizes N_A log q + N_a log(1-q) + alpha log q + beta log(1-q).
  N_A <- 7; N_a <- 13; alpha <- 3; beta <- 5
  update <- (N_A + alpha) / (N_A + N_a + alpha + beta)

  pen_obj <- function(q) N_A*log(q) + N_a*log(1 - q) + alpha*log(q) + beta*log(1 - q)
  argmax  <- stats::optimize(pen_obj, c(1e-9, 1 - 1e-9), maximum = TRUE, tol = 1e-12)$maximum
  expect_equal(update, argmax, tolerance = 1e-5)

  # mode of the Beta(N_A + alpha + 1, N_a + beta + 1) posterior == the update
  a_post <- N_A + alpha + 1; b_post <- N_a + beta + 1
  expect_equal(update, (a_post - 1) / (a_post + b_post - 2), tolerance = 1e-12)

  # treating alpha,beta as ordinary Beta shape params (penalty alpha-1,beta-1)
  # gives a different value -> the two interpretations are distinguishable
  wrong <- (N_A + alpha - 1) / (N_A + N_a + (alpha - 1) + (beta - 1))
  expect_false(isTRUE(all.equal(update, wrong)))
})

test_that("lambda = 20 equals alpha = beta = 10, and the default applies shrinkage", {
  RcppParallel::setThreadOptions(numThreads = 1)
  set.seed(7); Tm <- 30; r_true <- runif(Tm - 1, 0.01, 0.15)
  sim <- sim_multi_pop(T_markers = Tm, n_pops = 1, n_ind_per_pop = 200,
                       marker_intersection = 1, r_vec = r_true,
                       phase_mode = "random", repulsion_rate = 0.3,
                       maternal_geno_mode = "all_het",
                       paternal_pA_base = 0.4, error_rate = 0.01, seed = 7)
  dat <- make_dat(sim); mk <- sim$truth$markers_union
  oph <- oracle_phased(mk, 1L - sim$truth$v_true[[1]], "P1")
  m_20 <- hmm_map(dat, phased = oph, dam = 1, epsilon = 0.01, paternal_mode = "gametic",
                  lambda = 20, tol = 1e-6, maxit = 500)
  m_10 <- hmm_map(dat, phased = oph, dam = 1, epsilon = 0.01, paternal_mode = "gametic",
                  q_prior_in = list(alpha = 10, beta = 10), tol = 1e-6, maxit = 500)
  expect_identical(m_20$fit$r, m_10$fit$r)
  expect_identical(m_20$fit$q, m_10$fit$q)
  # the regularized default differs from no regularization (alpha = beta = 0)
  m_0 <- hmm_map(dat, phased = oph, dam = 1, epsilon = 0.01, paternal_mode = "gametic",
                 q_prior_in = list(alpha = 0, beta = 0), tol = 1e-6, maxit = 500)
  expect_false(isTRUE(all.equal(m_20$fit$q, m_0$fit$q)))
})

# ---- q_prior_list per-dam completeness in hmm_map_joint() ---------------------

mk_joint3 <- function(seed = 22, Tm = 30) {
  set.seed(seed); r_true <- runif(Tm - 1, 0.01, 0.15)
  sim <- sim_multi_pop(T_markers = Tm, n_pops = 3, n_ind_per_pop = c(80, 70, 60),
                       marker_intersection = 1, r_vec = r_true,
                       phase_mode = "random", repulsion_rate = 0.3,
                       maternal_geno_mode = "all_het",
                       paternal_pA_base = 0.4, paternal_pA_sd = 0.05,
                       error_rate = 0.01, seed = seed)
  mk <- sim$truth$markers_union
  list(dat = make_dat(sim), oph = oracle_multi(sim, mk), dn = names(sim$G_list))
}

test_that("q_prior_list: a single shared spec applies to all dams", {
  RcppParallel::setThreadOptions(numThreads = 1)
  j <- mk_joint3()
  m_shared <- hmm_map_joint(j$dat, phased = j$oph, dam = "all", epsilon = 0.01,
                            paternal_mode = "gametic", q_prior_list = list(alpha = 3, beta = 3),
                            tol = 1e-6, maxit = 300)
  per <- stats::setNames(rep(list(list(alpha = 3, beta = 3)), 3), j$dn)
  m_perdam <- hmm_map_joint(j$dat, phased = j$oph, dam = "all", epsilon = 0.01,
                            paternal_mode = "gametic", q_prior_list = per,
                            tol = 1e-6, maxit = 300)
  expect_equal(as.numeric(m_shared$fit$r), as.numeric(m_perdam$fit$r), tolerance = 1e-10)
  # a single numeric mean also applies to all dams
  m_num <- hmm_map_joint(j$dat, phased = j$oph, dam = "all", epsilon = 0.01,
                         paternal_mode = "gametic", q_prior_list = 0.3, tol = 1e-6, maxit = 200)
  expect_length(m_num$fit$q_list, 3L)
})

test_that("q_prior_list: complete named and positional per-dam lists work and agree", {
  RcppParallel::setThreadOptions(numThreads = 1)
  j <- mk_joint3()
  named <- stats::setNames(list(list(alpha = 2, beta = 8), list(alpha = 5, beta = 5),
                                list(alpha = 8, beta = 2)), j$dn)   # common concentration = 10
  m_named <- hmm_map_joint(j$dat, phased = j$oph, dam = "all", epsilon = 0.01,
                           paternal_mode = "gametic", q_prior_list = named, tol = 1e-6, maxit = 200)
  expect_length(m_named$fit$q_list, 3L)
  m_pos <- hmm_map_joint(j$dat, phased = j$oph, dam = "all", epsilon = 0.01,
                         paternal_mode = "gametic", q_prior_list = unname(named), tol = 1e-6, maxit = 200)
  expect_equal(as.numeric(m_named$fit$r), as.numeric(m_pos$fit$r), tolerance = 1e-10)
})

test_that("q_prior_list: named per-dam list missing a requested dam errors", {
  RcppParallel::setThreadOptions(numThreads = 1)
  j <- mk_joint3()
  bad <- stats::setNames(list(list(alpha = 3, beta = 3), list(alpha = 3, beta = 3)), j$dn[1:2])
  expect_error(
    hmm_map_joint(j$dat, phased = j$oph, dam = "all", paternal_mode = "gametic", q_prior_list = bad),
    "missing prior")
})

test_that("q_prior_list: named per-dam list with an unknown dam errors", {
  RcppParallel::setThreadOptions(numThreads = 1)
  j <- mk_joint3()
  # all requested dams ARE present, plus an extra unknown name -> unknown-dam error
  bad <- stats::setNames(rep(list(list(alpha = 3, beta = 3)), 4), c(j$dn, "NOPE"))
  expect_error(
    hmm_map_joint(j$dat, phased = j$oph, dam = "all", paternal_mode = "gametic", q_prior_list = bad),
    "unknown dam")
})

test_that("q_prior_list: unnamed list of the wrong length errors", {
  RcppParallel::setThreadOptions(numThreads = 1)
  j <- mk_joint3()
  too_few <- list(list(alpha = 3, beta = 3), list(alpha = 3, beta = 3))
  expect_error(
    hmm_map_joint(j$dat, phased = j$oph, dam = "all", paternal_mode = "gametic", q_prior_list = too_few),
    "one spec per requested dam")
  too_many <- rep(list(list(alpha = 3, beta = 3)), 4)
  expect_error(
    hmm_map_joint(j$dat, phased = j$oph, dam = "all", paternal_mode = "gametic", q_prior_list = too_many),
    "one spec per requested dam")
})
