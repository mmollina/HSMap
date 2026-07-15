# Regression tests for the paternal responsibility E-step in the HSMap EM.
#
# The single-locus paternal responsibility must normalize the emission PER
# maternal state h before averaging over the maternal-state posterior gamma:
#   R_g = sum_h gamma(h) * pi_g P(y|h,g) / b(h),   b(h) = sum_g pi_g P(y|h,g).
# The previous code normalized by the gamma-averaged emission, biasing the q
# update so the EM did not ascend the observed-data likelihood.
#
# Shared helpers (make_dat, oracle_phased, oracle_multi) live in helper-sim.R.

hwe_cols <- function(q) { q <- pmin(pmax(q, 1e-9), 1 - 1e-9); rbind(q^2, 2*q*(1-q), (1-q)^2) }
pen_q    <- function(q, lambda, q0 = 0.5) {   # penalty added to obs-LL: a log q + b log(1-q)
  a <- lambda * q0; b <- lambda * (1 - q0)
  sum(a * log(pmin(pmax(q, 1e-12), 1)) + b * log(pmin(pmax(1 - q, 1e-12), 1)))
}
one_dam <- function(seed, Tm = 60L, n = 300L, pA = 0.40, err = 0.01, rr = 0.3, mat = "all_het") {
  set.seed(seed); r_true <- runif(Tm - 1, 0.01, 0.15)
  sim <- sim_multi_pop(T_markers = Tm, n_pops = 1, n_ind_per_pop = n, marker_intersection = 1,
                       r_vec = r_true, phase_mode = "random", repulsion_rate = rr,
                       maternal_geno_mode = mat, maternal_pA = 0.5,
                       paternal_pA_base = pA, error_rate = err, seed = seed)
  mk <- sim$truth$markers_union
  dat <- make_dat(sim)
  oph <- oracle_phased(mk, 1L - sim$truth$v_true[[1]], "P1")
  G <- dat$G_list[[1]][, mk, drop = FALSE]; storage.mode(G) <- "integer"
  list(dat = dat, oph = oph, G = G, M = as.integer(dat$M_list[[1]][mk]),
       ph = as.integer(1L - sim$truth$v_true[[1]]), r_true = r_true, mk = mk)
}

test_that("paternal responsibilities sum to one (corrected per-state normalization)", {
  Rg <- function(ga, pi3, pyh0, pyh1) {
    b0 <- sum(pi3 * pyh0); b1 <- sum(pi3 * pyh1)
    r <- numeric(3)
    if (ga[1] > 0 && b0 > 0) r <- r + ga[1] * pi3 * pyh0 / b0
    if (ga[2] > 0 && b1 > 0) r <- r + ga[2] * pi3 * pyh1 / b1
    r
  }
  set.seed(1)
  for (i in 1:300) {
    ga  <- runif(2); ga <- ga / sum(ga)          # split maternal posterior
    pi3 <- runif(3); pi3 <- pi3 / sum(pi3)
    expect_equal(sum(Rg(ga, pi3, runif(3), runif(3))), 1, tolerance = 1e-12)
  }
  # degenerate gamma still sums to one
  expect_equal(sum(Rg(c(1, 0), c(.3,.4,.3), runif(3), runif(3))), 1, tolerance = 1e-12)
})

test_that("lambda=0 q-update equals the correct gamma-weighted paternal A count", {
  RcppParallel::setThreadOptions(numThreads = 1)
  d <- one_dam(101L, Tm = 40L, n = 400L, pA = 0.30, err = 1e-6)
  m <- hmm_map(d$dat, phased = d$oph, dam = 1, epsilon = 1e-6, paternal_mode = "HWE",
               lambda = 0, tol = 1e-9, maxit = 4000)
  q_fit <- as.numeric(m$fit$pi["AA", ] + 0.5 * m$fit$pi["Aa", ])
  g  <- gamma_cpp(d$G, d$M, d$ph, as.numeric(m$fit$r), m$fit$pi, 1e-6)  # [2, T, n]
  ha <- attr(g, "hap0_is_A")
  Tm <- length(d$mk); q_chk <- numeric(Tm)
  for (t in seq_len(Tm)) {
    p_a <- if (ha[t] == 1) g[2, t, ] else g[1, t, ]   # P(mom transmits a) = gamma(h=0)
    y <- d$G[, t]; obs <- !is.na(y)
    N_A <- sum(y[obs] == 2) + sum(p_a[obs][y[obs] == 1])   # dad-A count
    q_chk[t] <- N_A / sum(obs)
  }
  expect_equal(q_fit, q_chk, tolerance = 2e-3)
})

test_that("unregularized (lambda=0) EM warm-started from a regularized fit does not lower obs-LL", {
  RcppParallel::setThreadOptions(numThreads = 1)
  for (cfg in list(list(seed = 8001L, pA = 0.10), list(seed = 7L, pA = 0.40))) {
    d <- one_dam(cfg$seed, pA = cfg$pA)
    m1 <- hmm_map(d$dat, phased = d$oph, dam = 1, epsilon = 0.01, paternal_mode = "HWE",
                  lambda = 2, tol = 1e-8, maxit = 2000)
    q1 <- as.numeric(m1$fit$pi["AA", ] + 0.5 * m1$fit$pi["Aa", ]); r1 <- as.numeric(m1$fit$r)
    ll_start <- HSMap:::loglik_hs_cpp(d$G, d$M, d$ph, r1, hwe_cols(q1), 0.01)
    f <- HSMap:::hmm_hs_cpp_parallel(d$G, d$M, d$ph, r_start = 0.05, pi_mode = "HWE",
            pi_prior_in = hwe_cols(q1), lambda = 0, epsilon = 0.01, tol = 1e-10,
            maxit = 2000, paternal_mode = "HWE", Pi_prior_in = NULL, r_init = r1)
    ll_final <- HSMap:::loglik_hs_cpp(d$G, d$M, d$ph, as.numeric(f$r), f$pi, 0.01)
    expect_gte(ll_final, ll_start - 1e-6)
  }
})

test_that("penalized (lambda>0) objective does not decrease when warm-started", {
  RcppParallel::setThreadOptions(numThreads = 1)
  lam <- 4
  d <- one_dam(202L, pA = 0.20)
  # start from the unregularized solution, then run the penalized EM
  m0 <- hmm_map(d$dat, phased = d$oph, dam = 1, epsilon = 0.01, paternal_mode = "HWE",
                lambda = 0, tol = 1e-8, maxit = 2000)
  q0 <- as.numeric(m0$fit$pi["AA", ] + 0.5 * m0$fit$pi["Aa", ]); r0 <- as.numeric(m0$fit$r)
  obj_start <- HSMap:::loglik_hs_cpp(d$G, d$M, d$ph, r0, hwe_cols(q0), 0.01) + pen_q(q0, lam)
  f <- HSMap:::hmm_hs_cpp_parallel(d$G, d$M, d$ph, r_start = 0.05, pi_mode = "HWE",
          pi_prior_in = hwe_cols(rep(0.5, length(d$mk))), lambda = lam, epsilon = 0.01,
          tol = 1e-10, maxit = 2000, paternal_mode = "HWE", Pi_prior_in = NULL, r_init = r0)
  qf <- as.numeric(f$pi["AA", ] + 0.5 * f$pi["Aa", ])
  obj_final <- HSMap:::loglik_hs_cpp(d$G, d$M, d$ph, as.numeric(f$r), f$pi, 0.01) + pen_q(qf, lam)
  expect_gte(obj_final, obj_start - 1e-6)
})

test_that("impossible maternal state at epsilon=0 yields no NaN/Inf/division-by-zero", {
  RcppParallel::setThreadOptions(numThreads = 1)
  d <- one_dam(303L, Tm = 30L, n = 150L, err = 0)
  M <- d$M; G <- d$G
  M[10] <- 2L                    # force marker 10 maternal AA (homozygous)
  G[1:5, 10] <- 0L               # inject offspring 'aa' at marker 10: impossible (mom AA -> aa)
  f <- HSMap:::hmm_hs_cpp_parallel(G, M, d$ph, r_start = 0.05, pi_mode = "HWE",
          pi_prior_in = NULL, lambda = 0, epsilon = 0, tol = 1e-8, maxit = 500,
          paternal_mode = "HWE")
  expect_true(all(is.finite(as.numeric(f$r))))
  expect_true(all(is.finite(as.numeric(f$pi))))
  expect_true(is.finite(f$logLik))
})

test_that("D=1 joint equals the single-dam fit (shared corrected worker)", {
  RcppParallel::setThreadOptions(numThreads = 1)
  d <- one_dam(11L, Tm = 40L, n = 150L)
  ms <- hmm_map(d$dat, phased = d$oph, dam = 1, epsilon = 0.01, paternal_mode = "HWE",
                tol = 1e-6, maxit = 500)
  mj <- hmm_map_joint(d$dat, phased = d$oph, dam = "all", epsilon = 0.01,
                      paternal_mode = "HWE", tol = 1e-6, maxit = 500)
  expect_equal(as.numeric(mj$fit$r), as.numeric(ms$fit$r), tolerance = 1e-8)
})

test_that("final observed-data log-likelihood matches loglik_hs_cpp at returned params", {
  RcppParallel::setThreadOptions(numThreads = 1)
  d <- one_dam(404L, pA = 0.35)
  m <- hmm_map(d$dat, phased = d$oph, dam = 1, epsilon = 0.01, paternal_mode = "HWE",
               lambda = 2, tol = 1e-8, maxit = 2000)
  q <- as.numeric(m$fit$pi["AA", ] + 0.5 * m$fit$pi["Aa", ])
  ll <- HSMap:::loglik_hs_cpp(d$G, d$M, d$ph, as.numeric(m$fit$r), hwe_cols(q), 0.01)
  expect_true(is.finite(ll))
  expect_equal(ll, m$fit$logLik, tolerance = 1e-3)   # engine LL == recomputed LL at convergence
})
