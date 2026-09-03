# Tests for the DIRECT gametic paternal E-step / M-step (paternal_mode = "HWE").
#
# The public API and manuscript define a single gametic parameter q_k with the
# penalized (MAP) update
#   q = (N_A + alpha) / (N_A + N_a + alpha + beta),   alpha = lambda*q0, beta = lambda*(1-q0),
# where N_A, N_a are DIRECT expected transmitted paternal-gamete counts:
#   rho_A = sum_h gamma(h) q P(y|h,A) / b(h),   b(h) = q P(y|h,A) + (1-q) P(y|h,a).
# The engine previously updated q from expected diploid paternal-GENOTYPE
# responsibilities (N_A = Ng_AA + 0.5 Ng_Aa); a genotype carries two alleles, so
# with finite lambda that made the effective pseudocount strength ~2x the
# documented value. These tests pin the engine to the direct-gametic estimator.
#
# Shared helpers (make_dat, oracle_phased, one_dam, hwe_cols, pen_q) live in
# helper-sim.R.

test_that("one-marker analytic q: homozygous-dam gamete counts give the closed-form MAP", {
  RcppParallel::setThreadOptions(numThreads = 1)
  # Maternal AA at both markers => each offspring's paternal gamete is observed
  # directly (AA <=> paternal A, Aa <=> paternal a; y=0 impossible). With eps=0
  # the transmitted-gamete counts are exact and independent of the current q, so
  #   q_hat_t = (n_AA,t + alpha) / (n_AA,t + n_Aa,t + alpha + beta).
  M  <- c(2L, 2L)
  ph <- 1L
  y1 <- c(rep(2L, 70), rep(1L, 30))   # marker 1: 70 AA, 30 Aa
  y2 <- c(rep(2L, 40), rep(1L, 60))   # marker 2: 40 AA, 60 Aa
  G  <- cbind(y1, y2); storage.mode(G) <- "integer"
  q_of <- function(f) as.numeric(f$pi["AA", ] + 0.5 * f$pi["Aa", ])

  # lambda = 20, target q0 = 0.5  =>  alpha = beta = 10
  f20 <- HSMap:::hmm_hs_cpp_parallel(G, M, ph, r_start = 0.05, pi_mode = "HWE",
           pi_prior_in = hwe_cols(c(0.5, 0.5)), lambda = 20, epsilon = 0, tol = 1e-12,
           maxit = 500, paternal_mode = "HWE")
  expect_equal(q_of(f20), c((70 + 10)/(100 + 20), (40 + 10)/(100 + 20)), tolerance = 1e-9)

  # lambda = 0  =>  MLE = n_AA / n
  f0 <- HSMap:::hmm_hs_cpp_parallel(G, M, ph, r_start = 0.05, pi_mode = "HWE",
          pi_prior_in = hwe_cols(c(0.5, 0.5)), lambda = 0, epsilon = 0, tol = 1e-12,
          maxit = 500, paternal_mode = "HWE")
  expect_equal(q_of(f0), c(0.70, 0.40), tolerance = 1e-9)
})

test_that("engine q equals the coordinatewise maximizer of the penalized observed likelihood", {
  RcppParallel::setThreadOptions(numThreads = 1)
  lam <- 6; q0 <- 0.5
  d <- one_dam(2024L, Tm = 25L, n = 400L, pA = 0.30, err = 0.01)
  m <- hmm_map(d$dat, phased = d$oph, dam = 1, epsilon = 0.01, paternal_mode = "HWE",
               lambda = lam, tol = 1e-10, maxit = 3000)
  qhat <- as.numeric(m$fit$q); rhat <- as.numeric(m$fit$r)
  pen1 <- function(qt) lam * q0 * log(qt) + lam * (1 - q0) * log(1 - qt)  # marker-t penalty term
  # q_hat_t must maximize obsLL + penalty as a function of q_t, holding r_hat and
  # the other q's fixed (the EM stationarity condition == direct optimization).
  for (t in c(4L, 10L, 18L)) {
    f <- function(qt) {
      qv <- qhat; qv[t] <- qt
      HSMap:::loglik_hs_cpp(d$G, d$M, d$ph, rhat, hwe_cols(qv), 0.01) + pen1(qt)
    }
    opt <- optimize(f, c(1e-4, 1 - 1e-4), maximum = TRUE, tol = 1e-10)
    expect_equal(opt$maximum, qhat[t], tolerance = 2e-3)
  }
})

test_that("lambda=0 gametic fit matches the unregularized per_marker-induced q (and r)", {
  RcppParallel::setThreadOptions(numThreads = 1)
  d <- one_dam(55L, Tm = 30L, n = 400L, pA = 0.35, err = 0.01)
  # The observed likelihood depends on the paternal genotypes only through q, so
  # at lambda=0 the gametic (HWE) and free per_marker parameterizations must reach
  # the same identifiable q and the same shared r.
  fh <- HSMap:::hmm_hs_cpp_parallel(d$G, d$M, d$ph, r_start = 0.05, pi_mode = "HWE",
          pi_prior_in = NULL, lambda = 0, epsilon = 0.01, tol = 1e-10, maxit = 4000,
          paternal_mode = "HWE")
  fp <- HSMap:::hmm_hs_cpp_parallel(d$G, d$M, d$ph, r_start = 0.05, pi_mode = "per_marker",
          pi_prior_in = NULL, lambda = 0, epsilon = 0.01, tol = 1e-10, maxit = 4000,
          paternal_mode = "per_marker")
  qh <- as.numeric(fh$pi["AA", ] + 0.5 * fh$pi["Aa", ])
  qp <- as.numeric(fp$pi["AA", ] + 0.5 * fp$pi["Aa", ])
  expect_equal(qh, qp, tolerance = 5e-3)
  expect_equal(as.numeric(fh$r), as.numeric(fp$r), tolerance = 5e-3)
})

test_that("D=1 joint equals the single-dam gametic fit, including q (finite lambda)", {
  RcppParallel::setThreadOptions(numThreads = 1)
  d <- one_dam(11L, Tm = 35L, n = 200L, pA = 0.25, err = 0.01)
  ms <- hmm_map(d$dat, phased = d$oph, dam = 1, epsilon = 0.01, paternal_mode = "HWE",
                lambda = 6, tol = 1e-9, maxit = 1500)
  mj <- hmm_map_joint(d$dat, phased = d$oph, dam = "all", epsilon = 0.01,
                      paternal_mode = "HWE", lambda = 6, tol = 1e-9, maxit = 1500)
  expect_equal(as.numeric(mj$fit$r), as.numeric(ms$fit$r), tolerance = 1e-8)
  expect_equal(as.numeric(mj$fit$q_list[[1]]), as.numeric(ms$fit$q), tolerance = 1e-8)
})

test_that("penalized objective is non-decreasing across EM iterations", {
  RcppParallel::setThreadOptions(numThreads = 1)
  lam <- 10; q0 <- 0.5
  d <- one_dam(9001L, Tm = 40L, n = 300L, pA = 0.10, err = 0.01)  # extreme pA => q moves a lot
  obj <- numeric(20)
  for (k in 1:20) {
    f <- HSMap:::hmm_hs_cpp_parallel(d$G, d$M, d$ph, r_start = 0.05, pi_mode = "HWE",
           pi_prior_in = NULL, lambda = lam, epsilon = 0.01, tol = 0, maxit = k,  # tol=0 => exactly k iters
           paternal_mode = "HWE")
    q <- as.numeric(f$pi["AA", ] + 0.5 * f$pi["Aa", ])
    obj[k] <- HSMap:::loglik_hs_cpp(d$G, d$M, d$ph, as.numeric(f$r), hwe_cols(q), 0.01) +
              pen_q(q, lam, q0)
  }
  expect_true(all(diff(obj) >= -1e-6))
})
