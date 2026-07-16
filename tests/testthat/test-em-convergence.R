# Milestone 4: multipoint EM convergence and likelihood reporting.
#   * returned logLik is evaluated at the FINAL parameters,
#   * RELATIVE likelihood convergence + stable r (not the paternal decomposition),
#   * converged / iters / conv_reason / traces correct; iters never exceeds maxit.
# Shared helpers (one_dam, hwe_cols) live in helper-sim.R.

test_that("returned logLik matches an independent likelihood call at final parameters", {
  RcppParallel::setThreadOptions(numThreads = 1)
  d <- one_dam(11L, Tm = 30L, n = 250L, pA = 0.35)
  for (lam in c(0, 20)) {
    m <- hmm_map(d$dat, phased = d$oph, dam = 1, epsilon = 0.01,
                 paternal_mode = "HWE", lambda = lam, tol = 1e-7, maxit = 2000)
    q <- as.numeric(m$fit$q)
    ll <- HSMap:::loglik_hs_cpp(d$G, d$M, d$ph, as.numeric(m$fit$r), hwe_cols(q), 0.01)
    expect_equal(m$fit$logLik, ll, tolerance = 1e-6)   # engine logLik == recomputed at final params
    expect_true(is.finite(m$fit$logLik))
  }
})

test_that("convergence is relative (not absolute) and requires stable r", {
  RcppParallel::setThreadOptions(numThreads = 1)
  d <- one_dam(7L, Tm = 40L, n = 300L, pA = 0.40)
  m <- hmm_map(d$dat, phased = d$oph, dam = 1, epsilon = 0.01, paternal_mode = "HWE",
               lambda = 0, tol = 1e-6, maxit = 3000)
  expect_true(m$fit$converged)
  tr <- m$fit$loglik_trace
  n <- length(tr)
  rel_last <- abs(tr[n] - tr[n - 1]) / (1 + abs(tr[n - 1]))
  expect_lt(rel_last, 1e-6)                              # relative criterion satisfied
  # the last r change is below tol (stable r required for convergence)
  expect_lt(utils::tail(m$fit$max_dr_trace, 1), 1e-6)
  # a purely absolute |dLL| < tol would be MUCH stricter here (|logLik| ~ 1e4)
  expect_gt(1 + abs(tr[n - 1]), 100)
})

test_that("converged / iters / conv_reason are correct, and iters never exceeds maxit", {
  RcppParallel::setThreadOptions(numThreads = 1)
  d <- one_dam(3L, Tm = 25L, n = 200L)
  # well-behaved fit -> converges with the documented reason
  m <- hmm_map(d$dat, phased = d$oph, dam = 1, epsilon = 0.01, paternal_mode = "HWE",
               tol = 1e-6, maxit = 3000)
  expect_true(m$fit$converged)
  expect_identical(m$fit$conv_reason, "relative_loglik_and_r_stable")
  expect_lte(m$fit$iters, 3000L)
  # deliberately low maxit -> not converged, reason maxit_reached, iters == maxit
  m2 <- hmm_map(d$dat, phased = d$oph, dam = 1, epsilon = 0.01, paternal_mode = "HWE",
                tol = 1e-12, maxit = 3L)
  expect_false(m2$fit$converged)
  expect_identical(m2$fit$conv_reason, "maxit_reached")
  expect_identical(m2$fit$iters, 3L)                    # exactly maxit, no off-by-one
  expect_length(m2$fit$loglik_trace, 3L)
})

test_that("loglik trace terminates at (close to) the reported final logLik", {
  RcppParallel::setThreadOptions(numThreads = 1)
  d <- one_dam(9L, Tm = 30L, n = 250L, pA = 0.3)
  m <- hmm_map(d$dat, phased = d$oph, dam = 1, epsilon = 0.01, paternal_mode = "HWE",
               lambda = 0, tol = 1e-8, maxit = 3000)
  tr <- m$fit$loglik_trace
  # at convergence the last trace value ~ the final logLik (recomputed post-M-step)
  expect_lt(abs(utils::tail(tr, 1) - m$fit$logLik) / (1 + abs(m$fit$logLik)), 1e-5)
})

test_that("single-dam and D=1 joint agree on logLik and convergence", {
  RcppParallel::setThreadOptions(numThreads = 1)
  d <- one_dam(11L, Tm = 30L, n = 150L)
  ms <- hmm_map(d$dat, phased = d$oph, dam = 1, epsilon = 0.01, paternal_mode = "HWE",
                lambda = 6, tol = 1e-8, maxit = 2000)
  mj <- hmm_map_joint(d$dat, phased = d$oph, dam = "all", epsilon = 0.01,
                      paternal_mode = "HWE", lambda = 6, tol = 1e-8, maxit = 2000)
  expect_equal(as.numeric(mj$fit$r), as.numeric(ms$fit$r), tolerance = 1e-7)
  expect_equal(mj$fit$logLik, ms$fit$logLik, tolerance = 1e-5)
})

test_that("penalized objective is reported for finite-prior gametic fits and NA otherwise", {
  RcppParallel::setThreadOptions(numThreads = 1)
  d <- one_dam(5L, Tm = 20L, n = 200L, pA = 0.3)
  m20 <- hmm_map(d$dat, phased = d$oph, dam = 1, epsilon = 0.01, paternal_mode = "HWE",
                 lambda = 20, tol = 1e-7, maxit = 2000)
  m0  <- hmm_map(d$dat, phased = d$oph, dam = 1, epsilon = 0.01, paternal_mode = "HWE",
                 lambda = 0,  tol = 1e-7, maxit = 2000)
  expect_true(is.finite(m20$fit$penalized_obj))
  expect_lt(m20$fit$penalized_obj, m20$fit$logLik)      # penalty is negative
  expect_true(is.na(m0$fit$penalized_obj))              # no active penalty
})
