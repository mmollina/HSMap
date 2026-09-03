# Commit 1: EM active objective and convergence.
#   * active objective = obsLL (lambda=0) or obsLL + q-penalty (gametic, lambda>0),
#   * convergence on relative active objective AND stable identifiable r and q,
#   * traces end exactly at the returned quantities; monotone objective,
#   * non-convergence warns. Shared helpers (one_dam, hwe_cols) in helper-sim.R.

test_that("returned logLik matches an independent likelihood call at final parameters", {
  RcppParallel::setThreadOptions(numThreads = 1)
  d <- one_dam(11L, Tm = 30L, n = 250L, pA = 0.35)
  for (lam in c(0, 20)) {
    m <- hmm_map(d$dat, phased = d$oph, dam = 1, epsilon = 0.01,
                 paternal_mode = "HWE", lambda = lam, tol = 1e-7, maxit = 3000)
    q <- as.numeric(m$fit$q)
    ll <- HSMap:::loglik_hs_cpp(d$G, d$M, d$ph, as.numeric(m$fit$r), hwe_cols(q), 0.01)
    expect_equal(m$fit$logLik, ll, tolerance = 1e-6)
  }
})

test_that("active objective is obsLL for lambda=0 and penalized for lambda>0", {
  RcppParallel::setThreadOptions(numThreads = 1)
  d <- one_dam(7L, Tm = 25L, n = 220L, pA = 0.3)
  m0  <- hmm_map(d$dat, phased = d$oph, dam = 1, epsilon = 0.01, paternal_mode = "HWE",
                 lambda = 0,  tol = 1e-7, maxit = 3000)
  m20 <- hmm_map(d$dat, phased = d$oph, dam = 1, epsilon = 0.01, paternal_mode = "HWE",
                 lambda = 20, tol = 1e-7, maxit = 3000)
  expect_equal(m0$fit$objective,  m0$fit$logLik,        tolerance = 1e-9)  # lambda=0 -> obsLL
  expect_true(is.na(m0$fit$penalized_obj))
  expect_equal(m20$fit$objective, m20$fit$penalized_obj, tolerance = 1e-9) # lambda>0 -> penalized
  expect_lt(m20$fit$penalized_obj, m20$fit$logLik)                          # penalty < 0
})

test_that("all traces end exactly at the returned quantities", {
  RcppParallel::setThreadOptions(numThreads = 1)
  d <- one_dam(9L, Tm = 25L, n = 200L, pA = 0.3)
  m <- hmm_map(d$dat, phased = d$oph, dam = 1, epsilon = 0.01, paternal_mode = "HWE",
               lambda = 20, tol = 1e-8, maxit = 3000)
  expect_equal(utils::tail(m$fit$loglik_trace, 1),        m$fit$logLik,        tolerance = 1e-9)
  expect_equal(utils::tail(m$fit$penalized_obj_trace, 1), m$fit$penalized_obj, tolerance = 1e-9)
  expect_equal(utils::tail(m$fit$objective_trace, 1),     m$fit$objective,     tolerance = 1e-9)
  # trace length is iters + 1 (initial-through-final)
  expect_length(m$fit$objective_trace, m$fit$iters + 1L)
})

test_that("convergence requires stable r AND stable q", {
  RcppParallel::setThreadOptions(numThreads = 1)
  d <- one_dam(3L, Tm = 25L, n = 200L)
  m <- hmm_map(d$dat, phased = d$oph, dam = 1, epsilon = 0.01, paternal_mode = "HWE",
               lambda = 20, tol = 1e-6, maxit = 3000)
  expect_true(m$fit$converged)
  expect_identical(m$fit$conv_reason, "relative_objective_and_params_stable")
  expect_lt(utils::tail(m$fit$max_dr_trace, 1), 1e-6)   # r stable at convergence
  expect_lt(utils::tail(m$fit$max_dq_trace, 1), 1e-6)   # q stable at convergence
})

test_that("the active objective is monotone (no material decrease)", {
  RcppParallel::setThreadOptions(numThreads = 1)
  d <- one_dam(9001L, Tm = 30L, n = 250L, pA = 0.10)   # extreme pA
  for (lam in c(0, 20)) {
    m <- hmm_map(d$dat, phased = d$oph, dam = 1, epsilon = 0.01, paternal_mode = "HWE",
                 lambda = lam, tol = 1e-9, maxit = 3000)
    expect_true(all(diff(m$fit$objective_trace) >= -1e-7))   # non-decreasing
    expect_false(m$fit$objective_decreased)
  }
})

test_that("low maxit yields non-convergence with a warning and correct metadata", {
  RcppParallel::setThreadOptions(numThreads = 1)
  d <- one_dam(3L, Tm = 25L, n = 200L)
  expect_warning(
    m <- hmm_map(d$dat, phased = d$oph, dam = 1, epsilon = 0.01, paternal_mode = "HWE",
                 tol = 1e-12, maxit = 3),
    "did not converge"
  )
  expect_false(m$fit$converged)
  expect_identical(m$fit$conv_reason, "maxit_reached")
  expect_identical(m$fit$iters, 3L)                          # exactly maxit
  expect_length(m$fit$loglik_trace, 4L)                      # iters + 1
})

test_that("single-dam and D=1 joint agree (logLik, objective, convergence)", {
  RcppParallel::setThreadOptions(numThreads = 1)
  d <- one_dam(11L, Tm = 30L, n = 150L)
  ms <- hmm_map(d$dat, phased = d$oph, dam = 1, epsilon = 0.01, paternal_mode = "HWE",
                lambda = 6, tol = 1e-8, maxit = 3000)
  mj <- hmm_map_joint(d$dat, phased = d$oph, dam = "all", epsilon = 0.01,
                      paternal_mode = "HWE", lambda = 6, tol = 1e-8, maxit = 3000)
  expect_equal(as.numeric(mj$fit$r), as.numeric(ms$fit$r), tolerance = 1e-7)
  expect_equal(mj$fit$logLik,        ms$fit$logLik,        tolerance = 1e-5)
  expect_equal(mj$fit$objective,     ms$fit$objective,     tolerance = 1e-5)
})
