# Marker-specific emission-error model (epsilon_model = "marker_specific").
#
# Validation contract (m24-marker-specific-error-hmm, section 4):
#  A. global mode numerically unchanged (matches the independent R likelihood)
#  B. fixed equal epsilon_k likelihood == global likelihood
#  C. epsilon_k MAP update matches hand-computed expected error counts
#  D. missing genotypes contribute to neither N_k nor E_k
#  E. penalized objective non-decreasing under EM
#  F. deterministic under fixed inputs
#  G. no NaN / boundary failures

LAM <- 2

## independent R likelihood with PER-MARKER epsilon (extends the m20-validated
## global evaluator; gametic emission at dam-het markers, phase handled)
ll_ref <- function(Y, phase_vec, r, q, eps) {
  n <- nrow(Y); Tn <- ncol(Y)
  if (length(eps) == 1L) eps <- rep(eps, Tn)
  E0 <- matrix(1, n, Tn); E1 <- matrix(1, n, Tn)
  for (k in seq_len(Tn)) {
    y <- Y[, k]; obs <- !is.na(y)
    w <- function(hit) (1 - eps[k]) * hit + (eps[k] / 2) * (1 - hit)
    E0[obs, k] <- q[k] * w(y[obs] == 1) + (1 - q[k]) * w(y[obs] == 0)
    E1[obs, k] <- q[k] * w(y[obs] == 2) + (1 - q[k]) * w(y[obs] == 1)
  }
  a0 <- 0.5 * E0[, 1]; a1 <- 0.5 * E1[, 1]
  cc <- a0 + a1; ll <- log(cc); a0 <- a0 / cc; a1 <- a1 / cc
  for (k in seq_len(Tn - 1L)) {
    ps <- if (phase_vec[k] == 1L) 1 - r[k] else r[k]
    b0 <- (a0 * ps + a1 * (1 - ps)) * E0[, k + 1L]
    b1 <- (a1 * ps + a0 * (1 - ps)) * E1[, k + 1L]
    cc <- b0 + b1; cc[cc <= 0] <- 1e-300
    ll <- ll + log(cc); a0 <- b0 / cc; a1 <- b1 / cc
  }
  sum(ll)
}

## hand E-step for the error indicator at given params (all dam-het markers)
ek_ref <- function(Y, phase_vec, r, q, eps) {
  n <- nrow(Y); Tn <- ncol(Y)
  if (length(eps) == 1L) eps <- rep(eps, Tn)
  P0 <- matrix(1, n, Tn); EM0 <- matrix(1, n, Tn); EM1 <- matrix(1, n, Tn)
  Pc0 <- matrix(1, n, Tn); Pc1 <- matrix(1, n, Tn)
  for (k in seq_len(Tn)) {
    y <- Y[, k]; obs <- !is.na(y)
    Pc0[obs, k] <- q[k] * (y[obs] == 1) + (1 - q[k]) * (y[obs] == 0)  # clean
    Pc1[obs, k] <- q[k] * (y[obs] == 2) + (1 - q[k]) * (y[obs] == 1)
    EM0[obs, k] <- (1 - eps[k]) * Pc0[obs, k] + (eps[k] / 2) * (1 - Pc0[obs, k])
    EM1[obs, k] <- (1 - eps[k]) * Pc1[obs, k] + (eps[k] / 2) * (1 - Pc1[obs, k])
  }
  ## forward-backward state posteriors
  A0 <- matrix(0, n, Tn); A1 <- matrix(0, n, Tn)
  a0 <- 0.5 * EM0[, 1]; a1 <- 0.5 * EM1[, 1]; cc <- a0 + a1
  A0[, 1] <- a0 / cc; A1[, 1] <- a1 / cc
  for (k in seq_len(Tn - 1L)) {
    ps <- if (phase_vec[k] == 1L) 1 - r[k] else r[k]
    b0 <- (A0[, k] * ps + A1[, k] * (1 - ps)) * EM0[, k + 1L]
    b1 <- (A1[, k] * ps + A0[, k] * (1 - ps)) * EM1[, k + 1L]
    cc <- b0 + b1; A0[, k + 1L] <- b0 / cc; A1[, k + 1L] <- b1 / cc
  }
  B0 <- matrix(1, n, Tn); B1 <- matrix(1, n, Tn)
  for (k in seq(Tn - 1L, 1L)) {
    ps <- if (phase_vec[k] == 1L) 1 - r[k] else r[k]
    u0 <- EM0[, k + 1L] * B0[, k + 1L]; u1 <- EM1[, k + 1L] * B1[, k + 1L]
    v0 <- ps * u0 + (1 - ps) * u1; v1 <- ps * u1 + (1 - ps) * u0
    cc <- v0 + v1; B0[, k] <- v0 / cc; B1[, k] <- v1 / cc
  }
  g0 <- A0 * B0; g1 <- A1 * B1; gs <- g0 + g1
  g0 <- g0 / gs; g1 <- g1 / gs
  perr <- g0 * (eps[col(Y)] / 2) * (1 - Pc0) / EM0 +
          g1 * (eps[col(Y)] / 2) * (1 - Pc1) / EM1
  perr[is.na(Y)] <- 0
  colSums(perr)
}

sim_eps <- function(seed, Tm = 30L, n = 300L, r = 0.08, err = 0.01, miss = 0) {
  set.seed(seed)
  sim <- sim_multi_pop(T_markers = Tm, n_pops = 1, n_ind_per_pop = n,
                       marker_intersection = 1, r_vec = rep(r, Tm - 1L),
                       phase_mode = "random", repulsion_rate = 0.3,
                       maternal_geno_mode = "all_het", maternal_pA = 0.5,
                       paternal_pA_base = 0.4, error_rate = err, seed = seed,
                       miss_rate = miss)
  mk <- sim$truth$markers_union
  list(dat = make_dat(sim),
       oph = oracle_phased(mk, 1L - sim$truth$v_true[[1]], "P1"), mk = mk)
}
fit_mode <- function(s, model, eps = 0.01, tau = 100, maxit = 300L, ...) {
  suppressWarnings(hmm_map(s$dat, phased = s$oph, dam = 1, epsilon = eps,
                           epsilon_model = model, epsilon_tau = tau,
                           lambda = LAM, tol = 1e-6, maxit = maxit, ...))
}
## orient Y/phase exactly as the engine sees them
eng_view <- function(s) {
  Y <- s$dat$G_list[[1]][, s$mk, drop = FALSE]
  list(Y = Y, pv = as.integer(s$oph$phase_vec))
}

test_that("A/B: global mode matches the independent likelihood; equal epsilon_k == global", {
  s <- sim_eps(11L)
  v <- eng_view(s)
  mg <- fit_mode(s, "global")
  expect_null(mg$fit$epsilon_k)
  ## A: engine logLik at final params equals the R reference (global scalar)
  llA <- ll_ref(v$Y, v$pv, as.numeric(mg$fit$r), as.numeric(mg$fit$q), 0.01)
  expect_lt(abs(llA - as.numeric(mg$fit[["logLik"]])), 1e-6)
  ## B: the reference with a CONSTANT eps vector equals the same value
  llB <- ll_ref(v$Y, v$pv, as.numeric(mg$fit$r), as.numeric(mg$fit$q),
                rep(0.01, length(s$mk)))
  expect_identical(llA, llB)
  ## and a marker_specific fit's likelihood matches the reference at ITS eps_k
  ms <- fit_mode(s, "marker_specific")
  llC <- ll_ref(v$Y, v$pv, as.numeric(ms$fit$r), as.numeric(ms$fit$q),
                as.numeric(ms$fit$epsilon_k))
  expect_lt(abs(llC - as.numeric(ms$fit[["logLik"]])), 1e-6)
})

test_that("C: epsilon_k MAP update matches hand-computed expected error counts", {
  s <- sim_eps(12L, Tm = 12L, n = 150L)
  v <- eng_view(s)
  tau <- 50; eps0 <- 0.02
  ## one EM iteration from known initial params: r = r_start, q = 0.5 (flat prior)
  m1 <- fit_mode(s, "marker_specific", eps = eps0, tau = tau, maxit = 1L,
                 r_start = 0.07)
  Ek_hand <- ek_ref(v$Y, v$pv, rep(0.07, length(s$mk) - 1L),
                    rep(0.5, length(s$mk)), eps0)
  Nk <- colSums(!is.na(v$Y))
  eps_hand <- (Ek_hand + tau * eps0) / (Nk + tau)
  expect_lt(max(abs(as.numeric(m1$fit$epsilon_k) - eps_hand)), 1e-8)
  expect_identical(as.integer(m1$fit$N_k), as.integer(Nk))
})

test_that("D: missing genotypes contribute to neither N_k nor E_k", {
  s <- sim_eps(13L, Tm = 12L, n = 120L, miss = 0.15)
  ## force one marker fully missing
  s$dat$G_list[[1]][, s$mk[5L]] <- NA_integer_
  ms <- fit_mode(s, "marker_specific", eps = 0.01, tau = 80)
  expect_identical(as.integer(ms$fit$N_k),
                   as.integer(colSums(!is.na(s$dat$G_list[[1]][, s$mk]))))
  ## an all-missing marker has E_k = 0 and epsilon_k shrunk exactly to eps0
  expect_equal(unname(ms$fit$E_k[5L]), 0)
  expect_equal(unname(ms$fit$epsilon_k[5L]), 0.01, tolerance = 1e-12)
})

test_that("E: penalized objective is non-decreasing under EM", {
  s <- sim_eps(14L)
  ms <- fit_mode(s, "marker_specific")
  tr <- as.numeric(ms$fit$objective_trace)
  expect_false(isTRUE(ms$fit$objective_decreased))
  expect_true(all(diff(tr) > -1e-6 * (1 + abs(tr[-length(tr)]))))
  ## penalized_obj is reported (not NA) in marker_specific mode
  expect_false(is.na(ms$fit$penalized_obj))
})

test_that("F/G: deterministic; no NaN or boundary failures", {
  s <- sim_eps(15L, miss = 0.05)
  a <- fit_mode(s, "marker_specific")
  b <- fit_mode(s, "marker_specific")
  expect_identical(as.numeric(a$fit$epsilon_k), as.numeric(b$fit$epsilon_k))
  expect_identical(as.numeric(a$fit$r), as.numeric(b$fit$r))
  ek <- as.numeric(a$fit$epsilon_k)
  expect_true(all(is.finite(ek)))
  expect_true(all(ek > 0 & ek < 0.5))
})

test_that("a planted high-error marker earns a higher epsilon_k", {
  s <- sim_eps(16L, Tm = 30L, n = 300L, err = 0.002)
  bad <- s$mk[15L]
  G <- s$dat$G_list[[1]]
  set.seed(99)
  y <- G[, bad]; obs <- which(!is.na(y))
  hit <- obs[runif(length(obs)) < 0.06]
  y[hit] <- vapply(y[hit], function(g) sample(setdiff(0:2, g), 1L), 0L)
  G[, bad] <- y; s$dat$G_list[[1]] <- G
  ms <- fit_mode(s, "marker_specific", eps = 0.01, tau = 50)
  ek <- ms$fit$epsilon_k
  expect_gt(unname(ek[bad]), 2 * median(ek[setdiff(s$mk, bad)]))
})

test_that("input validation and joint guard", {
  s <- sim_eps(17L, Tm = 10L, n = 80L)
  expect_error(fit_mode(s, "marker_specific", tau = -1), "positive")
  expect_error(suppressWarnings(hmm_map(s$dat, phased = s$oph, dam = 1,
               epsilon_model = "marker_specific", paternal_mode = "two_locus",
               lambda = LAM)), "gametic")
})
