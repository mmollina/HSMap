# Regression suite for the exact (alpha, beta) M-step of the source-aware HMM.
#
# The initial distribution is the stationary law pi_H = alpha/(alpha+beta), so
# alpha and beta appear in the initial-state term of the expected complete-data
# log-likelihood as well as in the transitions. Maximizing only the transition
# part is the exact M-step iff D = g1H - n*pi_H is exactly zero, and it can
# decrease the observed likelihood. These tests pin the corrected update and
# keep the adversarial cases that exposed the original defect.

# local fixture (test files do not share helpers)
sa_fixture_ms <- function(n = 40L, Tn = 12L, seed = 35L) {
  set.seed(seed)
  h1 <- rbinom(Tn, 1L, 0.5)
  A <- rbind(h1, 1L - h1)
  q <- runif(Tn, 0.2, 0.8)
  K <- lapply(seq_len(Tn - 1L), function(k) {
    m <- matrix(c(0.8, 0.2, 0.3, 0.7), 2L, 2L, byrow = TRUE)
    m / rowSums(m)
  })
  M <- matrix(0L, n, Tn); cur <- rbinom(n, 1L, 0.5); M[, 1L] <- cur
  for (k in seq_len(Tn - 1L)) {
    cur <- ifelse(rbinom(n, 1L, 0.1) == 1L, 1L - cur, cur)
    M[, k + 1L] <- cur
  }
  Fm <- matrix(rbinom(n * Tn, 1L, rep(q, each = n)), n, Tn)
  Y <- matrix(0L, n, Tn)
  for (k in seq_len(Tn)) Y[, k] <- A[cbind(M[, k] + 1L, k)] + Fm[, k]
  Y[1L, 2L] <- NA_integer_
  list(Y = Y, A = A, q = q, K = K, Tn = Tn, n = n)
}

# full Q(alpha, beta), up to a constant free of alpha and beta
sa_Q <- function(a, b, N, g1H, n) {
  N[["UH"]] * log(a) + N[["UU"]] * log1p(-a) +
    N[["HU"]] * log(b) + N[["HH"]] * log1p(-b) +
    g1H * log(a) + (n - g1H) * log(b) - n * log(a + b)
}

sa_gradQ <- function(a, b, N, g1H, n) {
  s <- a + b
  c((N[["UH"]] + g1H) / a - N[["UU"]] / (1 - a) - n / s,
    (N[["HU"]] + n - g1H) / b - N[["HH"]] / (1 - b) - n / s)
}

# independent reference: bounded 2-D optimizer, never used in production
sa_mstep_2d <- function(N, g1H, n, start) {
  fn <- function(p) -sa_Q(p[1L], p[2L], N, g1H, n)
  gr <- function(p) -sa_gradQ(p[1L], p[2L], N, g1H, n)
  stats::optim(start, fn, gr, method = "L-BFGS-B",
               lower = c(1e-8, 1e-8), upper = c(0.5, 0.9),
               control = list(factr = 1e2, pgtol = 1e-12, maxit = 1000L))$par
}

sa_stat_cases <- function() {
  list(
    list(nm = "real-data-like",   N = c(UH = 274, UU = 236800, HU = 250, HH = 12000),
         g1H = 16.9, n = 325),
    list(nm = "high H occupancy", N = c(UH = 900, UU = 9000, HU = 850, HH = 4000),
         g1H = 30, n = 80),
    list(nm = "initial H excess", N = c(UH = 12, UU = 4000, HU = 10, HH = 300),
         g1H = 60, n = 80),
    list(nm = "initial H deficit", N = c(UH = 400, UU = 3000, HU = 380, HH = 2000),
         g1H = 0.2, n = 80),
    list(nm = "short chain",      N = c(UH = 5, UU = 250, HU = 4, HH = 60),
         g1H = 9, n = 30),
    list(nm = "tiny alpha",       N = c(UH = 0.4, UU = 500000, HU = 0.4, HH = 40),
         g1H = 0.3, n = 325),
    list(nm = "large beta",       N = c(UH = 300, UU = 2000, HU = 290, HH = 40),
         g1H = 25, n = 80)
  )
}

test_that("exact M-step zeroes the gradient of the FULL Q", {
  for (cs in sa_stat_cases()) {
    p <- .hsmap_sa_mstep_ab(cs$N, cs$g1H, cs$n)
    g <- sa_gradQ(p[["alpha"]], p[["beta"]], cs$N, cs$g1H, cs$n)
    sc <- max(abs(c((cs$N[["UH"]] + cs$g1H) / p[["alpha"]],
                    cs$N[["UU"]] / (1 - p[["alpha"]]))))
    expect_lt(max(abs(g)) / sc, 1e-9)
  }
})

test_that("the 1-D reduction agrees with an independent 2-D optimizer", {
  for (cs in sa_stat_cases()) {
    p1 <- .hsmap_sa_mstep_ab(cs$N, cs$g1H, cs$n)
    p2 <- sa_mstep_2d(cs$N, cs$g1H, cs$n, c(p1[["alpha"]], p1[["beta"]]))
    expect_equal(unname(p1[["alpha"]]), p2[1L], tolerance = 1e-7)
    expect_equal(unname(p1[["beta"]]), p2[2L], tolerance = 1e-7)
    expect_gte(sa_Q(p1[["alpha"]], p1[["beta"]], cs$N, cs$g1H, cs$n),
               sa_Q(p2[1L], p2[2L], cs$N, cs$g1H, cs$n) - 1e-8)
  }
})

test_that("the scalar root is uniquely bracketed across all cases", {
  for (cs in sa_stat_cases()) {
    A <- cs$N[["UH"]] + cs$g1H
    B <- cs$N[["HU"]] + cs$n - cs$g1H
    root01 <- function(cc, S, P) {
      if (P <= 0) return(0)
      s <- P + S + cc
      2 * P / (s + sqrt(max(s * s - 4 * cc * P, 0)))
    }
    g <- function(cc) {
      cc * (root01(cc, cs$N[["UU"]], A) + root01(cc, cs$N[["HH"]], B)) - cs$n
    }
    gv <- vapply(exp(seq(log(1e-8), log(1e14), length.out = 500)), g, 0)
    expect_true(all(diff(gv) > 0))
    expect_equal(sum(diff(sign(gv)) != 0), 1L)
    expect_lt(gv[1L], 0)
    expect_gt(gv[length(gv)], 0)
  }
})

test_that("the transition-only update is NOT the M-step unless D = 0", {
  cs <- sa_stat_cases()[[3L]]
  old_a <- cs$N[["UH"]] / (cs$N[["UH"]] + cs$N[["UU"]])
  old_b <- cs$N[["HU"]] / (cs$N[["HU"]] + cs$N[["HH"]])
  D <- cs$g1H - cs$n * old_a / (old_a + old_b)
  expect_gt(abs(D), 1)
  g <- sa_gradQ(old_a, old_b, cs$N, cs$g1H, cs$n)
  expect_equal(g[1L], D / old_a, tolerance = 1e-8)
  expect_equal(g[2L], -D / old_b, tolerance = 1e-8)
  p <- .hsmap_sa_mstep_ab(cs$N, cs$g1H, cs$n)
  expect_gt(sa_Q(p[["alpha"]], p[["beta"]], cs$N, cs$g1H, cs$n),
            sa_Q(old_a, old_b, cs$N, cs$g1H, cs$n))
})

# simulate from the implemented model so the H tracts are genuine; force_H1
# starts every offspring in the sharing mode, the configuration that made the
# pre-correction update decrease the observed likelihood
sa_sim <- function(n, Tn, r, alpha, beta, seed, force_H1 = FALSE) {
  set.seed(seed)
  h <- rep(c(1L, 0L), length.out = Tn)
  A <- rbind(h, 1L - h)
  q <- runif(Tn, 0.25, 0.75)
  piH <- alpha / (alpha + beta)
  Y <- matrix(0L, n, Tn)
  for (i in seq_len(n)) {
    M <- integer(Tn); P <- integer(Tn)
    M[1L] <- rbinom(1L, 1L, 0.5)
    P[1L] <- if (force_H1) 3L else
      sample.int(4L, 1L, prob = c((1 - piH) * (1 - q[1L]), (1 - piH) * q[1L],
                                  piH / 2, piH / 2))
    for (k in seq_len(Tn - 1L)) {
      M[k + 1L] <- if (runif(1) < r[k]) 1L - M[k] else M[k]
      P[k + 1L] <- if (P[k] <= 2L) {
        if (runif(1) < alpha) sample(3:4, 1L) else
          sample.int(2L, 1L, prob = c(1 - q[k + 1L], q[k + 1L]))
      } else if (runif(1) < beta) {
        sample.int(2L, 1L, prob = c(1 - q[k + 1L], q[k + 1L]))
      } else P[k]
    }
    Y[i, ] <- vapply(seq_len(Tn), function(k)
      A[M[k] + 1L, k] + c(0L, 1L, A[1L, k], A[2L, k])[P[k]], 0L)
  }
  K <- lapply(seq_len(Tn - 1L), function(k)
    rbind(c(1 - q[k + 1L], q[k + 1L]), c(1 - q[k + 1L], q[k + 1L])))
  list(Y = Y, A = A, q = q, K = K, n = n, Tn = Tn)
}

test_that("EM is monotone on the adversarial cases that broke the old update", {
  cases <- list(list(Tn = 25L, n = 80L, a = 0.002, b = 0.10, fH = TRUE,  s = 1L),
                list(Tn = 60L, n = 80L, a = 0.002, b = 0.01, fH = TRUE,  s = 2L),
                list(Tn = 25L, n = 30L, a = 0.002, b = 0.10, fH = TRUE,  s = 3L),
                list(Tn = 10L, n = 30L, a = 0.050, b = 0.01, fH = FALSE, s = 4L),
                list(Tn = 15L, n = 80L, a = 0.002, b = 0.10, fH = TRUE,  s = 5L))
  for (cc in cases) {
    f <- sa_sim(cc$n, cc$Tn, rep(0.06, cc$Tn - 1L), cc$a, cc$b,
                seed = 4000L + cc$s, force_H1 = cc$fH)
    for (st in list(c(0.005, 0.02), c(0.2, 0.6))) {
      r <- rep(0.05, f$Tn - 1L); a <- st[1L]; b <- st[2L]
      lls <- numeric(0)
      for (it in 1:30) {
        fb <- .hsmap_sa_forward_backward(f$Y, f$A, f$K, f$q, r, a, b, 0.01)
        lls <- c(lls, fb$loglik)
        r <- pmin(0.5, pmax(1e-6, fb$switches / f$n))
        p <- .hsmap_sa_mstep_ab(fb$N, fb$g1H, f$n)
        expect_gte(sa_Q(p[["alpha"]], p[["beta"]], fb$N, fb$g1H, f$n),
                   sa_Q(a, b, fb$N, fb$g1H, f$n) - 1e-10)
        a <- p[["alpha"]]; b <- p[["beta"]]
      }
      expect_true(all(diff(lls) >= 0))
    }
  }
})

test_that("EM converges to the self-consistent stationary solution", {
  f <- sa_sim(120L, 40L, rep(0.05, 39L), 0.06, 0.15, seed = 9001L)
  r <- rep(0.05, 39L); a <- 0.005; b <- 0.02
  for (it in 1:2000) {
    fb <- .hsmap_sa_forward_backward(f$Y, f$A, f$K, f$q, r, a, b, 0.01)
    r_new <- pmin(0.5, pmax(1e-6, fb$switches / f$n))
    p <- .hsmap_sa_mstep_ab(fb$N, fb$g1H, f$n)
    done <- max(abs(r_new - r)) < 1e-12 &&
      max(abs(c(p[["alpha"]] - a, p[["beta"]] - b))) < 1e-14
    r <- r_new; a <- p[["alpha"]]; b <- p[["beta"]]
    if (done) break
  }
  # the fixed point must satisfy BOTH conditions with statistics taken at the
  # same parameters: applying the M-step reproduces (alpha, beta), and the full
  # Q gradient vanishes there.
  fb <- .hsmap_sa_forward_backward(f$Y, f$A, f$K, f$q, r, a, b, 0.01)
  p <- .hsmap_sa_mstep_ab(fb$N, fb$g1H, f$n)
  expect_equal(unname(p[["alpha"]]), a, tolerance = 1e-10)
  expect_equal(unname(p[["beta"]]), b, tolerance = 1e-10)
  g <- sa_gradQ(p[["alpha"]], p[["beta"]], fb$N, fb$g1H, f$n)
  sc <- (fb$N[["UH"]] + fb$g1H) / p[["alpha"]]
  expect_lt(max(abs(g)) / sc, 1e-9)
  expect_true(a > 1e-8 && a < 0.5 && b > 1e-8 && b < 0.9)
})

test_that("g1H is a valid mass and needs no posterior arrays", {
  f <- sa_fixture_ms(n = 50L, Tn = 18L, seed = 21L)
  r <- rep(0.08, f$Tn - 1L)
  lean <- .hsmap_sa_forward_backward(f$Y, f$A, f$K, f$q, r, 0.05, 0.1, 0.01)
  full <- .hsmap_sa_forward_backward(f$Y, f$A, f$K, f$q, r, 0.05, 0.1, 0.01,
                                     want_post = TRUE)
  expect_true(lean$g1H >= 0 && lean$g1H <= f$n)
  expect_null(lean$gammaH)
  expect_equal(lean$g1H, sum(full$gammaH[, 1L]), tolerance = 1e-10)
})

test_that("alpha held at 0 still reproduces the population-mode model", {
  f <- sa_fixture_ms(n = 50L, Tn = 14L, seed = 11L)
  r <- runif(f$Tn - 1L, 0.01, 0.2)
  sa <- .hsmap_sa_forward_backward(f$Y, f$A, f$K, f$q, r, 0, 0.02, 0.01,
                                   want_post = TRUE)
  expect_equal(max(sa$gammaH), 0)
  expect_equal(sa$g1H, 0)
  p <- .hsmap_sa_mstep_ab(sa$N, sa$g1H, f$n)
  expect_equal(unname(p[["alpha"]]), 1e-8)
})
