# Tests for the source-aware multipoint estimator (8-state (M, P) HMM).
# The statistical core is tested directly; the alpha = 0 nesting test is a
# permanent regression guard.

# --- fixtures ---------------------------------------------------------------
sa_fixture <- function(n = 40L, Tn = 12L, seed = 35L) {
  set.seed(seed)
  h1 <- rbinom(Tn, 1L, 0.5)
  A <- rbind(h1, 1L - h1)                       # dam homolog alleles
  q <- runif(Tn, 0.2, 0.8)
  K <- lapply(seq_len(Tn - 1L), function(k) {
    m <- matrix(c(0.8, 0.2, 0.3, 0.7), 2L, 2L, byrow = TRUE); m / rowSums(m)
  })
  M <- matrix(0L, n, Tn); cur <- rbinom(n, 1L, 0.5); M[, 1L] <- cur
  for (k in seq_len(Tn - 1L)) {
    cur <- ifelse(rbinom(n, 1L, 0.1) == 1L, 1L - cur, cur); M[, k + 1L] <- cur
  }
  Fm <- matrix(rbinom(n * Tn, 1L, rep(q, each = n)), n, Tn)
  Y <- matrix(0L, n, Tn)
  for (k in seq_len(Tn)) Y[, k] <- A[cbind(M[, k] + 1L, k)] + Fm[, k]
  Y[1L, 2L] <- NA_integer_                       # exercise missing handling
  colnames(Y) <- sprintf("m%02d", seq_len(Tn))
  rownames(Y) <- sprintf("i%02d", seq_len(n))
  list(Y = Y, A = A, q = q, K = K, Tn = Tn, n = n)
}

# --- A. paternal transition matrices ----------------------------------------
test_that("paternal transition is a valid 4x4 stochastic matrix", {
  f <- sa_fixture()
  for (a in c(0, 0.001, 0.05, 0.5)) for (b in c(1e-8, 0.02, 0.5, 0.9)) {
    Kp <- .hsmap_sa_paternal_K(f$K[[1L]], f$q[2L], a, b)
    expect_equal(dim(Kp), c(4L, 4L))
    expect_true(all(Kp >= 0))
    expect_equal(rowSums(Kp), rep(1, 4L), tolerance = 1e-12)
  }
})

test_that("alpha = 0 makes the sharing states unreachable from U", {
  f <- sa_fixture()
  Kp <- .hsmap_sa_paternal_K(f$K[[1L]], f$q[2L], 0, 0.02)
  expect_equal(Kp[1:2, 3:4], matrix(0, 2L, 2L))
  expect_equal(Kp[1:2, 1:2], f$K[[1L]], tolerance = 1e-12)
})

test_that("beta = 0 makes the sharing mode absorbing and homolog-preserving", {
  f <- sa_fixture()
  Kp <- .hsmap_sa_paternal_K(f$K[[1L]], f$q[2L], 0.01, 0)
  expect_equal(Kp[3:4, 1:2], matrix(0, 2L, 2L))
  expect_equal(Kp[3:4, 3:4], diag(2L))          # H1 -> H1, H2 -> H2
})

test_that("exit distribution follows the local paternal marginal", {
  f <- sa_fixture()
  b <- 0.3; qn <- f$q[2L]
  Kp <- .hsmap_sa_paternal_K(f$K[[1L]], qn, 0.01, b)
  expect_equal(Kp[3L, 1:2], b * c(1 - qn, qn), tolerance = 1e-12)
  expect_equal(Kp[4L, 1:2], b * c(1 - qn, qn), tolerance = 1e-12)
})

# --- B. joint transitions ----------------------------------------------------
test_that("joint transition is 8x8, stochastic, with the right maternal marginal", {
  f <- sa_fixture()
  Kp <- .hsmap_sa_paternal_K(f$K[[1L]], f$q[2L], 0.01, 0.02)
  r <- 0.13
  T_k <- .hsmap_sa_transition(r, Kp)
  expect_equal(dim(T_k), c(8L, 8L))
  expect_true(all(T_k >= 0))
  expect_equal(rowSums(T_k), rep(1, 8L), tolerance = 1e-12)
  # marginalizing the paternal component must return R_k(r)
  mat <- matrix(0, 2L, 2L)
  for (m1 in 0:1) for (m2 in 0:1)
    mat[m1 + 1L, m2 + 1L] <- sum(T_k[4L * m1 + 1L, 4L * m2 + 1:4])
  expect_equal(mat, matrix(c(1 - r, r, r, 1 - r), 2L, 2L), tolerance = 1e-12)
})

# --- C. emissions ------------------------------------------------------------
test_that("emissions encode the correct latent dosage for every state", {
  f <- sa_fixture(); eps <- 0.01
  E <- .hsmap_sa_emission(f$Y, f$A, eps)
  k <- 3L
  pat <- c(0L, 1L, f$A[1L, k], f$A[2L, k])
  for (m in 0:1) for (p in 1:4) {
    d <- f$A[m + 1L, k] + pat[p]
    got <- unname(E[[k]][4L * m + p, ])
    want <- unname(ifelse(is.na(f$Y[, k]), 1,
                          ifelse(f$Y[, k] == d, 1 - eps, eps / 2)))
    expect_equal(got, want, tolerance = 1e-12)
  }
})

test_that("missing genotypes emit 1 for every state", {
  f <- sa_fixture(); E <- .hsmap_sa_emission(f$Y, f$A, 0.01)
  expect_equal(unname(E[[2L]][, 1L]), rep(1, 8L))  # i01 is NA at marker 2
})

test_that("sharing states emit the phased dam homolog allele", {
  f <- sa_fixture(); E <- .hsmap_sa_emission(f$Y, f$A, 0.01)
  k <- 5L
  # H1 (index 3 within the M = 0 block) must behave like paternal allele a_1k
  a1 <- f$A[1L, k]
  u_row <- if (a1 == 0L) 1L else 2L                 # matching U state
  expect_equal(E[[k]][3L, ], E[[k]][u_row, ], tolerance = 1e-12)
})

# --- D. forward-backward -----------------------------------------------------
test_that("forward-backward returns finite likelihood and valid posteriors", {
  f <- sa_fixture()
  r <- rep(0.1, f$Tn - 1L)
  fb <- .hsmap_sa_forward_backward(f$Y, f$A, f$K, f$q, r, 0.01, 0.02, 0.01,
                                   want_post = TRUE)
  expect_true(is.finite(fb$loglik))
  expect_equal(length(fb$loglik_i), f$n)
  expect_true(all(fb$gammaM >= 0 & fb$gammaM <= 1))
  expect_true(all(fb$gammaH >= 0 & fb$gammaH <= 1))
  expect_true(all(fb$xi >= 0 & fb$xi <= 1))
  expect_equal(fb$loglik, sum(fb$loglik_i), tolerance = 1e-9)
})

# --- E. maternal xi ----------------------------------------------------------
test_that("paternal-only transitions contribute no maternal recombination", {
  f <- sa_fixture()
  # r = 0 forbids maternal switching structurally; xi must be exactly 0 even
  # though the paternal process is free to move between U and H states.
  fb <- .hsmap_sa_forward_backward(f$Y, f$A, f$K, f$q, rep(0, f$Tn - 1L),
                                   0.05, 0.05, 0.01, want_post = TRUE)
  expect_equal(max(fb$xi), 0)
  expect_equal(max(fb$switches), 0)
})

# --- F. EM behaviour ---------------------------------------------------------
test_that("EM increases the likelihood and respects parameter bounds", {
  f <- sa_fixture(n = 60L, Tn = 15L, seed = 7L)
  lls <- numeric(0); r <- rep(0.05, f$Tn - 1L); a <- 0.005; b <- 0.02
  for (it in 1:8) {
    fb <- .hsmap_sa_forward_backward(f$Y, f$A, f$K, f$q, r, a, b, 0.01)
    lls <- c(lls, fb$loglik)
    r <- pmin(0.5, pmax(1e-6, fb$switches / f$n))
    a <- max(1e-8, min(0.5, fb$N[["UH"]] / (fb$N[["UH"]] + fb$N[["UU"]])))
    b <- max(1e-8, min(0.9, fb$N[["HU"]] / (fb$N[["HU"]] + fb$N[["HH"]])))
  }
  expect_true(all(diff(lls) > -1e-6))            # nondecreasing to tolerance
  expect_true(all(r >= 1e-6 & r <= 0.5))
  expect_true(a >= 0 && a < 1 && b > 0 && b < 1)
})

# --- G. exact nesting on alpha = 0 ------------------------------------------
test_that("alpha exactly 0 reproduces the population-mode model exactly", {
  f <- sa_fixture(n = 50L, Tn = 14L, seed = 11L)
  r <- runif(f$Tn - 1L, 0.01, 0.2)
  sa <- .hsmap_sa_forward_backward(f$Y, f$A, f$K, f$q, r, 0, 0.02, 0.01,
                                   want_post = TRUE)
  # reference: the same chain restricted to the four population states
  ref_ll <- local({
    n <- f$n; Tn <- f$Tn
    E <- .hsmap_sa_emission(f$Y, f$A, 0.01)[seq_len(Tn)]
    E4 <- lapply(E, function(e) e[c(1L, 2L, 5L, 6L), , drop = FALSE])
    pi0 <- 0.5 * c(1 - f$q[1L], f$q[1L], 1 - f$q[1L], f$q[1L])
    a <- pi0 * E4[[1L]]; cc <- colSums(a); ll <- sum(log(cc))
    al <- a / rep(cc, each = 4L)
    for (k in seq_len(Tn - 1L)) {
      Kk <- f$K[[k]]; rr <- r[k]
      T4 <- matrix(0, 4L, 4L)
      T4[1:2, 1:2] <- (1 - rr) * Kk; T4[1:2, 3:4] <- rr * Kk
      T4[3:4, 1:2] <- rr * Kk;       T4[3:4, 3:4] <- (1 - rr) * Kk
      a <- crossprod(T4, al) * E4[[k + 1L]]
      cc <- colSums(a); ll <- ll + sum(log(cc)); al <- a / rep(cc, each = 4L)
    }
    ll
  })
  expect_equal(sa$loglik, ref_ll, tolerance = 1e-10)
  expect_equal(max(sa$gammaH), 0)                # sharing states unreachable
})

test_that("alpha = 1e-12 is NOT a valid nesting test", {
  # Documented rationale: with a strictly positive entry probability the
  # sharing states remain reachable, and on real data the likelihood strongly
  # prefers them, so the tiny prior is overwhelmed. Only alpha exactly 0
  # removes the states from the chain.
  f <- sa_fixture()
  z <- .hsmap_sa_paternal_K(f$K[[1L]], f$q[2L], 1e-12, 0.02)
  expect_gt(sum(z[1:2, 3:4]), 0)                 # still reachable
  z0 <- .hsmap_sa_paternal_K(f$K[[1L]], f$q[2L], 0, 0.02)
  expect_equal(sum(z0[1:2, 3:4]), 0)             # only alpha = 0 removes them
})

# --- helpers -----------------------------------------------------------------
test_that("state table documents the state order used by the engine", {
  s <- .hsmap_sa_states()
  expect_equal(nrow(s), 8L)
  expect_equal(s$maternal, rep(0:1, each = 4L))
  expect_equal(s$paternal[1:4], c("U0", "U1", "H1", "H2"))
})
