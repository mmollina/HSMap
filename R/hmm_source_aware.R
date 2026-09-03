# ---------------------------------------------------------------------------
# Source-aware multipoint maternal HMM.
#
# The hidden state pairs the transmitted maternal homolog with a paternal
# source state:
#
#     Z_ik = (M_ik, P_ik),   M in {0,1},  P in {U0, U1, H1, H2}
#
# U0/U1 are the ordinary pollen-population modes transmitting allele 0/1;
# H1/H2 are a paternal HAPLOTYPE-SHARING mode in which the transmitted
# paternal allele coincides with the allele carried by the dam's homolog 1 or
# homolog 2. The H states are a statistical description of a shared-haplotype
# tract; they are not a claim of identity by descent.
#
# The maternal process R_k(r_k) is shared by every offspring; only the
# paternal side varies. Setting alpha = 0 makes the H states unreachable and
# recovers the population-mode (four-state) estimator exactly.
# ---------------------------------------------------------------------------

# State order used throughout: s = 4 * M + P, with M in {0,1} and
# P in {U0 = 1, U1 = 2, H1 = 3, H2 = 4}. So states 1:4 carry M = 0 and states
# 5:8 carry M = 1, each block ordered U0, U1, H1, H2.
.hsmap_sa_states <- function() {
  data.frame(
    state = 1:8,
    maternal = rep(0:1, each = 4L),
    paternal = rep(c("U0", "U1", "H1", "H2"), 2L),
    stringsAsFactors = FALSE
  )
}

# Dam homolog allele matrix (2 x T) implied by the fixed phase: row j gives
# the allele count carried by homolog j at each marker.
.hsmap_sa_homolog_alleles <- function(order, phase_vec) {
  h <- .hsmap_marker_phase(order, phase_vec)[order]
  rbind(as.integer(h == 1L), as.integer(h != 1L))
}

# Paternal 4 x 4 transition for one interval.
#   U -> U : (1 - alpha) K_k        population LD kernel
#   U -> H : alpha * (1/2, 1/2)     enter a shared-haplotype tract
#   H -> U : beta  * (1 - q, q)     leave, re-drawing from the local marginal
#   H -> H : (1 - beta) I           stay, retaining homolog identity
.hsmap_sa_paternal_K <- function(K_k, q_next, alpha, beta) {
  M <- matrix(0, 4L, 4L)
  M[1L, 1:2] <- (1 - alpha) * K_k[1L, ]
  M[2L, 1:2] <- (1 - alpha) * K_k[2L, ]
  M[1:2, 3:4] <- alpha * 0.5
  M[3L, 1:2] <- beta * c(1 - q_next, q_next)
  M[4L, 1:2] <- beta * c(1 - q_next, q_next)
  M[3L, 3L] <- 1 - beta
  M[4L, 4L] <- 1 - beta
  M
}

# Joint 8 x 8 transition, T_k = R_k(r_k) (x) Kp_k in the state order above.
.hsmap_sa_transition <- function(r_k, Kp) {
  T_k <- matrix(0, 8L, 8L)
  T_k[1:4, 1:4] <- (1 - r_k) * Kp
  T_k[1:4, 5:8] <- r_k * Kp
  T_k[5:8, 1:4] <- r_k * Kp
  T_k[5:8, 5:8] <- (1 - r_k) * Kp
  T_k
}

# Emissions: list over markers of 8 x n matrices. The latent dosage is
# D = a[M, k] + paternal_allele(P, k) with paternal allele 0, 1, a[1, k],
# a[2, k] for U0, U1, H1, H2. Missing calls emit 1.
.hsmap_sa_emission <- function(Y, A, epsilon) {
  n <- nrow(Y); Tn <- ncol(Y)
  lapply(seq_len(Tn), function(k) {
    y <- Y[, k]
    pat <- c(0L, 1L, A[1L, k], A[2L, k])
    E_k <- matrix(1, 8L, n)
    for (m in 0:1) for (p in 1:4) {
      d <- A[m + 1L, k] + pat[p]
      E_k[4L * m + p, ] <- ifelse(is.na(y), 1,
                                  ifelse(y == d, 1 - epsilon, epsilon / 2))
    }
    E_k
  })
}

# Initial state distribution: maternal uniform; paternal mode at its
# stationary probability pi_H = alpha / (alpha + beta), population allele from
# q_1 and homolog uniform inside the sharing mode.
.hsmap_sa_initial <- function(q1, alpha, beta) {
  piH <- if (alpha + beta > 0) alpha / (alpha + beta) else 0
  p <- c((1 - piH) * (1 - q1), (1 - piH) * q1, piH / 2, piH / 2)
  0.5 * c(p, p)
}

# Scaled forward-backward. Returns the log-likelihood, the pooled expected
# maternal switch counts per interval, the paternal mode counts needed for the
# alpha/beta updates and, optionally, per-offspring posteriors.
.hsmap_sa_forward_backward <- function(Y, A, K_list, q, r, alpha, beta,
                                       epsilon, want_post = FALSE) {
  n <- nrow(Y); Tn <- ncol(Y)
  E <- .hsmap_sa_emission(Y, A, epsilon)
  pi0 <- .hsmap_sa_initial(q[1L], alpha, beta)

  alp <- vector("list", Tn); cs <- matrix(0, Tn, n)
  T_list <- vector("list", Tn - 1L)
  a <- pi0 * E[[1L]]
  cc <- colSums(a); cc[cc <= 0] <- 1e-300
  alp[[1L]] <- a / rep(cc, each = 8L); cs[1L, ] <- cc
  for (k in seq_len(Tn - 1L)) {
    Kp <- .hsmap_sa_paternal_K(K_list[[k]], q[k + 1L], alpha, beta)
    T_k <- .hsmap_sa_transition(r[k], Kp); T_list[[k]] <- T_k
    a <- crossprod(T_k, alp[[k]]) * E[[k + 1L]]
    cc <- colSums(a); cc[cc <= 0] <- 1e-300
    alp[[k + 1L]] <- a / rep(cc, each = 8L); cs[k + 1L, ] <- cc
  }
  ll_i <- colSums(log(cs))

  bet <- matrix(1, 8L, n)
  switches <- numeric(Tn - 1L)
  N <- c(UH = 0, UU = 0, HU = 0, HH = 0)
  gammaM <- gammaH <- NULL; xi <- NULL
  if (want_post) {
    gammaM <- matrix(0, n, Tn); gammaH <- matrix(0, n, Tn)
    xi <- matrix(0, n, Tn - 1L)
    g <- alp[[Tn]] * bet; g <- g / rep(colSums(g), each = 8L)
    gammaM[, Tn] <- colSums(g[5:8, , drop = FALSE])
    gammaH[, Tn] <- colSums(g[c(3L, 4L, 7L, 8L), , drop = FALSE])
  }
  for (k in seq(Tn - 1L, 1L)) {
    B_k <- E[[k + 1L]] * bet
    T_k <- T_list[[k]]
    TB <- T_k %*% B_k
    nrm <- colSums(alp[[k]] * TB); nrm[nrm <= 0] <- 1e-300
    Aw <- alp[[k]] / rep(nrm, each = 8L)
    Xi <- T_k * tcrossprod(Aw, B_k)
    # maternal recombination counts transitions with M' != M only; paternal
    # state changes on their own contribute nothing.
    switches[k] <- sum(Xi[1:4, 5:8]) + sum(Xi[5:8, 1:4])
    for (m1 in 0:1) for (m2 in 0:1) {
      blk <- Xi[4L * m1 + 1:4, 4L * m2 + 1:4]
      N[["UU"]] <- N[["UU"]] + sum(blk[1:2, 1:2])
      N[["UH"]] <- N[["UH"]] + sum(blk[1:2, 3:4])
      N[["HU"]] <- N[["HU"]] + sum(blk[3:4, 1:2])
      N[["HH"]] <- N[["HH"]] + sum(blk[3:4, 3:4])
    }
    if (want_post) {
      s01 <- rep(0, n)
      for (s in 1:4) for (sp in 5:8)
        s01 <- s01 + alp[[k]][s, ] * T_k[s, sp] * B_k[sp, ]
      for (s in 5:8) for (sp in 1:4)
        s01 <- s01 + alp[[k]][s, ] * T_k[s, sp] * B_k[sp, ]
      xi[, k] <- s01 / nrm
    }
    bet <- TB / rep(cs[k + 1L, ], each = 8L)
    if (want_post) {
      g <- alp[[k]] * bet; g <- g / rep(colSums(g), each = 8L)
      gammaM[, k] <- colSums(g[5:8, , drop = FALSE])
      gammaH[, k] <- colSums(g[c(3L, 4L, 7L, 8L), , drop = FALSE])
    }
  }
  # Posterior sharing-mode mass at the first marker. The initial distribution is
  # the stationary law of the U/H process, so alpha and beta enter the expected
  # complete-data log-likelihood through the initial state as well as through
  # the transitions; g1H is the sufficient statistic that term needs. The
  # backward pass has already reached marker 1, so this costs one 8 x n product
  # and never builds the full posterior arrays.
  g1 <- alp[[1L]] * bet
  g1 <- g1 / rep(pmax(colSums(g1), 1e-300), each = 8L)
  g1H <- sum(g1[c(3L, 4L, 7L, 8L), , drop = FALSE])

  list(loglik = sum(ll_i), loglik_i = ll_i, switches = switches, N = N,
       g1H = g1H, gammaM = gammaM, gammaH = gammaH, xi = xi)
}

# Exact M-step for (alpha, beta).
#
# With A = N_UH + g1H and B = N_HU + n - g1H, the terms of the expected
# complete-data log-likelihood that involve alpha and beta are
#
#   Q(a,b) = A log a + N_UU log(1-a) + B log b + N_HH log(1-b) - n log(a+b),
#
# the last term coming from the stationary initial law pi_H = a/(a+b). Setting
# both stationarity equations equal to their common term c = n/(a+b) decouples
# them into quadratics
#
#   c a^2 - (A + N_UU + c) a + A = 0,    c b^2 - (B + N_HH + c) b + B = 0.
#
# For the first, f(0) = A > 0 and f(1) = -N_UU <= 0, so exactly one root lies in
# (0, 1] and it is the smaller one; it is evaluated in the cancellation-free
# form 2P / (S + sqrt(S^2 - 4cP)). The discriminant is non-negative for every
# c >= 0 because its own discriminant in c is -16 A N_UU <= 0.
#
# The M-step is then the scalar root of g(c) = c[a(c) + b(c)] - n, solved in
# log c so that precision is relative. g(0+) = -n < 0 and g(c) -> N_UH + N_HU
# as c -> Inf, so a sign change is guaranteed whenever any U <-> H mass exists.
# Dropping the -n log(a+b) term recovers the transition-only closed forms
# a = N_UH/(N_UH + N_UU) and b = N_HU/(N_HU + N_HH), which are NOT the M-step of
# this model.
.hsmap_sa_mstep_ab <- function(N, g1H, n,
                               a_lo = 1e-8, a_hi = 0.5,
                               b_lo = 1e-8, b_hi = 0.9) {
  A <- N[["UH"]] + g1H
  B <- N[["HU"]] + n - g1H
  UU <- N[["UU"]]; HH <- N[["HH"]]
  clamp <- function(a, b) c(alpha = min(a_hi, max(a_lo, a)),
                            beta  = min(b_hi, max(b_lo, b)))
  # smaller root in (0, 1) of c x^2 - (P + S + c) x + P = 0
  root01 <- function(cc, S, P) {
    if (P <= 0) return(0)
    s <- P + S + cc
    2 * P / (s + sqrt(max(s * s - 4 * cc * P, 0)))
  }
  if (!is.finite(A) || !is.finite(B) || A <= 0)      # no sharing mass at all
    return(clamp(a_lo, N[["HU"]] / max(N[["HU"]] + HH, 1e-9)))

  h <- function(u) {
    cc <- exp(u)
    cc * (root01(cc, UU, A) + root01(cc, HH, B)) - n
  }
  ulo <- log(1e-12); uhi <- log(1e12)
  while (h(uhi) < 0 && uhi < log(1e30)) uhi <- uhi + log(1e6)
  if (h(uhi) < 0) return(clamp(a_lo, b_lo))         # optimum at a + b -> 0
  while (h(ulo) > 0 && ulo > log(1e-30)) ulo <- ulo - log(1e6)
  if (h(ulo) > 0) return(clamp(a_lo, b_lo))
  cc <- exp(stats::uniroot(h, c(ulo, uhi), tol = 1e-12)$root)
  clamp(root01(cc, UU, A), root01(cc, HH, B))
}

# Flank-inferred transmitted paternal alleles, used only to estimate the
# population kernel K. The maternal homolog at marker k is assigned from
# guarded flanking windows, so that the paternal draw at k never uses markers
# whose own paternal alleles are most likely to be dependent with it.
.hsmap_sa_paternal_draws <- function(Y, q, h_al, epsilon = 0.01,
                                     window = 10L, guard = 6L) {
  n <- nrow(Y); Tn <- ncol(Y)
  w <- function(hit) (1 - epsilon) * hit + (epsilon / 2) * (1 - hit)
  LL <- matrix(0, n, Tn)
  for (k in seq_len(Tn)) {
    y <- Y[, k]; obs <- !is.na(y)
    e0 <- q[k] * w(y[obs] == 1) + (1 - q[k]) * w(y[obs] == 0)
    e1 <- q[k] * w(y[obs] == 2) + (1 - q[k]) * w(y[obs] == 1)
    LL[obs, k] <- if (h_al[k]) log(e1) - log(e0) else log(e0) - log(e1)
  }
  CS <- cbind(0, t(apply(LL, 1, cumsum)))
  PAT <- matrix(NA_real_, n, Tn)
  for (k in seq_len(Tn)) {
    lo <- max(1L, k - guard - window); l2 <- k - guard - 1L
    r1 <- k + guard + 1L; hi <- min(Tn, k + guard + window)
    if (l2 - lo < 3L || hi - r1 < 3L) next
    FL <- CS[, l2 + 1L] - CS[, lo]; FR <- CS[, hi + 1L] - CS[, r1]
    keep <- sign(FL) == sign(FR) & abs(FL) > 3 & abs(FR) > 3
    H <- ifelse(FL > 0, 1, 0)
    md <- ifelse(rep(h_al[k], n), H, 1 - H)
    p <- Y[, k] - md
    p[!keep | p < 0 | p > 1] <- NA
    PAT[, k] <- p
  }
  PAT
}

# Population kernel K_k from the paternal draws: smoothed empirical
# conditionals, falling back to independence where support is thin.
.hsmap_sa_population_K <- function(PAT, q, nmin = 60L, smooth = 0.5) {
  Tn <- ncol(PAT)
  K <- vector("list", Tn - 1L)
  npairs <- integer(Tn - 1L); rho <- rep(0, Tn - 1L)
  for (k in seq_len(Tn - 1L)) {
    a <- PAT[, k]; b <- PAT[, k + 1L]
    ok <- !is.na(a) & !is.na(b)
    npairs[k] <- sum(ok)
    if (sum(ok) >= nmin && stats::sd(a[ok]) > 0 && stats::sd(b[ok]) > 0) {
      tb <- table(factor(a[ok], 0:1), factor(b[ok], 0:1)) + smooth
      K[[k]] <- tb / rowSums(tb)
      rho[k] <- suppressWarnings(stats::cor(a[ok], b[ok]))
    } else {
      K[[k]] <- matrix(c(1 - q[k + 1L], q[k + 1L],
                         1 - q[k + 1L], q[k + 1L]), 2L, 2L, byrow = TRUE)
    }
  }
  list(K = K, npairs = npairs, rho = rho)
}

#' Multipoint maternal map with a source-aware paternal process
#'
#' Fits the maternal recombination map of one dam by a hidden Markov model
#' whose state pairs the transmitted maternal homolog with a paternal source
#' state. Alongside the ordinary pollen-population modes the paternal process
#' may occupy a \emph{haplotype-sharing} mode in which the transmitted
#' paternal allele coincides with an allele carried by the dam herself. This
#' matters in open-pollinated families, where some pollen parents are related
#' to the dam: without such a mode the shared-haplotype tracts are absorbed
#' into the maternal process and inflate the estimated map.
#'
#' The maternal recombination fractions are shared by all offspring. The
#' paternal side contributes exactly two estimated scalars, the probability of
#' entering (\code{alpha}) and of leaving (\code{beta}) the sharing mode;
#' within the mode the associated homolog is retained until the mode is left.
#' Setting \code{alpha = 0} makes the sharing states unreachable and recovers
#' the population-mode estimator exactly.
#'
#' The initial state distribution places the paternal process at its stationary
#' sharing probability \eqn{\pi_H = \alpha/(\alpha+\beta)}, so \code{alpha} and
#' \code{beta} enter the expected complete-data log-likelihood through the
#' initial state as well as through the transitions. Their M-step maximizes that
#' full objective, not only its transition part: the two stationarity equations
#' decouple at \eqn{c = n/(\alpha+\beta)} into quadratics with closed-form roots,
#' and the update is the scalar root of \eqn{g(c) = c[\alpha(c)+\beta(c)] - n}.
#' Each EM iteration therefore has the usual ascent guarantee.
#'
#' The haplotype-sharing states are a statistical description of a
#' shared-haplotype tract. They are not evidence of identity by descent.
#'
#' @section Selfed offspring:
#' If an offspring is the product of self-pollination, both gametes come from
#' the same phased diploid dam and the observed dosage
#' \eqn{D_k = a_{M_k,k} + a_{M'_k,k}} is symmetric in the two transmitted
#' paths. Maternal-specific recombination is then not identifiable, whatever
#' estimator is used. Confirmed selfs should therefore be removed by sample
#' quality control before mapping. This is an identifiability statement, not a
#' claim that selfing cannot occur under open pollination; no attempt is made
#' to detect selfs inside the model.
#'
#' @param x An \code{HSMap.data} object.
#' @param phased An \code{HSMap.phased} object giving the marker order and the
#'   dam's phase.
#' @param dam Dam index or name (single dam).
#' @param epsilon Genotyping error rate for the hard-call emission
#'   (default \code{0.01}); not estimated.
#' @param lambda Total pseudocount for the paternal gametic frequency
#'   \code{q} in the population-mode pre-fit (default \code{2}).
#' @param alpha_start,beta_start Starting values for the paternal mode entry
#'   and exit probabilities.
#' @param r_start Starting recombination fraction.
#' @param tol,maxit EM convergence tolerance and iteration cap.
#' @param fit_mode If \code{FALSE}, \code{alpha} and \code{beta} are held at
#'   their starting values instead of being estimated.
#' @param q,K Optional pre-computed paternal marginals and population
#'   kernels; if omitted they are estimated internally.
#' @param nmin,smooth Support threshold and cell smoothing for the population
#'   kernel.
#'
#' @return An object of class \code{HSMap.sourceaware} with elements
#'   \code{r} (interval recombination fractions), \code{pos} (Haldane
#'   positions in cM), \code{alpha}, \code{beta}, \code{loglik},
#'   \code{converged}, \code{iters}, \code{gammaM} (posterior probability that
#'   homolog 1 was transmitted), \code{xi} (per-interval posterior maternal
#'   recombination), \code{gammaH} (posterior occupancy of the sharing mode),
#'   \code{q}, \code{K}, \code{rho} and \code{settings}.
#'
#' @seealso \code{\link{hmm_map}} for the population-mode estimator.
#' @export
hmm_map_source_aware <- function(x, phased, dam = 1,
                                 epsilon = 0.01, lambda = 2,
                                 alpha_start = 0.005, beta_start = 0.02,
                                 r_start = 0.05, tol = 1e-6, maxit = 2000L,
                                 fit_mode = TRUE, q = NULL, K = NULL,
                                 nmin = 60L, smooth = 0.5) {
  if (!inherits(x, "HSMap.data"))
    stop("`x` must be HSMap.data (see read_HSMap_data).", call. = FALSE)
  if (inherits(phased, "HSMap.phased.multi")) phased <- phased[[1L]]
  if (!inherits(phased, "HSMap.phased"))
    stop("`phased` must be an HSMap.phased object.", call. = FALSE)
  if (!is.numeric(epsilon) || length(epsilon) != 1L || epsilon < 0 || epsilon >= 1)
    stop("`epsilon` must be a single number in [0, 1).", call. = FALSE)
  if (alpha_start < 0 || alpha_start >= 1 || beta_start <= 0 || beta_start >= 1)
    stop("`alpha_start` must be in [0, 1) and `beta_start` in (0, 1).",
         call. = FALSE)

  order <- phased$order
  pv <- as.integer(phased$phase_vec)
  A <- .hsmap_sa_homolog_alleles(order, pv)
  Y <- x$G_list[[dam]][, order, drop = FALSE]
  Tn <- length(order)

  # Population-mode pre-fit supplies the paternal marginals q.
  if (is.null(q)) {
    pre <- suppressWarnings(hmm_map(x, phased = phased, dam = dam,
                                    epsilon = epsilon, lambda = lambda,
                                    r_start = r_start, tol = tol,
                                    maxit = maxit))
    q <- as.numeric(pre$fit$q)
  }
  if (length(q) != Tn)
    stop("`q` must have one entry per marker in the order.", call. = FALSE)
  if (is.null(K)) {
    PAT <- .hsmap_sa_paternal_draws(Y, q, A[1L, ] == 1L, epsilon = epsilon)
    pop <- .hsmap_sa_population_K(PAT, q, nmin = nmin, smooth = smooth)
  } else {
    pop <- list(K = K, npairs = rep(NA_integer_, Tn - 1L),
                rho = rep(NA_real_, Tn - 1L))
  }

  alpha <- alpha_start; beta <- beta_start
  r <- rep(r_start, Tn - 1L)
  n <- nrow(Y)
  ll_old <- -Inf; converged <- FALSE; it_done <- 0L
  for (it in seq_len(maxit)) {
    it_done <- it
    fb <- .hsmap_sa_forward_backward(Y, A, pop$K, q, r, alpha, beta, epsilon)
    r_new <- pmin(0.5, pmax(1e-6, fb$switches / n))
    dr <- max(abs(r_new - r)); r <- r_new
    dm <- 0
    if (fit_mode) {
      ab <- .hsmap_sa_mstep_ab(fb$N, fb$g1H, n)
      dm <- max(abs(c(ab[["alpha"]] - alpha, ab[["beta"]] - beta)))
      alpha <- ab[["alpha"]]; beta <- ab[["beta"]]
    }
    if (it > 2L && dr < tol && dm < tol &&
        abs(fb$loglik - ll_old) < tol * (1 + abs(ll_old))) {
      converged <- TRUE; ll_old <- fb$loglik; break
    }
    ll_old <- fb$loglik
  }
  if (!converged)
    warning("source-aware EM did not converge in ", maxit, " iterations.",
            call. = FALSE)

  post <- .hsmap_sa_forward_backward(Y, A, pop$K, q, r, alpha, beta, epsilon,
                                     want_post = TRUE)
  d <- inv_haldane(r)                # package Haldane map function (cM)
  d[r >= 0.499] <- NA_real_
  out <- list(
    r = r, pos = c(0, cumsum(ifelse(is.na(d), 0, d))), dist = d,
    alpha = alpha, beta = beta, loglik = post$loglik,
    loglik_i = post$loglik_i, converged = converged, iters = it_done,
    gammaM = 1 - post$gammaM, xi = post$xi, gammaH = post$gammaH,
    q = q, K = pop$K, rho = pop$rho, npairs = pop$npairs,
    markers = order, phase_vec = pv, dam = dam,
    settings = list(epsilon = epsilon, lambda = lambda, r_start = r_start,
                    alpha_start = alpha_start, beta_start = beta_start,
                    tol = tol, maxit = maxit, fit_mode = fit_mode,
                    nmin = nmin, smooth = smooth, map_function = "haldane")
  )
  class(out) <- "HSMap.sourceaware"
  out
}

#' @export
print.HSMap.sourceaware <- function(x, ...) {
  cat("HSMap source-aware maternal map\n")
  cat(sprintf("  markers   : %d\n", length(x$markers)))
  cat(sprintf("  length    : %.1f cM (Haldane)\n",
              sum(x$dist, na.rm = TRUE)))
  cat(sprintf("  paternal  : alpha = %.5f, beta = %.5f, sharing-mode occupancy = %.3f\n",
              x$alpha, x$beta, mean(x$gammaH)))
  cat(sprintf("  logLik    : %.2f (%d iterations, converged = %s)\n",
              x$loglik, x$iters, x$converged))
  invisible(x)
}
