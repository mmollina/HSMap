#' Test for heterogeneity of the recombination map across dams
#'
#' @description
#' After fitting a joint shared map (\code{\link{hmm_map_joint}} or
#' \code{hmm_map(method = "joint")}), assess whether individual dams depart from the
#' shared map by per-dam scalings of Haldane distances. For dam \eqn{d} the scaled
#' alternative sets \eqn{m_k^{(d)} = \eta^{(d)} m_k}, where \eqn{m_k} is the Haldane
#' distance of the shared map, so that
#' \deqn{r_k^{(d)} = \tfrac12\bigl\{1 - \exp(-2\,\eta^{(d)} m_k)\bigr\}.}
#' Each \eqn{\eta^{(d)}} is profiled by a one-dimensional search conditional on the
#' shared \eqn{\mathbf r}. \eqn{\eta^{(d)} = 1} for all dams is the homogeneous null.
#'
#' @details
#' The test is a nested likelihood-ratio test: the null allows a single common
#' scale \eqn{\eta} (so the overall map length is free but shared), the alternative
#' allows a separate \eqn{\eta^{(d)}} per dam, giving
#' \deqn{\mathrm{LR} = 2\Bigl[\textstyle\sum_d \ell_d(\hat\eta^{(d)}) -
#'       \textstyle\sum_d \ell_d(\hat\eta_{\mathrm{common}})\Bigr]}
#' on \eqn{D-1} degrees of freedom. (When \eqn{\mathbf r} is the joint MLE,
#' \eqn{\hat\eta_{\mathrm{common}} \approx 1}; estimating it explicitly removes the
#' overall-scale degree of freedom, the manuscript's "one dam fixed for
#' identifiability".) The chi-square calibration is asymptotic; for small \eqn{D} or
#' a definitive p-value, calibrate by parametric bootstrap. The per-dam
#' \eqn{\hat\eta^{(d)}} are the primary diagnostic: \eqn{\hat\eta^{(d)} > 1} means
#' dam \eqn{d} prefers a longer map than the shared one.
#'
#' @param dat An \code{HSMap.data} object.
#' @param map A joint shared map of class \code{"HSMap.map.joint"}.
#' @param epsilon Optional override for the emission error rate; defaults to the
#'   value stored in the fit.
#' @param eta_range Length-2 positive, increasing search interval for \eqn{\eta}.
#'   Default \code{c(0.1, 10)}.
#'
#' @return An object of class \code{"HSMap.hetero"}: a list with \code{eta} (named
#'   \eqn{\hat\eta^{(d)}}), \code{eta_common}, \code{loglik_null}, \code{loglik_alt},
#'   \code{per_dam} (data frame), and \code{LR}, \code{df}, \code{p_value}.
#'
#' @seealso \code{\link{hmm_map_joint}}
#' @export
test_map_heterogeneity <- function(dat, map, epsilon = NULL, eta_range = c(0.1, 10)) {
  `%||%` <- function(a, b) if (is.null(a)) b else a
  if (!inherits(dat, "HSMap.data")) stop("`dat` must be an HSMap.data object.")
  if (!inherits(map, "HSMap.map.joint"))
    stop("`map` must be a joint shared map (class 'HSMap.map.joint'); ",
         "fit with hmm_map(method = 'joint') or hmm_map_joint().")
  if (!is.numeric(eta_range) || length(eta_range) != 2L ||
      eta_range[1] <= 0 || eta_range[1] >= eta_range[2])
    stop("`eta_range` must be a length-2 increasing positive interval.")

  ord  <- map$order
  fit  <- map$fit
  r    <- as.numeric(fit$r)
  eps  <- epsilon %||% fit$epsilon %||% 1e-3
  dams <- map$dams %||% names(map$phase_list)
  D    <- length(dams)
  T    <- length(ord)
  if (D < 2L) stop("Heterogeneity test needs at least 2 dams.")
  if (length(r) != T - 1L) stop("`map$fit$r` length must equal length(order) - 1.")

  # Haldane distances of the shared map, and the scaled-r map for a given eta
  rc <- pmin(pmax(r, 1e-9), 0.5 - 1e-9)
  m  <- -0.5 * log(1 - 2 * rc)
  r_scaled <- function(eta) pmin(pmax(0.5 * (1 - exp(-2 * eta * m)), 1e-9), 0.5 - 1e-9)

  # per-dam emission table (3 x T), aligned to the shared order
  emis_for_dam <- function(d) {
    em <- if (!is.null(fit$pi_emission_list)) fit$pi_emission_list[[d]]
          else if (!is.null(fit$pi_list))     fit$pi_list[[d]]
          else stop("Joint fit must contain `pi_list` or `pi_emission_list`.")
    em <- as.matrix(em)
    if (is.null(colnames(em))) colnames(em) <- ord else em <- em[, ord, drop = FALSE]
    em
  }

  # cache each dam's aligned inputs, then a likelihood evaluator at a given eta
  dam_data <- lapply(dams, function(d) {
    Gmat <- dat$G_list[[d]]; Mvec <- dat$M_list[[d]]
    if (is.null(Gmat)) stop("Dam '", d, "' not found in dat$G_list.")
    G <- matrix(NA_integer_, nrow(Gmat), T, dimnames = list(rownames(Gmat), ord))
    cg <- intersect(ord, colnames(Gmat)); if (length(cg)) G[, cg] <- Gmat[, cg, drop = FALSE]
    storage.mode(G) <- "integer"
    Mv <- rep(NA_integer_, T); names(Mv) <- ord
    cm <- intersect(ord, names(Mvec)); if (length(cm)) Mv[cm] <- as.integer(Mvec[cm])
    list(G = G, M = Mv, ph = as.integer(map$phase_list[[d]]), pe = emis_for_dam(d))
  })
  ll_d <- function(i, eta) {
    dd <- dam_data[[i]]
    loglik_hs_cpp(dd$G, dd$M, dd$ph, r_scaled(eta), dd$pe, eps)
  }

  # alternative: profile eta separately per dam
  alt <- lapply(seq_len(D), function(i)
    stats::optimize(function(e) ll_d(i, e), interval = eta_range, maximum = TRUE))
  eta_hat  <- vapply(alt, function(o) o$maximum,   numeric(1))
  ll_alt_d <- vapply(alt, function(o) o$objective, numeric(1))

  # null: a single common eta shared by all dams
  oc <- stats::optimize(function(e) sum(vapply(seq_len(D), function(i) ll_d(i, e), numeric(1))),
                        interval = eta_range, maximum = TRUE)
  eta_common <- oc$maximum
  ll_null_d  <- vapply(seq_len(D), function(i) ll_d(i, eta_common), numeric(1))

  per_dam <- data.frame(
    dam         = dams,
    n_offspring = vapply(dam_data, function(x) nrow(x$G), integer(1)),
    eta_hat     = eta_hat,
    ll_common   = ll_null_d,
    ll_alt      = ll_alt_d,
    lr          = 2 * (ll_alt_d - ll_null_d),
    stringsAsFactors = FALSE
  )
  LR <- sum(per_dam$lr)
  df <- D - 1L
  out <- list(
    eta         = stats::setNames(eta_hat, dams),
    eta_common  = eta_common,
    loglik_null = sum(ll_null_d),
    loglik_alt  = sum(ll_alt_d),
    per_dam     = per_dam,
    LR          = LR,
    df          = df,
    p_value     = stats::pchisq(LR, df = df, lower.tail = FALSE)
  )
  class(out) <- "HSMap.hetero"
  out
}

#' @export
print.HSMap.hetero <- function(x, ...) {
  cat("HSMap recombination-map heterogeneity test (per-dam Haldane scaling)\n")
  cat(sprintf("  %d dams | common scale eta = %.3f\n", nrow(x$per_dam), x$eta_common))
  cat(sprintf("  LR = %.2f on %d df,  p = %.3g\n", x$LR, x$df, x$p_value))
  pd <- x$per_dam
  pd$eta_hat <- round(pd$eta_hat, 3)
  pd$lr      <- round(pd$lr, 2)
  print(pd[, c("dam", "n_offspring", "eta_hat", "lr")], row.names = FALSE)
  invisible(x)
}
