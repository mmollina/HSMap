#' Conditional global-scale test for map heterogeneity across dams
#'
#' @description
#' A \strong{conditional global-scale likelihood-ratio test} (LRT). After fitting a
#' joint shared map, it asks whether dams differ by a single per-dam multiplier of
#' the \emph{whole} map, \strong{not} whether individual intervals are heterogeneous.
#' For dam \eqn{d} the alternative scales every Haldane distance by one factor,
#' \eqn{m_k^{(d)} = \eta^{(d)} m_k}, so
#' \deqn{r_k^{(d)} = \tfrac12\bigl\{1 - \exp(-2\,\eta^{(d)} m_k)\bigr\}.}
#'
#' \strong{Hypotheses.}
#' \itemize{
#'   \item \strong{Null:} one common global scale \eqn{\eta} shared by all dams
#'     (a single free parameter).
#'   \item \strong{Alternative:} a separate global scale \eqn{\eta^{(d)}} per dam
#'     (\eqn{D} free parameters).
#' }
#' This is a \strong{global}-scale contrast: it does \strong{not} test
#' interval-specific heterogeneity (the interval \emph{shape} of the shared map is
#' held fixed; only its overall length is rescaled per dam).
#'
#' @details
#' \strong{Conditional statement.} Each \eqn{\eta^{(d)}} is profiled by a
#' one-dimensional search with the shared interval shape \eqn{m_k}, the per-dam
#' linkage phase, and the per-dam paternal nuisance parameters (\eqn{q_k^{(d)}}, via
#' the fitted emission tables) \strong{held fixed}. The likelihood-ratio statistic is
#' therefore \emph{conditional} on the fitted phase and paternal parameters; they are
#' not re-estimated under the scaled alternative. The statistic is
#' \deqn{\mathrm{LR} = 2\Bigl[\textstyle\sum_d \ell_d(\hat\eta^{(d)}) -
#'       \textstyle\sum_d \ell_d(\hat\eta_{\mathrm{common}})\Bigr]}
#' on \eqn{D-1} degrees of freedom (\eqn{D} scales under the alternative minus the
#' one common scale under the null).
#'
#' \strong{Calibration.} The chi-square p-value is \emph{asymptotic}. For small
#' \eqn{D}, small families, or \eqn{\hat\eta} near the search boundary, the
#' asymptotic calibration can be inaccurate and a \strong{parametric bootstrap}
#' should be used for a definitive p-value; a bootstrap is not performed here. The
#' per-dam \eqn{\hat\eta^{(d)}} are the primary diagnostic (\eqn{>1} = longer map than
#' shared).
#'
#' @param dat An \code{HSMap.data} object.
#' @param map A joint shared map of class \code{"HSMap.map.joint"}.
#' @param epsilon Optional override for the emission error rate; defaults to the
#'   value stored in the fit.
#' @param eta_range Length-2 positive, increasing search interval for \eqn{\eta}.
#'   Default \code{c(0.1, 10)}.
#' @param gap_r No-linkage threshold (default \code{0.499}). The test \strong{rejects}
#'   a map with any \code{r >= gap_r} (or non-finite \code{r}, unresolved phase, or a
#'   non-converged fit): the global Haldane scaling is only defined for a fully
#'   linked, resolved, converged map/block. \code{r = 0.5} is not silently replaced.
#'
#' @return An object of class \code{"HSMap.hetero"}: a list with \code{test_type}
#'   (\code{"conditional global-scale LRT"}), \code{eta} (named \eqn{\hat\eta^{(d)}}),
#'   \code{eta_common}, \code{loglik_null}, \code{loglik_alt}, \code{per_dam} (data
#'   frame incl. an \code{at_boundary} flag), \code{LR}, \code{df}, \code{p_value},
#'   \code{n_params_null} (1), \code{n_params_alt} (\eqn{D}), \code{calibration}
#'   (\code{"asymptotic"}), \code{conditional_on}, \code{interval_specific}
#'   (\code{FALSE}), and \code{any_boundary}.
#'
#' @seealso \code{\link{hmm_map_joint}}
#' @export
test_map_heterogeneity <- function(dat, map, epsilon = NULL, eta_range = c(0.1, 10),
                                   gap_r = 0.499) {
  `%||%` <- function(a, b) if (is.null(a)) b else a
  if (!inherits(dat, "HSMap.data")) stop("`dat` must be an HSMap.data object.")
  if (!inherits(map, "HSMap.map.joint"))
    stop("`map` must be a joint shared map (class 'HSMap.map.joint'); ",
         "fit with hmm_map(method = 'joint') or hmm_map_joint().")
  if (!is.numeric(eta_range) || length(eta_range) != 2L ||
      eta_range[1] <= 0 || eta_range[1] >= eta_range[2])
    stop("`eta_range` must be a length-2 increasing positive interval.")
  if (!is.numeric(gap_r) || length(gap_r) != 1L || !is.finite(gap_r) || gap_r <= 0 || gap_r > 0.5)
    stop("`gap_r` must be a single number in (0, 0.5].")

  ord  <- map$order
  fit  <- map$fit
  r    <- as.numeric(fit$r)
  eps  <- epsilon %||% fit$epsilon %||% 1e-3
  dams <- map$dams %||% names(map$phase_list)
  D    <- length(dams)
  T    <- length(ord)
  if (D < 2L) stop("Heterogeneity test needs at least 2 dams.")
  if (length(r) != T - 1L) stop("`map$fit$r` length must equal length(order) - 1.")

  # The conditional global-scale test scales Haldane distances, which is only defined
  # for a fully LINKED, RESOLVED, CONVERGED map. Reject anything else rather than
  # silently transforming a no-linkage interval (r = 0.5 is NOT replaced by 0.5-eps).
  msg_tail <- paste0(" Apply the heterogeneity test only to valid resolved, linked, ",
                     "converged joint blocks (e.g. a single hmm_map_blocks() block).")
  if (any(!is.finite(r)))
    stop("`map` has non-finite recombination fractions (r).", msg_tail, call. = FALSE)
  if (any(r >= gap_r))
    stop("`map` contains at least one no-linkage interval (r >= gap_r = ", gap_r,
         "); the global Haldane scaling is undefined there.", msg_tail, call. = FALSE)
  if (!is.null(map$resolved_interval) && any(!as.logical(map$resolved_interval)))
    stop("`map` contains unresolved-phase intervals.", msg_tail, call. = FALSE)
  if (isFALSE(fit$converged))
    stop("`map` was fit by an EM that did not converge (converged = FALSE).", msg_tail,
         call. = FALSE)
  if (isTRUE(fit$objective_decreased))
    stop("`map` was fit by an EM whose active objective decreased materially ",
         "(objective_decreased = TRUE); the fit is unreliable.", msg_tail, call. = FALSE)

  # Haldane distances of the shared map (finite, since the input r < gap_r < 0.5).
  # The scaled map uses the biological formula r_scaled(eta) = 0.5*(1 - exp(-2*eta*m));
  # it is NOT capped at gap_r (gap_r is only an input-map validity/reporting threshold).
  # A scaled r may exceed gap_r when eta > 1 --- that is the whole point of the test.
  # Only a numerical guard at 0.5 - 1e-12 (independent of gap_r) is applied.
  m  <- -0.5 * log(1 - 2 * r)
  r_scaled <- function(eta) pmin(0.5 * (1 - exp(-2 * eta * m)), 0.5 - 1e-12)

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

  # Boundary / convergence diagnostics: an eta at the search boundary is unreliable
  # (tolerance generous enough for optimize()'s own precision).
  btol <- 0.01 * (eta_range[2] - eta_range[1])
  at_boundary <- (eta_hat <= eta_range[1] + btol) | (eta_hat >= eta_range[2] - btol)
  common_boundary <- (eta_common <= eta_range[1] + btol) | (eta_common >= eta_range[2] - btol)
  per_dam$at_boundary <- at_boundary

  out <- list(
    test_type      = "conditional global-scale LRT",
    eta            = stats::setNames(eta_hat, dams),
    eta_common     = eta_common,
    loglik_null    = sum(ll_null_d),
    loglik_alt     = sum(ll_alt_d),
    per_dam        = per_dam,
    LR             = LR,
    df             = df,
    p_value        = stats::pchisq(LR, df = df, lower.tail = FALSE),
    n_params_null  = 1L,          # one common global scale
    n_params_alt   = D,           # one global scale per dam
    calibration    = "asymptotic",
    conditional_on = c("fitted linkage phase", "fitted paternal q (emission tables)",
                       "shared interval shape (m_k)"),
    interval_specific = FALSE,    # this is a GLOBAL-scale test, not interval-specific
    eta_range      = eta_range,
    any_boundary   = any(at_boundary) || common_boundary,
    note           = paste0("Asymptotic chi-square p-value, conditional on the fitted ",
                            "phase and paternal parameters; use a parametric bootstrap ",
                            "for small D or boundary eta.")
  )
  class(out) <- "HSMap.hetero"
  out
}

#' @export
print.HSMap.hetero <- function(x, ...) {
  cat("HSMap conditional global-scale heterogeneity test (per-dam Haldane scaling)\n")
  cat("  Null: one common global map scale;  Alt: a per-dam global scale.\n")
  cat("  Global-scale test (NOT interval-specific); conditional on fitted phase and paternal q.\n")
  cat(sprintf("  %d dams | common scale eta = %.3f\n", nrow(x$per_dam), x$eta_common))
  cat(sprintf("  LR = %.2f on %d df,  p = %.3g  (%s calibration)\n",
              x$LR, x$df, x$p_value, x$calibration))
  if (isTRUE(x$any_boundary))
    cat("  NOTE: an eta estimate is at the search boundary; use a bootstrap.\n")
  pd <- x$per_dam
  pd$eta_hat <- round(pd$eta_hat, 3)
  pd$lr      <- round(pd$lr, 2)
  print(pd[, c("dam", "n_offspring", "eta_hat", "lr")], row.names = FALSE)
  invisible(x)
}
