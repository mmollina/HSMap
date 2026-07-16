#' Multipoint HMM map per half-sib family (maternal-only)
#'
#' @description
#' Fit a maternal-only multipoint HMM along a supplied marker order for one or
#' more **selected** dams. The model is always fit **separately for each dam**.
#' When several dams are requested, the function iterates over them, fits one
#' HMM per family, and returns the per-dam results together. If
#' `return_consensus = TRUE`, a consensus recombination fraction (`r`) per
#' interval is also computed across the selected dams (see *Consensus r*).
#'
#' @param x An \code{HSMap.data} object containing \code{G_list} and \code{M_list}.
#' @param phased An object (e.g., from \code{hsmap_phase_from_pairwise}) with
#'   fields \code{order} and \code{phase_vec}. It can represent one dam
#'   (\code{HSMap.phased}) or multiple dams (\code{HSMap.phased.multi}).
#' @param dam Which dam(s) to fit. One of:
#'   \itemize{
#'     \item a single integer index into \code{x$G_list},
#'     \item a single character name matching \code{names(x$G_list)},
#'     \item a vector of such indices and/or names, or
#'     \item \code{"all"} to fit every dam present in \code{x}.
#'   }
#'   Default is \code{1} (the most common use: one family).
#' @param threads Integer number of threads for \pkg{RcppParallel}. If
#'   \code{NULL}, the current setting is not changed.
#' @param epsilon Genotyping error rate used in emissions (in \code{[0,1)}).
#'   Default \code{0.01}.
#' @param tol Convergence tolerance for EM: relative change in the active objective
#'   (observed log-likelihood, or penalized objective when \code{lambda > 0}) and the
#'   maximum change in the identifiable parameters \code{r} and \code{q}. Default \code{1e-6}.
#' @param pi_mode Retained for backward compatibility only and **ignored** by the
#'   identifiable paternal model (the paternal contribution is a single per-marker
#'   allele frequency; see \code{paternal_mode}).
#' @param paternal_mode Paternal model for the HMM. One of
#'   \code{c("gametic","HWE","per_marker","two_locus")}, default \code{"gametic"}.
#'   \describe{
#'     \item{\code{"gametic"}}{The identifiable model. The paternal contribution
#'       is a single per-marker (and per-dam) paternal gametic frequency
#'       \eqn{q_k = \Pr(\text{paternal gamete transmits } A) = \pi_{AA}+\tfrac12\pi_{Aa}},
#'       the only paternal quantity the half-sib offspring likelihood identifies.
#'       Returned in \code{fit$q}.}
#'     \item{\code{"HWE"}}{Accepted compatibility alias for \code{"gametic"}:
#'       mathematically identical (under Hardy--Weinberg the allele frequency
#'       \eqn{p_k \equiv q_k}), so it yields byte-identical results. It is
#'       \emph{not} a separately identifiable biological model.}
#'     \item{\code{"per_marker"}}{Deprecated. Accepted with a warning and routed
#'       to \code{"gametic"}; any supplied \code{pi_prior_in} is collapsed to its
#'       induced \eqn{q}. The three free genotype frequencies are not identifiable
#'       (only \eqn{q_k} is). The legacy engine remains reachable via the internal
#'       \code{hmm_hs_cpp_parallel()} for historical reproduction.}
#'     \item{\code{"two_locus"}}{Disabled (informative error). The 10-class
#'       two-locus mixture does not identify paternal linkage disequilibrium.}
#'   }
#' @param q_prior_in Optional gametic prior on the paternal gametic frequency
#'   \eqn{q_k}, supplied as \strong{pseudocounts}. One of: \code{NULL} (default;
#'   uses the historical strength \code{lambda = 20}, i.e. \eqn{\alpha=\beta=10},
#'   pseudocount target \eqn{0.5}); a numeric scalar or length-\eqn{T} vector of
#'   pseudocount targets \eqn{q^{(0)}_k} (total pseudocount \code{lambda}); or
#'   \code{list(alpha=, beta=)} giving explicit non-negative pseudocounts (scalars
#'   or length-\eqn{T}). Per-marker targets are supported; the total pseudocount
#'   \eqn{\alpha+\beta} must be constant across markers. The paternal M-step is the
#'   penalized (MAP) update
#'   \deqn{\hat q_k = (N_{A,k}+\alpha_k)/(N_{A,k}+N_{a,k}+\alpha_k+\beta_k),}
#'   where \eqn{N_{A,k}, N_{a,k}} are the expected paternal A-/a-gamete counts;
#'   it maximizes \eqn{\log L + \sum_k[\alpha_k\log q_k + \beta_k\log(1-q_k)]}.
#'   \strong{\eqn{\alpha,\beta} are pseudocounts, not Beta shape parameters}: the
#'   equivalent probability prior is \eqn{\mathrm{Beta}(\alpha_k+1,\beta_k+1)}, so
#'   \eqn{\alpha=\beta=0} is the uniform \eqn{\mathrm{Beta}(1,1)} prior (MAP =
#'   unregularized fit). A numeric \eqn{q^{(0)}_k} is the \strong{pseudocount
#'   target (the prior mode} \eqn{q^{(0)}_k=\alpha_k/(\alpha_k+\beta_k)}, \strong{not
#'   the Beta-prior mean} \eqn{(\alpha_k+1)/(\alpha_k+\beta_k+2)}): it sets
#'   \eqn{\alpha_k=\lambda q^{(0)}_k} and \eqn{\beta_k=\lambda(1-q^{(0)}_k)}. Takes
#'   precedence over \code{pi_prior_in}.
#' @param r_start Initial recombination fraction for all intervals.
#'   Default \code{0.05}.
#' @param lambda Default total pseudocount \eqn{\alpha+\beta} for the gametic
#'   paternal prior when \code{q_prior_in} is \code{NULL} or a numeric target.
#'   Default \code{20} (\eqn{\alpha=\beta=10}, pseudocount target \eqn{0.5}). This
#'   value is the historical default and is kept for backward compatibility only;
#'   it is \strong{not} a statistically validated recommendation. It can over-shrink
#'   \eqn{q} at extreme paternal gametic frequencies (a separate simulation study on
#'   the corrected engine is needed to choose a default). See \code{q_prior_in};
#'   set \eqn{\alpha=\beta=0} for no regularization.
#' @param maxit Maximum EM iterations. Default \code{1000}. The R wrapper warns if
#'   the EM has not converged within \code{maxit}.
#' @param pi_prior_in Optional 3 x T paternal prior (legacy). Only its implied
#'   gametic allele frequency \eqn{q^{(0)}_k = \pi^{(0)}_{AA} +
#'   \tfrac12\pi^{(0)}_{Aa}} is used (the genotype split is not identifiable);
#'   prefer \code{q_prior_in}, which supersedes it.
#' @param Pi_prior_in Optional 10 x (T-1) prior for interval classes when
#'   \code{paternal_mode == "two_locus"}.
#' @param method Multi-dam estimator. \code{"joint"} (default) fits one shared
#'   recombination map by pooling the expected recombination counts across dams
#'   within a single EM. This is the likelihood-based estimator implied by the
#'   model (a proper multilocus likelihood, compatible with likelihood-based
#'   inference); see \code{\link{hmm_map_joint}}. \code{"consensus"} is the legacy
#'   per-dam fit followed by an offspring-weighted average of the per-dam
#'   estimates; it behaves like a shrinkage toward \code{r_start} and is retained
#'   as a diagnostic/comparison only. Ignored when a single dam is requested.
#' @param return_consensus Logical. Only used when \code{method = "consensus"}.
#'   If \code{TRUE} and multiple dams are requested, also compute a consensus
#'   \code{r} across the selected dams (see *Consensus r*). Default \code{TRUE}.
#'
#' @return
#' If a single dam is requested, an object of class \code{"HSMap.map"} with
#' \code{order}, \code{phase_vec}, \code{fit} (the C++ HMM result) and
#' \code{dam}. The canonical paternal output is \code{fit$q}, the per-marker sire
#' gametic allele frequency \eqn{q_k}. \strong{\code{fit$pi} is deprecated for
#' direct interpretation}: it is retained only as the derived HWE-form emission
#' table (\eqn{[q^2, 2q(1-q), (1-q)^2]}) consumed by downstream decoders, and is
#' \emph{not} an estimate of paternal genotype frequencies. If multiple dams are requested
#' with \code{method = "joint"}
#' (default), a single shared map of class \code{c("HSMap.map.joint","HSMap.map")}
#' (see \code{\link{hmm_map_joint}}), directly usable by \code{\link{get_map}} and
#' \code{\link{plot_map_list}}. With \code{method = "consensus"}, a named list of
#' \code{"HSMap.map"} objects (class \code{"HSMap.map.multi"}); when
#' \code{return_consensus = TRUE} the list also includes \code{$consensus} with
#' fields \code{order}, \code{r}, and \code{weights}.
#'
#' @section Paternal identifiability:
#' The offspring likelihood depends on the paternal contribution \emph{only}
#' through the paternal gametic frequency \eqn{q_k = \pi_{AA} + \tfrac12
#' \pi_{Aa}}: any two genotype-frequency vectors with the same \eqn{q} give the
#' same likelihood and the same fitted \eqn{\mathbf r}. The \emph{regularized}
#' estimate does depend on the prior. The \code{"gametic"} model applies a
#' penalty \eqn{\alpha\log q + \beta\log(1-q)} (a \eqn{\mathrm{Beta}(\alpha+1,
#' \beta+1)} prior; see \code{q_prior_in}), which \emph{need not} coincide with
#' the legacy \code{"per_marker"} Dirichlet-on-genotypes regularizer even though
#' both share the same unpenalized objective. Setting \eqn{\alpha=\beta=0} removes
#' the penalty (uniform \eqn{\mathrm{Beta}(1,1)} prior), giving the unregularized
#' fit; because the objective may have multiple optima this is not asserted to be
#' the global MLE.
#'
#' @section Consensus r:
#' For several dams fitted on the same marker order, we form a consensus
#' recombination fraction per interval by a simple weighted average:
#'
#' \preformatted{r_hat_t = sum_d w_t(d) * r_hat_t(d) / sum_d w_t(d)}
#'
#' where \code{r_hat_t(d)} is the interval estimate for dam \code{d} and
#' \code{w_t(d)} is an information weight for that dam at that interval.
#' A practical proxy right now: weight by the dam’s LOD support per interval
#' or by the number of offspring with observed genotypes at the two loci; both
#' correlate well with information. The implementation here uses the latter
#' (number of offspring with data at both markers) as \code{w_t(d)}. This is a
#' standard composite-likelihood-style estimator and behaves well in practice.
#'
#' @examples
#' \dontrun{
#' dat <- read_HSMap_data("ped.csv", "geno.csv")
#' # ph can be a single-dam HSMap.phased or a multi-dam HSMap.phased.multi
#'
#' # Single family (default dam = 1)
#' map1 <- hmm_map(dat, phased = ph_single, dam = 1)
#'
#' # Several families by name
#' maps <- hmm_map(dat, phased = ph_multi, dam = c("MOM1","MOM3"))
#'
#' # All dams with consensus r
#' maps_all <- hmm_map(dat, phased = ph_multi, dam = "all", return_consensus = TRUE)
#' }
#'
#' @export
hmm_map <- function(
    x,
    phased,
    dam = 1,
    threads = NULL,
    epsilon = 0.01,
    tol = 1e-6,
    pi_mode = c("per_marker", "HWE"),
    paternal_mode = c("gametic", "HWE", "per_marker", "two_locus"),
    r_start = 0.05,
    lambda = 20,
    maxit = 1000,
    q_prior_in = NULL,
    pi_prior_in = NULL,
    Pi_prior_in = NULL,
    method = c("joint", "consensus"),
    return_consensus = TRUE
) {
  `%||%` <- function(a, b) if (is.null(a)) b else a

  pi_mode <- match.arg(pi_mode)
  paternal_mode <- match.arg(paternal_mode)
  method <- match.arg(method)

  # M1: resolve the public paternal model to the internal engine. `gametic`
  # (default) and `HWE` are the same identifiable estimator (parameterized by
  # q_k); `per_marker` warns and routes to `gametic`, while `two_locus` is disabled.
  eff_paternal <- .hsmap_paternal_engine(paternal_mode)

  if (!inherits(x, "HSMap.data"))
    stop("`x` must be HSMap.data (see read_HSMap_data).")

  if (!is.null(threads)) {
    if (!is.numeric(threads) || length(threads) != 1L || threads < 1)
      stop("`threads` must be a positive integer.")
    RcppParallel::setThreadOptions(numThreads = as.integer(threads))
  }

  # Normalize phased input to a named list of per-dam objects
  phased_list <- if (inherits(phased, "HSMap.phased.multi")) {
    phased
  } else if (inherits(phased, "HSMap.phased")) {
    nm <- phased$dam %||% "Dam1"
    stats::setNames(list(phased), nm)
  } else {
    stop("`phased` must be HSMap.phased (single) or HSMap.phased.multi (many).")
  }

  dam_names <- names(x$G_list)
  if (is.null(dam_names)) dam_names <- as.character(seq_along(x$G_list))

  # Resolve requested dams -> names
  if (identical(dam, "all")) {
    dams_req <- dam_names
  } else {
    dam_idx <- dam
    # allow mixture of indices and names
    if (is.numeric(dam_idx)) {
      dam_idx <- as.integer(dam_idx)
      if (any(dam_idx < 1L | dam_idx > length(dam_names)))
        stop("`dam` indices out of range.")
      dams_req <- dam_names[dam_idx]
    } else if (is.character(dam)) {
      if (!all(dam %in% dam_names))
        stop("Unknown dam name(s): ", paste(setdiff(dam, dam_names), collapse = ", "))
      dams_req <- dam
    } else {
      stop("`dam` must be an index (or indices), name(s), or 'all'.")
    }
  }

  # Keep only requested dams from phased_list
  if (!all(dams_req %in% names(phased_list))) {
    missing_ph <- setdiff(dams_req, names(phased_list))
    stop("`phased` does not contain phase/order for dam(s): ",
         paste(missing_ph, collapse = ", "))
  }
  phased_list <- phased_list[dams_req]

  run_one <- function(ph) {
    dam_id <- ph$dam %||% names(phased_list)[[1]]
    if (!(dam_id %in% dam_names))
      stop("Dam '", dam_id, "' not found in HSMap.data.")

    Gmat <- x$G_list[[dam_id]]
    Mvec <- x$M_list[[dam_id]]
    if (is.null(Gmat) || !is.matrix(Gmat))
      stop("x$G_list[['", dam_id, "']] must be a matrix: offspring x markers.")
    if (is.null(Mvec) || is.null(names(Mvec)))
      stop("x$M_list[['", dam_id, "']] must be a named vector.")

    ord <- ph$order
    pv  <- ph$phase_vec
    if (length(pv) != length(ord) - 1L)
      stop("phase_vec length must be length(order)-1 for dam '", dam_id, "'.")
    # Never run the HMM through unresolved (NA) phase intervals: the model requires a
    # known relative phase for every fitted interval. Split the map first.
    if (anyNA(pv))
      stop("phase_vec for dam '", dam_id, "' has ", sum(is.na(pv)),
           " unresolved interval(s) (NA phase). hmm_map() requires a fully resolved ",
           "phase; use hmm_map_blocks() to fit resolved phase blocks separately, or ",
           "restrict `order` to a single resolved block.", call. = FALSE)

    missG <- setdiff(ord, colnames(Gmat))
    missM <- setdiff(ord, names(Mvec))
    if (length(missG))
      stop("Markers not in x$G_list[['", dam_id, "']]: ",
           paste(utils::head(missG, 8), collapse = ", "),
           if (length(missG) > 8) " ...")
    if (length(missM))
      stop("Markers not in x$M_list[['", dam_id, "']]: ",
           paste(utils::head(missM, 8), collapse = ", "),
           if (length(missM) > 8) " ...")

    G_sub <- Gmat[, ord, drop = FALSE]  # offspring x markers
    storage.mode(G_sub) <- "integer"
    M_sub <- as.integer(Mvec[ord])

    # Resolve the paternal prior. A gametic pseudocount prior (`q_prior_in`) takes
    # precedence over a legacy 3xT `pi_prior_in` (which the HWE engine collapses
    # to its induced q anyway).
    if (!is.null(q_prior_in)) {
      eng <- .hsmap_q_prior_to_engine(q_prior_in, lambda, length(ord))
      pi_prior_eff <- eng$pi_prior; lambda_eff <- eng$lambda
      if (!is.null(pi_prior_in))
        warning("Both `q_prior_in` and `pi_prior_in` supplied; using `q_prior_in`.", call. = FALSE)
    } else {
      pi_prior_eff <- pi_prior_in; lambda_eff <- lambda
    }

    fit <- hmm_hs_cpp_parallel(
      G = G_sub,
      M = M_sub,
      phase_vec = as.integer(pv),
      r_start = r_start,
      pi_mode = "HWE",
      pi_prior_in = pi_prior_eff,
      lambda = lambda_eff,
      epsilon = epsilon,
      tol = tol,
      maxit = maxit,
      paternal_mode = eff_paternal,
      Pi_prior_in = if (is.null(Pi_prior_in)) NULL else Pi_prior_in
    )

    if (isFALSE(fit$converged))
      warning("hmm_map(): EM did not converge for dam '", dam_id, "' in ", fit$iters,
              " iterations (reason: ", fit$conv_reason, "). Increase `maxit` or relax ",
              "`tol`; the returned estimates are the last iterate.", call. = FALSE)
    if (isTRUE(fit$objective_decreased))
      warning("hmm_map(): the active EM objective decreased materially for dam '",
              dam_id, "'; results may be unreliable.", call. = FALSE)

    # Canonical identifiable paternal output: q_k = P(paternal gamete transmits A)
    # = pi_AA + 0.5*pi_Aa. `fit$pi` is retained only as the derived HWE-form
    # emission table used by downstream decoders; it is NOT an estimate of
    # paternal genotype frequencies.
    fit$q <- stats::setNames(as.numeric(fit$pi["AA", ] + 0.5 * fit$pi["Aa", ]), ord)
    fit$paternal_mode <- paternal_mode

    out <- list(
      dam       = dam_id,
      order     = ord,
      phase_vec = as.integer(pv),
      fit       = fit
    )
    class(out) <- "HSMap.map"
    out
  }

  if (length(phased_list) == 1L) {
    return(run_one(phased_list[[1L]]))
  }

  # Default multi-dam behavior: the joint, likelihood-based shared-map estimator.
  if (identical(method, "joint")) {
    ph_multi <- phased_list
    class(ph_multi) <- "HSMap.phased.multi"
    return(hmm_map_joint(
      x, phased = ph_multi, dam = "all", threads = NULL,
      epsilon = epsilon, tol = tol, pi_mode = pi_mode,
      paternal_mode = paternal_mode, r_start = r_start,
      lambda = lambda, maxit = maxit, q_prior_list = q_prior_in
    ))
  }

  # method == "consensus": legacy per-dam fits + offspring-weighted average (diagnostic)
  per_dam <- lapply(phased_list, run_one)

  # Optional consensus r using weights = #offspring with observed genotypes at both loci
  out <- list(per_dam = per_dam)
  if (isTRUE(return_consensus)) {
    ord0 <- per_dam[[1]]$order
    K <- length(ord0) - 1L
    r_sum   <- numeric(K)
    weights <- numeric(K)

    for (nm in names(per_dam)) {
      Gd <- x$G_list[[nm]][, ord0, drop = FALSE]
      w_d <- vapply(seq_len(ncol(Gd) - 1L), function(t)
        sum(!is.na(Gd[, t]) & !is.na(Gd[, t + 1L])), numeric(1))
      rd <- as.numeric(per_dam[[nm]]$fit$r)
      r_sum   <- r_sum   + w_d * rd
      weights <- weights + w_d
    }
    r_cons <- ifelse(weights > 0, r_sum / weights, NA_real_)
    out$consensus <- list(order = ord0, r = r_cons, weights = weights)
  }

  class(out) <- "HSMap.map.multi"
  out
}

# Resolve the public `paternal_mode` to the internal C++ engine mode (M1).
#
# The paternal contribution enters the likelihood only through the sire gametic
# allele frequency q_k = P(paternal gamete transmits A) = pi_AA + 0.5*pi_Aa.
# `gametic` (default) and `HWE` are the same identifiable estimator and both use
# the HWE engine, which is parameterized by a single per-marker allele frequency
# p_k identical to q_k. `per_marker` (free 3 genotype frequencies) is deprecated: it
# warns and routes to `gametic`. `two_locus` (10-class interval mixture) is disabled.
# Both are non-identifiable; the legacy engines remain reachable via the internal
# `hmm_hs_cpp_parallel()` for reproducing historical results.
.hsmap_paternal_engine <- function(paternal_mode) {
  if (identical(paternal_mode, "two_locus"))
    stop("paternal_mode = 'two_locus' is disabled: the two-locus paternal mixture does ",
         "not identify paternal linkage disequilibrium (coupling and repulsion ",
         "double-heterozygous paternal genotypes are likelihood-indistinguishable). ",
         "Use paternal_mode = 'gametic'.", call. = FALSE)
  if (identical(paternal_mode, "per_marker"))
    warning("paternal_mode = 'per_marker' is deprecated: the three paternal genotype ",
            "frequencies are not identifiable from half-sib offspring (only ",
            "q_k = P(paternal gamete transmits A) = pi_AA + 0.5*pi_Aa is identified). ",
            "Routing to the identifiable 'gametic' model; any supplied pi_prior is ",
            "collapsed to its induced q. The legacy engine remains available via ",
            "HSMap:::hmm_hs_cpp_parallel(paternal_mode = 'per_marker') for historical ",
            "reproduction.", call. = FALSE)
  # gametic, HWE, and (deprecated) per_marker all route to the HWE engine (p_k == q_k).
  "HWE"
}

# Convert a gametic pseudocount prior on q_k (alpha, beta) into the internal HWE
# engine's (pi_prior 3xT, lambda) representation. The engine's per-marker paternal
# M-step
#
#   q_k^{new} = (N_A,k + alpha_k) / (N_A,k + N_a,k + alpha_k + beta_k)
#
# maximizes the penalized objective logLik + sum_k[alpha_k log q_k +
# beta_k log(1 - q_k)] over the expected paternal A-/a-gamete counts (N_A, N_a);
# equivalently it is the posterior MODE under a Beta(alpha_k + 1, beta_k + 1)
# prior (alpha, beta are pseudocounts, not Beta shape parameters). Here
# alpha_k = lambda * q0_k, beta_k = lambda * (1 - q0_k), and the total pseudocount
# lambda = alpha_k + beta_k. Here q0_k = alpha_k/(alpha_k + beta_k) is the pseudocount
# target (the prior MODE), not the Beta-prior mean (alpha_k+1)/(alpha_k+beta_k+2).
# Per-marker targets q0_k are supported; the total pseudocount (alpha + beta) must
# be common across markers (the engine takes a scalar lambda).
#
# Accepts: NULL (use lambda_default, engine's uniform prior); a numeric scalar or
# length-T vector of pseudocount targets q0 (total pseudocount = lambda_default); or
# list(alpha=, beta=) with scalar or length-T entries (explicit pseudocounts).
.hsmap_q_prior_to_engine <- function(q_prior, lambda_default, T) {
  if (is.null(q_prior)) return(list(pi_prior = NULL, lambda = lambda_default))

  hwe_cols <- function(q0) {
    q0 <- pmin(pmax(as.numeric(q0), 1e-9), 1 - 1e-9)
    matrix(rbind(q0^2, 2 * q0 * (1 - q0), (1 - q0)^2), nrow = 3,
           dimnames = list(c("AA", "Aa", "aa"), NULL))
  }

  if (is.numeric(q_prior)) {
    q0 <- if (length(q_prior) == 1L) rep(q_prior, T) else q_prior
    if (length(q0) != T) stop("numeric `q_prior` must have length 1 or T = ", T, ".", call. = FALSE)
    if (any(!is.finite(q0) | q0 <= 0 | q0 >= 1))
      stop("`q_prior` pseudocount targets must be in (0, 1).", call. = FALSE)
    return(list(pi_prior = hwe_cols(q0), lambda = lambda_default))
  }

  if (is.list(q_prior) && all(c("alpha", "beta") %in% names(q_prior))) {
    a <- if (length(q_prior$alpha) == 1L) rep(as.numeric(q_prior$alpha), T) else as.numeric(q_prior$alpha)
    b <- if (length(q_prior$beta)  == 1L) rep(as.numeric(q_prior$beta),  T) else as.numeric(q_prior$beta)
    if (length(a) != T || length(b) != T)
      stop("`alpha`/`beta` must have length 1 or T = ", T, ".", call. = FALSE)
    if (any(!is.finite(a) | !is.finite(b) | a < 0 | b < 0))
      stop("`alpha`/`beta` must be finite and non-negative.", call. = FALSE)
    conc <- a + b
    if (all(conc == 0)) return(list(pi_prior = hwe_cols(rep(0.5, T)), lambda = 0))  # no shrinkage
    if (any(conc == 0))
      stop("mixed zero / non-zero total pseudocount (alpha + beta) across markers is not supported.", call. = FALSE)
    if (diff(range(conc)) > 1e-8 * max(conc))
      stop("per-marker total pseudocount (alpha + beta) must be constant across markers ",
           "in this version; vary the target alpha/(alpha+beta) but keep alpha+beta fixed.",
           call. = FALSE)
    return(list(pi_prior = hwe_cols(a / conc), lambda = conc[1]))
  }

  stop("`q_prior` must be NULL, a numeric vector of pseudocount targets, or list(alpha=, beta=).",
       call. = FALSE)
}
