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
#' @param tol Convergence tolerance for EM. Default \code{1e-4}.
#' @param pi_mode Emission parameterization for the paternal mixture when
#'   \code{paternal_mode != "two_locus"}. One of \code{c("per_marker","HWE")}.
#'   Default \code{"per_marker"}.
#' @param paternal_mode Paternal model for the HMM. One of
#'   \code{c("per_marker","HWE","two_locus")}. Default \code{"per_marker"}.
#' @param r_start Initial recombination fraction for all intervals.
#'   Default \code{0.05}.
#' @param lambda Dirichlet shrinkage strength for paternal parameters.
#'   Default \code{20}.
#' @param maxit Maximum EM iterations. Default \code{200}.
#' @param pi_prior_in Optional 3 x T prior for paternal genotype frequencies
#'   when \code{paternal_mode != "two_locus"}.
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
#' \code{dam}. If multiple dams are requested with \code{method = "joint"}
#' (default), a single shared map of class \code{c("HSMap.map.joint","HSMap.map")}
#' (see \code{\link{hmm_map_joint}}), directly usable by \code{\link{get_map}} and
#' \code{\link{plot_map_list}}. With \code{method = "consensus"}, a named list of
#' \code{"HSMap.map"} objects (class \code{"HSMap.map.multi"}); when
#' \code{return_consensus = TRUE} the list also includes \code{$consensus} with
#' fields \code{order}, \code{r}, and \code{weights}.
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
    tol = 1e-4,
    pi_mode = c("per_marker", "HWE"),
    paternal_mode = c("per_marker", "HWE", "two_locus"),
    r_start = 0.05,
    lambda = 20,
    maxit = 200,
    pi_prior_in = NULL,
    Pi_prior_in = NULL,
    method = c("joint", "consensus"),
    return_consensus = TRUE
) {
  `%||%` <- function(a, b) if (is.null(a)) b else a

  pi_mode <- match.arg(pi_mode)
  paternal_mode <- match.arg(paternal_mode)
  method <- match.arg(method)

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

    fit <- hmm_hs_cpp_parallel(
      G = G_sub,
      M = M_sub,
      phase_vec = as.integer(pv),
      r_start = r_start,
      pi_mode = pi_mode,
      pi_prior_in = if (is.null(pi_prior_in)) NULL else pi_prior_in,
      lambda = lambda,
      epsilon = epsilon,
      tol = tol,
      maxit = maxit,
      paternal_mode = paternal_mode,
      Pi_prior_in = if (is.null(Pi_prior_in)) NULL else Pi_prior_in
    )

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
      lambda = lambda, maxit = maxit
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
