#' Joint multi-dam HMM map (one shared recombination map across dams)
#'
#' @description
#' Estimate a **single recombination map** shared across several dams by maximum
#' likelihood, pooling information over all offspring of all dams while keeping
#' **phase and paternal parameters dam-specific**. This is the joint
#' Expectation-Maximization described in the manuscript (one shared \eqn{\mathbf r};
#' expected recombination counts summed over all offspring across dams), and it
#' replaces the per-dam-then-average ("consensus") behavior of [hmm_map()] for the
#' multi-dam case.
#'
#' All selected dams must share the **same marker order** (same markers, same
#' sequence) — typically the order produced by [phase_from_pairwise()] when called
#' with a common `order`. Markers absent in a given dam are allowed (filled with
#' `NA` and treated as uninformative emissions).
#'
#' @param x An \code{HSMap.data} object with \code{G_list} and \code{M_list}.
#' @param phased An \code{HSMap.phased.multi} (one entry per dam) or a single
#'   \code{HSMap.phased}; each must carry \code{order} and \code{phase_vec}.
#' @param dam Which dam(s) to include: \code{"all"} (default), a vector of names
#'   or indices into \code{phased}.
#' @param threads Integer threads for \pkg{RcppParallel}; \code{NULL} leaves the
#'   current setting unchanged.
#' @param epsilon Genotyping error rate in emissions. Default \code{0.01}.
#' @param tol EM convergence tolerance (on log-likelihood and \code{r}). Default \code{1e-4}.
#' @param paternal_mode Paternal model, one of
#'   \code{c("gametic","HWE","per_marker","two_locus")}, default \code{"gametic"}.
#'   \code{"gametic"} (and its identical alias \code{"HWE"}) parameterizes the
#'   paternal contribution by the single identifiable per-marker, per-dam sire
#'   gametic allele frequency \eqn{q_k = \pi_{AA} + \tfrac12\pi_{Aa}}, returned in
#'   \code{fit$q_list}. \code{"per_marker"} is deprecated (accepted with a warning
#'   and routed to \code{"gametic"}); \code{"two_locus"} is disabled with an
#'   informative error. See \code{\link{hmm_map}}.
#' @param pi_mode Retained for backward compatibility only and ignored by the
#'   identifiable paternal model.
#' @param q_prior_list Optional gametic pseudocount prior on \eqn{q_k}, shared
#'   across dams or per dam. Either a single spec applied to all dams (a numeric
#'   mean vector, or \code{list(alpha=, beta=)} pseudocounts; see
#'   \code{\link{hmm_map}}'s \code{q_prior_in}) or a per-dam list of such specs.
#'   A per-dam list must supply \strong{exactly one spec per requested dam}: a
#'   named list must name every requested dam and no others (missing or unknown
#'   names are an error); an unnamed list must have one entry per requested dam,
#'   in the requested-dam order. Per-dam prior means are supported; the total
#'   pseudocount \eqn{\alpha+\beta} must be common across dams. Takes precedence
#'   over \code{pi_prior_list}.
#' @param r_start Initial recombination fraction for all intervals. Default \code{0.05}.
#' @param lambda Default total pseudocount \eqn{\alpha+\beta} for the gametic
#'   paternal prior when \code{q_prior_list} is \code{NULL} or a numeric mean.
#'   Default \code{20} (\eqn{\alpha=\beta=10} at prior mean \eqn{0.5}); the
#'   historical default, kept for backward compatibility only and \strong{not} a
#'   statistically validated recommendation. See \code{\link{hmm_map}}'s
#'   \code{q_prior_in}.
#' @param maxit Maximum EM iterations. Default \code{200}.
#' @param pi_prior_list,Pi_prior_list Optional named lists of per-dam priors
#'   (\code{3 x T} for Model A, \code{10 x (T-1)} for \code{two_locus}).
#'
#' @return An object of class \code{c("HSMap.map.joint","HSMap.map")} with
#'   \code{order} (shared), \code{dams}, \code{phase_list} (per dam),
#'   \code{fit} (the C++ result: shared \code{r}, per-dam \code{pi_list}, etc.) and
#'   a top-level \code{r}. The canonical per-dam paternal output is
#'   \code{fit$q_list} (sire gametic allele frequencies \eqn{q_k}).
#'   \strong{\code{fit$pi_list} is deprecated for direct interpretation}: it is
#'   retained only as the derived HWE-form emission table and is \emph{not} an
#'   estimate of paternal genotype frequencies. Because it carries \code{fit$r}
#'   and \code{order}, it works directly with [get_map()] and [plot_map_list()].
#'
#' @seealso [hmm_map()] for single-dam fitting and the legacy consensus path.
#' @export
hmm_map_joint <- function(
    x,
    phased,
    dam = "all",
    threads = NULL,
    epsilon = 0.01,
    tol = 1e-4,
    pi_mode = c("per_marker", "HWE"),
    paternal_mode = c("gametic", "HWE", "per_marker", "two_locus"),
    r_start = 0.05,
    lambda = 20,
    maxit = 200,
    q_prior_list = NULL,
    pi_prior_list = NULL,
    Pi_prior_list = NULL
) {
  `%||%` <- function(a, b) if (is.null(a)) b else a
  pi_mode <- match.arg(pi_mode)
  paternal_mode <- match.arg(paternal_mode)

  # M1: resolve the public paternal model to the internal engine (see hmm_map).
  eff_paternal <- .hsmap_paternal_engine(paternal_mode)

  if (!inherits(x, "HSMap.data"))
    stop("`x` must be HSMap.data (see read_HSMap_data).")
  if (!is.null(threads)) {
    if (!is.numeric(threads) || length(threads) != 1L || threads < 1)
      stop("`threads` must be a positive integer.")
    RcppParallel::setThreadOptions(numThreads = as.integer(threads))
  }

  phased_list <- if (inherits(phased, "HSMap.phased.multi")) {
    phased
  } else if (inherits(phased, "HSMap.phased")) {
    stats::setNames(list(phased), phased$dam %||% "Dam1")
  } else {
    stop("`phased` must be HSMap.phased (single) or HSMap.phased.multi (many).")
  }

  dam_names <- names(x$G_list) %||% as.character(seq_along(x$G_list))

  # resolve requested dams -> names present in phased_list
  ph_names <- names(phased_list)
  if (identical(dam, "all")) {
    dams_req <- ph_names
  } else if (is.numeric(dam)) {
    idx <- as.integer(dam)
    if (any(idx < 1L | idx > length(ph_names))) stop("`dam` index out of range.")
    dams_req <- ph_names[idx]
  } else if (is.character(dam)) {
    if (!all(dam %in% ph_names))
      stop("Unknown dam(s) in `phased`: ", paste(setdiff(dam, ph_names), collapse = ", "))
    dams_req <- dam
  } else stop("`dam` must be 'all', name(s), or index/indices.")
  phased_list <- phased_list[dams_req]

  # common order across dams (required)
  ord <- phased_list[[1]]$order
  if (is.null(ord)) stop("phased[[1]]$order is missing.")
  for (nm in dams_req) {
    if (!identical(phased_list[[nm]]$order, ord))
      stop("All dams must share the same marker order for joint estimation. ",
           "Phase with a common `order` (e.g. phase_from_pairwise(..., order=common)) first.")
  }
  K <- length(ord) - 1L

  # build per-dam aligned G, M, phase
  G_list <- vector("list", length(dams_req))
  M_list <- vector("list", length(dams_req))
  phase_list <- vector("list", length(dams_req))
  names(G_list) <- names(M_list) <- names(phase_list) <- dams_req

  for (nm in dams_req) {
    if (!(nm %in% dam_names))
      stop("Dam '", nm, "' not found in HSMap.data (G_list).")
    Gmat <- x$G_list[[nm]]
    Mvec <- x$M_list[[nm]]
    if (is.null(Gmat) || !is.matrix(Gmat)) stop("x$G_list[['", nm, "']] must be a matrix.")
    if (is.null(names(Mvec))) stop("x$M_list[['", nm, "']] must be a named vector.")

    Gsub <- matrix(NA_integer_, nrow(Gmat), length(ord),
                   dimnames = list(rownames(Gmat), ord))
    commonG <- intersect(ord, colnames(Gmat))
    if (length(commonG)) Gsub[, commonG] <- Gmat[, commonG, drop = FALSE]
    storage.mode(Gsub) <- "integer"

    Msub <- rep(NA_integer_, length(ord)); names(Msub) <- ord
    commonM <- intersect(ord, names(Mvec))
    if (length(commonM)) Msub[commonM] <- as.integer(Mvec[commonM])

    pv <- phased_list[[nm]]$phase_vec
    if (length(pv) != K)
      stop("phase_vec length for dam '", nm, "' must be length(order)-1.")

    G_list[[nm]] <- Gsub
    M_list[[nm]] <- Msub
    phase_list[[nm]] <- as.integer(pv)
  }

  # Resolve gametic pseudocount priors into per-dam engine priors (pi_prior_list)
  # and a shared total pseudocount (lambda). A single spec (numeric mean, or
  # list(alpha=, beta=)) is applied to all dams; a per-dam list supplies exactly
  # one spec for every requested dam. The total pseudocount (alpha + beta) must be
  # common across dams (the engine takes one lambda).
  Tn <- length(ord)
  if (is.null(q_prior_list)) {
    pi_prior_list_eff <- pi_prior_list
    lambda_eff <- lambda
  } else {
    is_single <- is.numeric(q_prior_list) ||
      (is.list(q_prior_list) && all(c("alpha", "beta") %in% names(q_prior_list)))
    if (is_single) {
      eng <- .hsmap_q_prior_to_engine(q_prior_list, lambda, Tn)
      pi_prior_list_eff <- stats::setNames(rep(list(eng$pi_prior), length(dams_req)), dams_req)
      lambda_eff <- eng$lambda
    } else {
      if (!is.list(q_prior_list))
        stop("`q_prior_list` must be NULL, a single spec (numeric mean or ",
             "list(alpha=, beta=)), or a per-dam list of such specs.", call. = FALSE)
      # per-dam list: exactly one spec per requested dam (no silent defaults)
      nm <- names(q_prior_list)
      if (!is.null(nm)) {
        miss <- setdiff(dams_req, nm)
        if (length(miss))
          stop("`q_prior_list` is a per-dam list but is missing prior(s) for dam(s): ",
               paste(miss, collapse = ", "),
               ". Provide one spec per requested dam, or a single shared spec.", call. = FALSE)
        unknown <- setdiff(nm, dams_req)
        if (length(unknown))
          stop("`q_prior_list` has unknown dam name(s): ", paste(unknown, collapse = ", "),
               ". Requested dams: ", paste(dams_req, collapse = ", "), ".", call. = FALSE)
        specs <- q_prior_list[dams_req]
      } else {
        if (length(q_prior_list) != length(dams_req))
          stop("unnamed `q_prior_list` must have exactly one spec per requested dam (",
               length(dams_req), "); got ", length(q_prior_list), ".", call. = FALSE)
        specs <- q_prior_list
      }
      if (any(vapply(specs, is.null, logical(1))))
        stop("`q_prior_list` must supply a non-NULL prior for every requested dam.", call. = FALSE)
      engs <- lapply(specs, function(spec) .hsmap_q_prior_to_engine(spec, lambda, Tn))
      lams <- vapply(engs, function(e) e$lambda, numeric(1))
      if (diff(range(lams)) > 1e-8 * max(c(lams, 1)))
        stop("joint fits require a common total pseudocount (alpha + beta) across dams.", call. = FALSE)
      pi_prior_list_eff <- stats::setNames(lapply(engs, function(e) e$pi_prior), dams_req)
      lambda_eff <- lams[[1]]
    }
    if (!is.null(pi_prior_list))
      warning("Both `q_prior_list` and `pi_prior_list` supplied; using `q_prior_list`.", call. = FALSE)
  }

  fit <- hmm_hs_joint_cpp(
    G_list = G_list,
    M_list = M_list,
    phase_list = phase_list,
    r_start = r_start,
    pi_mode = "HWE",
    pi_prior_list_in = pi_prior_list_eff,
    lambda = lambda_eff,
    epsilon = epsilon,
    tol = tol,
    maxit = maxit,
    paternal_mode = eff_paternal,
    Pi_prior_list_in = Pi_prior_list
  )

  # Canonical identifiable paternal output per dam: q_k = pi_AA + 0.5*pi_Aa.
  # fit$pi_list is retained only as the derived HWE-form emission table used by
  # downstream decoders; it is NOT an estimate of paternal genotype frequencies.
  fit$q_list <- stats::setNames(
    lapply(fit$pi_list, function(p)
      stats::setNames(as.numeric(p["AA", ] + 0.5 * p["Aa", ]), ord)),
    names(fit$pi_list))
  fit$paternal_mode <- paternal_mode

  out <- list(
    order      = ord,
    dams       = dams_req,
    phase_list = phase_list,
    fit        = fit,
    r          = fit$r
  )
  class(out) <- c("HSMap.map.joint", "HSMap.map")
  out
}
