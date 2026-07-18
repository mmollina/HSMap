#' Phase from pairwise for one or more dams (HSMap)
#'
#' @description
#' Convert dam-specific pairwise maternal phase calls and their LOD support into an
#' ordered phase assignment along a user-supplied marker order, honestly separating
#' resolved from unresolved phase.
#'
#' The signed, LOD-weighted phase graph is built **separately for each dam** from
#' that dam's own `lod_ph_list[[d]]` (not the pooled `lod_ph`). Only edges with an
#' available phase sign, finite LOD, and LOD above the support threshold are used.
#' The graph is split into connected components; within each component a spectral
#' initialization and greedy coordinate ascent orient the markers (a heuristic, not
#' guaranteed globally optimal). The relative orientation of different components is
#' **unidentified** and is reported as unresolved rather than forced to coupling.
#'
#' @param tpt An object of class \code{HSMap.tpt}. Must contain
#'   \code{tpt$fit$mom_phase_list} (per-dam 0/1 phase matrices) and per-dam LOD
#'   support in \code{tpt$fit$lod_ph_list} (falling back to the pooled
#'   \code{tpt$fit$lod_ph} with a warning if the per-dam list is absent).
#' @param order Character vector with marker IDs in the desired linear order.
#'   Markers absent for a dam are dropped with a warning.
#' @param dam Which dam(s): a dam name/index, a vector of them, or \code{"all"}
#'   (default).
#' @param min_phase_lod Minimum phase LOD for an edge to be treated as supported
#'   (default \code{0}). Even at the default, edges with zero/near-zero evidence are
#'   excluded (see \code{tie_tol}); an edge counts only when its LOD strictly exceeds
#'   \code{max(min_phase_lod, tie_tol)}. \strong{The default \code{0} is a backward-
#'   compatibility choice, not a statistically validated threshold}; a suitable value
#'   awaits a simulation study.
#' @param tie_tol Numerical tolerance below which a coupling-vs-repulsion LOD is
#'   treated as an unresolved tie (default \code{1e-8}).
#' @param anchor_idx Integer index in \code{order} whose component is sign-anchored
#'   (breaks that component's sign symmetry only). Default \code{1}.
#' @param anchor_label Integer in \code{c(0,1)} for the anchor marker. Default \code{0}.
#' @param max_passes Maximum greedy passes per component. Default \code{50}.
#' @param tol Non-negative improvement tolerance for greedy ascent. Default \code{1e-9}.
#' @param verbose Logical; print progress. Default \code{FALSE}.
#'
#' @return For a single dam, an \code{HSMap.phased} with (among others): \code{order};
#' \code{clusters} (0/1 homolog labels, meaningful only within a component);
#' \code{component} (component ID per marker); \code{phase_vec} (length
#' \code{T-1}; 1=coupling, 0=repulsion for resolved adjacent intervals, \code{NA}
#' for unresolved ones); \code{resolved_interval}; \code{direct_lod} (the \emph{raw}
#' finite LOD of the adjacent \eqn{(t, t+1)} edge, \code{NA} if that edge is absent or
#' non-finite -- reported even when it is below threshold, and NOT to be read as the
#' total support for a path-derived phase); \code{direct_supported} (whether that
#' adjacent edge passes the phase-present, finite-value, tie, and threshold rules);
#' \code{resolved_via} (\code{"direct"} if the phase comes from the supported adjacent
#' edge, \code{"path"} if it is derived through the component -- possibly with zero
#' direct adjacent LOD -- or \code{"unresolved"}); \code{interval_support} (a
#' compatibility alias of \code{direct_lod});
#' \code{n_components}, \code{component_sizes}; \code{unresolved_markers},
#' \code{unresolved_intervals}; \code{objective}, \code{component_objective},
#' \code{converged}, \code{n_flips}. Multiple dams return an
#' \code{HSMap.phased.multi} (a named list).
#'
#' @export
phase_from_pairwise <- function(
    tpt,
    order,
    dam           = "all",
    min_phase_lod = 0,
    tie_tol       = 1e-8,
    anchor_idx    = 1L,
    anchor_label  = 0L,
    max_passes    = 50L,
    tol           = 1e-9,
    verbose       = FALSE
) {
  if (!inherits(tpt, "HSMap.tpt")) stop("`tpt` must be HSMap.tpt.")
  if (!is.character(order) || !length(order))
    stop("`order` must be a non-empty character vector of marker IDs.")
  if (!is.numeric(min_phase_lod) || length(min_phase_lod) != 1L ||
      !is.finite(min_phase_lod) || min_phase_lod < 0)
    stop("`min_phase_lod` must be a single finite, non-negative number.")

  mom_list <- tpt$fit$mom_phase_list
  if (!is.list(mom_list) || !length(mom_list))
    stop("tpt$fit$mom_phase_list is missing or empty.")
  dam_names <- names(mom_list)
  if (is.null(dam_names)) dam_names <- paste0("Dam", seq_along(mom_list))

  # Per-dam LOD support (preferred); fall back to pooled lod_ph with a warning.
  lod_list <- tpt$fit$lod_ph_list
  use_pooled <- FALSE
  if (!is.list(lod_list) || !length(lod_list)) {
    if (!is.matrix(tpt$fit$lod_ph))
      stop("tpt$fit needs `lod_ph_list` (per dam) or a pooled `lod_ph` matrix.")
    warning("tpt$fit$lod_ph_list absent; falling back to the pooled lod_ph for every ",
            "dam. Per-dam edge weights are preferred.", call. = FALSE)
    use_pooled <- TRUE
  }

  sel <- if (identical(dam, "all")) dam_names
    else if (is.character(dam)) {
      if (!all(dam %in% dam_names)) stop("Unknown dam(s): ", paste(setdiff(dam, dam_names), collapse = ", "))
      dam
    } else if (is.numeric(dam)) {
      idx <- as.integer(dam)
      if (any(idx < 1L | idx > length(dam_names))) stop("Dam index out of range.")
      dam_names[idx]
    } else stop("`dam` must be 'all', dam name(s), or dam index/indices.")

  thr <- max(min_phase_lod, tie_tol)

  phase_one <- function(dam_label) {
    A_all <- mom_list[[dam_label]]
    if (!is.matrix(A_all)) stop("Phase matrix for dam '", dam_label, "' is not a matrix.")
    L_all <- if (use_pooled) tpt$fit$lod_ph else lod_list[[dam_label]]
    if (!is.matrix(L_all)) stop("LOD matrix for dam '", dam_label, "' is not a matrix.")

    common <- intersect(rownames(A_all), rownames(L_all))
    if (length(common) < 2L) stop("Not enough common markers for dam '", dam_label, "'.")
    o <- intersect(order, common)
    if (!length(o)) stop("None of the SNPs in `order` are present for dam '", dam_label, "'.")
    if (length(o) < length(order)) {
      miss <- setdiff(order, o)
      warning("Dropped ", length(miss), " missing SNP(s) for dam '", dam_label,
              "'. Example: ", paste(utils::head(miss, 5), collapse = ", "), call. = FALSE)
    }
    A   <- A_all[o, o, drop = FALSE]
    LOD <- L_all[o, o, drop = FALSE]
    Tn  <- length(o)

    # Supported edges: phase sign present, finite LOD strictly above threshold.
    Sup <- (!is.na(A)) & is.finite(LOD) & (LOD > thr)
    diag(Sup) <- FALSE
    Sup <- Sup | t(Sup)

    # Signed, LOD-weighted graph on supported edges only.
    S <- matrix(0, Tn, Tn)
    S[!is.na(A) & A >= 0.5] <-  1
    S[!is.na(A) & A <  0.5] <- -1
    W <- LOD; W[!is.finite(W)] <- 0
    S[!Sup] <- 0; W[!Sup] <- 0
    diag(S) <- 0; diag(W) <- 0
    J <- (S * W); J <- (J + t(J)) / 2

    comp <- .pf_components(Sup)
    x <- rep(1, Tn)
    comp_obj <- numeric(0); comp_conv <- logical(0); tot_flips <- 0L; tot_iters <- 0L

    for (cc in sort(unique(comp))) {
      idx <- which(comp == cc)
      if (length(idx) == 1L) { x[idx] <- 1; next }
      Jc <- J[idx, idx, drop = FALSE]
      if (all(Jc == 0)) { x[idx] <- 1; next }
      v1 <- .pf_leading_eigvec(Jc)
      xc <- ifelse(v1 >= 0, 1, -1)
      gr <- pf_greedy_cpp(Jc, as.integer(xc), max_passes, tol)
      xc <- gr$x
      if (xc[1] != 1) xc <- -xc                       # anchor component's first marker to +1
      x[idx] <- xc
      comp_obj[as.character(cc)]  <- gr$objective
      comp_conv[as.character(cc)] <- gr$converged
      tot_flips <- tot_flips + gr$n_flips
      tot_iters <- max(tot_iters, gr$iters)
      if (isTRUE(verbose))
        message("dam=", dam_label, " comp=", cc, " size=", length(idx),
                " conv=", gr$converged, " flips=", gr$n_flips)
    }

    # Global anchor: flip only the anchor's component so order[anchor_idx] == anchor_label.
    if (anchor_idx >= 1L && anchor_idx <= Tn) {
      want_plus <- (anchor_label == 0L)
      if ((x[anchor_idx] == 1) != want_plus) {
        ai <- which(comp == comp[anchor_idx]); x[ai] <- -x[ai]
      }
    }
    clusters <- ifelse(x == 1, 0L, 1L)

    # Intervals: resolved iff both markers share a component (path of supported edges).
    resolved <- if (Tn >= 2) comp[-Tn] == comp[-1] else logical(0)
    phase_vec <- if (Tn >= 2)
      ifelse(resolved, as.integer(clusters[-Tn] == clusters[-1]), NA_integer_) else integer(0)
    # Direct adjacent-edge phase LOD (the support of the (t, t+1) edge ITSELF), and
    # whether that adjacent edge is a supported edge. An interval can be resolved
    # INDIRECTLY (both markers in one component, via a longer path) even when the
    # direct adjacent LOD is 0/absent: `resolved_via` distinguishes the two. Do not
    # read a zero direct adjacent LOD as the total support for a path-resolved interval.
    idx2 <- if (Tn >= 2) cbind(seq_len(Tn - 1L), seq_len(Tn - 1L) + 1L) else NULL
    # direct_lod: the RAW LOD of the adjacent (t, t+1) edge (NA if non-finite/absent),
    # reported even when below threshold. direct_supported: whether that edge passes
    # the phase-present, finite-value, tie, and threshold rules. resolved_via records
    # whether the interval is resolved by that direct edge, indirectly via a path in
    # the component, or not at all.
    direct_lod <- if (Tn >= 2) { s <- LOD[idx2]; s[!is.finite(s)] <- NA_real_; s } else numeric(0)
    direct_supported <- if (Tn >= 2) as.logical(Sup[idx2]) else logical(0)
    resolved_via <- if (Tn >= 2)
      ifelse(!resolved, "unresolved", ifelse(direct_supported, "direct", "path")) else character(0)

    csz <- as.integer(table(comp))
    unresolved_markers   <- o[comp %in% which(csz == 1L)]
    unresolved_intervals <- if (Tn >= 2) which(!resolved) else integer(0)

    if (!any(Sup))
      warning("No supported phase edges for dam '", dam_label,
              "': every marker is its own component and all intervals are unresolved. ",
              "Phase is NOT resolved (not returned as all-coupling).", call. = FALSE)

    agree <- (S == 1 & outer(clusters, clusters, "==")) |
             (S == -1 & outer(clusters, clusters, "!="))
    num <- sum(W[upper.tri(W)] * agree[upper.tri(agree)])
    den <- sum(W[upper.tri(W)])
    w_agree <- if (den > 0) num / den else NA_real_

    out <- list(
      dam                = dam_label,
      order              = o,
      clusters           = as.integer(clusters),
      component          = as.integer(comp),
      phase_vec          = as.integer(phase_vec),
      resolved_interval  = resolved,
      direct_lod         = direct_lod,        # raw finite LOD of the adjacent (t,t+1) edge (NA if absent)
      direct_supported   = direct_supported,  # edge passes phase/finite/tie/threshold rules
      resolved_via       = resolved_via,       # "direct" | "path" | "unresolved"
      interval_support   = direct_lod,         # compatibility alias of direct_lod (raw adjacent LOD)
      n_components       = length(unique(comp)),
      component_sizes    = csz,
      unresolved_markers = unresolved_markers,
      unresolved_intervals = as.integer(unresolved_intervals),
      objective          = if (length(comp_obj)) sum(comp_obj) else 0,
      component_objective = comp_obj,
      weighted_agreement = w_agree,
      iters              = tot_iters,
      n_flips            = tot_flips,
      converged          = if (length(comp_conv)) all(comp_conv) else TRUE,
      min_phase_lod      = min_phase_lod,
      fully_resolved     = (length(unique(comp)) == 1L)
    )
    class(out) <- "HSMap.phased"
    out
  }

  res <- lapply(sel, phase_one)
  names(res) <- sel
  if (length(res) == 1L) return(res[[1L]])
  class(res) <- "HSMap.phased.multi"
  res
}

# ---- internal helpers -------------------------------------------------------

# Connected components of a symmetric logical adjacency matrix (BFS/DFS).
.pf_components <- function(Sup) {
  n <- nrow(Sup); comp <- integer(n); cid <- 0L
  for (s in seq_len(n)) {
    if (comp[s] != 0L) next
    cid <- cid + 1L
    stack <- s
    while (length(stack)) {
      v <- stack[length(stack)]; stack <- stack[-length(stack)]
      if (comp[v] != 0L) next
      comp[v] <- cid
      nb <- which(Sup[v, ] & comp == 0L)
      if (length(nb)) stack <- c(stack, nb)
    }
  }
  comp
}

# Leading (largest-algebraic) eigenvector of a component's signed graph, used only to
# seed the greedy ascent. Only this one eigenvector is needed, so for components large
# enough to benefit we use RSpectra::eigs_sym() (a partial solver) instead of a full
# eigen() decomposition. Base eigen() is the fallback for small components, when
# RSpectra is unavailable, or if the iterative solver fails. The sign is fixed
# deterministically (largest-magnitude entry made non-negative); note that a global sign
# flip of a whole component is an equivalent phase solution and leaves phase_vec
# unchanged, so this only pins the raw cluster labels for reproducibility.
.pf_leading_eigvec <- function(Jc) {
  n <- nrow(Jc)
  v <- NULL
  if (n >= 8L && requireNamespace("RSpectra", quietly = TRUE)) {
    v <- tryCatch({
      e <- RSpectra::eigs_sym(Jc, k = 1L, which = "LA")
      vv <- if (!is.null(e$vectors) && NCOL(e$vectors) >= 1L) e$vectors[, 1] else NULL
      if (is.null(vv) || length(vv) != n || anyNA(vv)) NULL else vv
    }, error = function(err) NULL)
  }
  if (is.null(v)) v <- eigen(Jc, symmetric = TRUE, only.values = FALSE)$vectors[, 1]
  k <- which.max(abs(v))
  if (length(k) && is.finite(v[k]) && v[k] < 0) v <- -v
  v
}

# Greedy coordinate ascent on the signed quadratic sum_{i<j} J[i,j] x_i x_j (x in +/-1).
# R REFERENCE implementation, retained for tests only; the normal workflow uses the
# equivalent C++ pf_greedy_cpp() (same objective and best-improvement rule).
.pf_greedy <- function(J, x, max_passes, tol) {
  Tn <- length(x); improve <- TRUE; pass <- 0L; nflips <- 0L
  while (improve && pass < max_passes) {
    pass <- pass + 1L; improve <- FALSE
    g <- as.numeric(J %*% x)
    for (iter in seq_len(5L * Tn)) {
      delta <- -2 * x * g
      best  <- which.max(delta)
      if (!length(best) || delta[best] <= tol) break
      improve <- TRUE; nflips <- nflips + 1L
      xo <- x[best]; x[best] <- -xo; g <- g - 2 * xo * J[, best]
    }
  }
  obj <- if (Tn >= 2) sum(J[upper.tri(J)] * outer(x, x)[upper.tri(J)]) else 0
  list(x = x, objective = obj, converged = !improve, iters = pass, n_flips = nflips)
}
