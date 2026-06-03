#' Phase from pairwise for one or more dams (HSMap)
#'
#' @description
#' Convert pairwise maternal phase calls and their LOD support into a robust,
#' ordered phase string along a user-supplied marker order. Supports one dam
#' or more dams in a single call. Internally it builds a signed, LOD-weighted
#' pairwise matrix and optimizes a 2-cluster assignment (homolog 0 vs 1)
#' with a spectral initialization followed by greedy coordinate ascent.
#'
#' @param tpt An object of class \code{HSMap.tpt}. Must contain
#'   \code{tpt$fit$mom_phase_list} (a list of per-dam 0/1 matrices:
#'   1=coupling, 0=repulsion) and \code{tpt$fit$lod_ph} (non-negative
#'   LOD support matrix), with identical dimnames.
#' @param order Character vector with SNP marker IDs in the desired linear order.
#'   Only markers present in \code{tpt} are used; missing names are dropped with
#'   a warning.
#' @param dam Which dam(s) to phase. One of:
#'   \itemize{
#'     \item a single dam name (matching \code{names(tpt$fit$mom_phase_list)})
#'           or integer index,
#'     \item a vector of dam names or indices,
#'     \item \code{"all"} to process every dam available.
#'   }
#'   Default \code{"all"}.
#' @param anchor_idx Integer index within \code{order} to fix an absolute label
#'   (breaks the global sign symmetry). Default \code{1}.
#' @param anchor_label Integer in \code{c(0,1)} to fix the label at
#'   \code{order[anchor_idx]} (0 means homolog 0, 1 means homolog 1). Default \code{0}.
#' @param max_passes Maximum greedy passes. Default \code{50}.
#' @param tol Non-negative improvement tolerance. Default \code{1e-9}.
#' @param verbose Logical; print progress. Default \code{FALSE}.
#'
#' @details
#' Let \code{A} be the dam-specific pairwise phase matrix (1=coupling, 0=repulsion,
#' NA uninformative) and \code{W} be the \code{lod_ph} matrix (non-negative). The
#' algorithm builds \code{S = sign(2*A - 1)} (so coupling=+1, repulsion=-1, NA=0)
#' and sets diagonal to 0. It then forms a symmetric signed weight matrix
#' \code{J = (S * W)}, symmetrized. A spectral initialization uses the leading
#' eigenvector of \code{J} to seed a +/-1 assignment; the assignment is flipped
#' so that \code{order[anchor_idx]} matches \code{anchor_label}. A greedy
#' coordinate ascent flips markers that improve the quadratic objective
#' \code{sum_{i<j} J\[i,j\] * x\[i\] * x\[j\]} until no flip exceeds \code{tol} or
#' \code{max_passes} passes are reached. The final 0/1 clusters define the
#' adjacent phase vector: phase\[t\]=1 (coupling) when clusters\[t\]==clusters\[t+1\],
#' otherwise 0 (repulsion).
#'
#' @return
#' If a single dam is requested, returns an \code{HSMap.phased}:
#' \itemize{
#'   \item \code{dam} Dam label.
#'   \item \code{order} Marker order used (after dropping missing).
#'   \item \code{clusters} Integer vector in \code{c(0,1)} (homolog assignment).
#'   \item \code{phase_vec} Integer vector length \code{length(order)-1}, 1=coupling, 0=repulsion.
#'   \item \code{objective}, \code{weighted_agreement}, \code{iters}, \code{converged}.
#' }
#' If multiple dams are requested (including \code{"all"}), returns an
#' \code{HSMap.phased.multi}: a named list of \code{HSMap.phased} objects (one per dam).
#'
#' @export
phase_from_pairwise <- function(
    tpt,
    order,
    dam = "all",
    anchor_idx   = 1L,
    anchor_label = 0L,
    max_passes   = 50L,
    tol          = 1e-9,
    verbose      = FALSE
) {
  if (!inherits(tpt, "HSMap.tpt")) stop("`tpt` must be HSMap.tpt.")
  if (!is.character(order) || !length(order))
    stop("`order` must be a non-empty character vector of marker IDs.")

  # available dams
  mom_list <- tpt$fit$mom_phase_list
  if (!is.list(mom_list) || !length(mom_list))
    stop("tpt$fit$mom_phase_list is missing or empty.")
  dam_names <- names(mom_list)
  if (is.null(dam_names)) dam_names <- paste0("Dam", seq_along(mom_list))

  # resolve selection
  sel <- if (identical(dam, "all")) {
    dam_names
  } else if (is.character(dam)) {
    if (!all(dam %in% dam_names))
      stop("Unknown dam(s): ", paste(setdiff(dam, dam_names), collapse = ", "))
    dam
  } else if (is.numeric(dam)) {
    idx <- as.integer(dam)
    if (any(idx < 1L | idx > length(dam_names)))
      stop("Dam index out of range.")
    dam_names[idx]
  } else {
    stop("`dam` must be 'all', dam name(s), or dam index/indices.")
  }

  # common LOD (shared across dams in current HSMap.tpt)
  LOD_all <- tpt$fit$lod_ph
  if (!is.matrix(LOD_all)) stop("tpt$fit$lod_ph must be a matrix.")

  # align to common names then subset to order
  # helper to phase one dam
  phase_one <- function(dam_label) {
    A_all <- mom_list[[dam_label]]
    if (!is.matrix(A_all)) stop("Phase matrix for dam '", dam_label, "' is not a matrix.")

    common <- intersect(rownames(A_all), rownames(LOD_all))
    if (length(common) < 2L) stop("Not enough common markers for dam '", dam_label, "'.")
    A0   <- A_all[common, common, drop = FALSE]
    LOD0 <- LOD_all[common, common, drop = FALSE]

    o <- intersect(order, common)
    if (!length(o)) stop("None of the SNPs in `order` are present for dam '", dam_label, "'.")
    if (length(o) < length(order)) {
      missing <- setdiff(order, o)
      warning("Dropped ", length(missing), " missing SNP(s) for dam '", dam_label,
              "'. Example: ", paste(utils::head(missing, 5), collapse = ", "))
    }
    A   <- A0[o, o, drop = FALSE]
    LOD <- LOD0[o, o, drop = FALSE]

    Tn <- nrow(A)
    if (Tn < 2) {
      out <- list(dam = dam_label, order = o, clusters = integer(Tn),
                  phase_vec = if (Tn) integer(Tn - 1L) else integer(0),
                  objective = 0, weighted_agreement = NA_real_,
                  iters = 0L, converged = TRUE)
      class(out) <- "HSMap.phased"
      return(out)
    }

    # Signed matrix S from pairwise calls, weighted by LOD
    S <- matrix(0, Tn, Tn)
    S[A >= 0.5 & !is.na(A)] <-  1
    S[A <  0.5 & !is.na(A)] <- -1
    W <- LOD; W[is.na(W)] <- 0
    diag(S) <- 0; diag(W) <- 0
    S <- (S + t(S)) / 2
    W <- (W + t(W)) / 2
    J <- (S * W + t(S * W)) / 2

    # Degenerate: no information
    if (all(J == 0)) {
      clusters  <- integer(Tn)
      phase_vec <- as.integer(rep(1L, Tn - 1L))
      out <- list(dam = dam_label, order = o, clusters = clusters, phase_vec = phase_vec,
                  objective = 0, weighted_agreement = NA_real_,
                  iters = 0L, converged = TRUE)
      class(out) <- "HSMap.phased"
      return(out)
    }

    # Spectral init
    ev <- eigen(J, symmetric = TRUE, only.values = FALSE)
    v1 <- ev$vectors[, 1]
    x  <- ifelse(v1 >= 0, 1, -1)

    # Anchor orientation
    if (anchor_idx < 1L || anchor_idx > Tn)
      stop("`anchor_idx` must be within 1..length(order).")
    want_plus <- (anchor_label == 0L)
    if ((x[anchor_idx] == 1) != want_plus) x <- -x

    # Greedy coordinate ascent
    improve <- TRUE; pass <- 0L
    while (improve && pass < max_passes) {
      pass <- pass + 1L
      improve <- FALSE
      g <- as.numeric(J %*% x)
      for (iter in seq_len(5L * Tn)) {
        delta <- -2 * x * g
        best  <- which.max(delta)
        if (length(best) == 0 || delta[best] <= tol) break
        improve <- TRUE
        x_old <- x[best]
        x[best] <- -x_old
        g <- g - 2 * x_old * J[, best]
      }
      if (isTRUE(verbose)) message("dam=", dam_label, " pass=", pass, " improved=", improve)
    }

    clusters  <- ifelse(x == 1, 0L, 1L)
    phase_vec <- as.integer(clusters[-Tn] == clusters[-1])

    # Agreement diagnostic
    agree <- (S ==  1 & outer(clusters, clusters, "==")) |
      (S == -1 & outer(clusters, clusters, "!="))
    agree[is.na(agree)] <- TRUE
    num <- sum(W[upper.tri(W)] * agree[upper.tri(agree)])
    den <- sum(W[upper.tri(W)])
    w_agree <- if (den > 0) num / den else NA_real_

    obj <- sum(J[upper.tri(J)] * (outer(x, x, "*")[upper.tri(J)]))

    out <- list(
      dam        = dam_label,
      order      = o,
      clusters   = as.integer(clusters),
      phase_vec  = as.integer(phase_vec),
      objective  = obj,
      weighted_agreement = w_agree,
      iters      = pass,
      converged  = !improve
    )
    class(out) <- "HSMap.phased"
    out
  }

  res <- lapply(sel, phase_one)
  names(res) <- sel

  if (length(res) == 1L) {
    return(res[[1L]])
  } else {
    class(res) <- "HSMap.phased.multi"
    return(res)
  }
}
