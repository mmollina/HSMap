#' Blockwise multipoint fitting across unresolved phase
#'
#' @description
#' Fit the maternal multipoint HMM only within \strong{resolved phase blocks},
#' splitting the marker order wherever the relative phase of an adjacent interval is
#' unresolved (\code{NA} in \code{phase_vec}). Phase is never imputed across a
#' disconnected boundary, and no recombination estimate is produced across it.
#'
#' For a single dam the map is split wherever an adjacent interval is unresolved.
#' For several dams (an \code{HSMap.phased.multi}) a \emph{conservative common-block}
#' rule is used: an interval is a block boundary if it is unresolved for \strong{any}
#' dam that carries phase information adjacent to that interval (i.e. at least one of
#' the two markers is in a multi-marker phase component for that dam). Within a block,
#' the joint fit uses only the dams whose phase is fully resolved across that block;
#' the contributing dams are reported per block.
#'
#' @param x An \code{HSMap.data} object.
#' @param phased An \code{HSMap.phased} (single dam) or \code{HSMap.phased.multi}
#'   (several dams), typically from [phase_from_pairwise()].
#' @param ... Passed to [hmm_map()] (single dam) or [hmm_map_joint()] (several dams),
#'   e.g. \code{epsilon}, \code{paternal_mode}, \code{lambda}, \code{q_prior_in},
#'   \code{tol}, \code{maxit}.
#' @param min_block_markers Minimum markers for a block to be fitted (default 2).
#'   Smaller blocks are retained in metadata but not estimated.
#'
#' @return An \code{HSMap.map.blocks}: a list with \code{blocks} (one entry per block,
#' each with \code{markers}, \code{block}, \code{contributing_dams}, and the block
#' \code{fit} or \code{NULL}), \code{block_id} (per marker), \code{unresolved_boundaries}
#' (interval indices that separate blocks), \code{n_blocks}, and \code{interval_table}
#' (per adjacent interval: markers, block, status, and fitted r where available).
#'
#' @seealso [hmm_map()], [get_map()]
#' @export
hmm_map_blocks <- function(x, phased, ..., min_block_markers = 2L) {
  if (!inherits(x, "HSMap.data")) stop("`x` must be an HSMap.data object.")
  if (inherits(phased, "HSMap.phased.multi"))
    return(.hmm_map_blocks_joint(x, phased, ..., min_block_markers = min_block_markers))
  if (!inherits(phased, "HSMap.phased"))
    stop("`phased` must be an HSMap.phased or HSMap.phased.multi object.")

  ord <- as.character(phased$order)
  pv  <- phased$phase_vec
  if (length(pv) != length(ord) - 1L)
    stop("phase_vec length must be length(order) - 1.")
  ranges <- .blocks_from_phase(length(ord), pv)

  blocks <- vector("list", length(ranges))
  block_id <- integer(length(ord))
  for (b in seq_along(ranges)) {
    idx <- ranges[[b]]
    block_id[idx] <- b
    mk <- ord[idx]
    if (length(mk) < min_block_markers) {
      blocks[[b]] <- list(block = b, markers = mk, contributing_dams = phased$dam,
                          fit = NULL, note = "block with < min_block_markers; r not estimated")
      next
    }
    sub_pv <- as.integer(pv[idx[-length(idx)]])
    sub_ph <- structure(list(dam = phased$dam, order = mk, phase_vec = sub_pv),
                         class = "HSMap.phased")
    fit <- hmm_map(x, phased = sub_ph, ...)
    blocks[[b]] <- list(block = b, markers = mk, contributing_dams = phased$dam, fit = fit)
  }

  itab <- .block_interval_table(ord, pv, block_id, blocks)
  out <- list(
    blocks                = blocks,
    block_id              = stats::setNames(block_id, ord),
    unresolved_boundaries = which(is.na(pv)),
    n_blocks              = length(ranges),
    order                 = ord,
    interval_table        = itab,
    dams                  = phased$dam
  )
  class(out) <- "HSMap.map.blocks"
  out
}

# Joint (multi-dam) blockwise fitting with a conservative common-block rule.
.hmm_map_blocks_joint <- function(x, phased_multi, ..., min_block_markers = 2L) {
  dams <- names(phased_multi)
  if (is.null(dams)) dams <- paste0("Dam", seq_along(phased_multi))
  orders <- lapply(phased_multi, function(p) as.character(p$order))
  ord <- orders[[1]]
  if (!all(vapply(orders, function(o) identical(o, ord), logical(1))))
    stop("hmm_map_blocks() for several dams requires an identical marker `order` in ",
         "every dam (phase with a common order first).")
  Ti <- length(ord) - 1L

  # per-dam phase_vec and 'informative at interval t' (a marker of the interval sits
  # in a multi-marker phase component for that dam)
  pv_mat  <- vapply(phased_multi, function(p) as.integer(p$phase_vec), integer(Ti))
  inform  <- matrix(FALSE, Ti, length(dams))
  for (d in seq_along(dams)) {
    comp <- phased_multi[[d]]$component
    csz  <- table(comp)
    big  <- as.integer(names(csz)[csz > 1L])
    if (Ti >= 1) inform[, d] <- (comp[-length(comp)] %in% big) | (comp[-1] %in% big)
  }
  # boundary at t if any dam is informative-but-unresolved there
  boundary <- logical(Ti)
  for (t in seq_len(Ti)) boundary[t] <- any(inform[t, ] & is.na(pv_mat[t, ]))
  pv_common <- ifelse(boundary, NA_integer_, 0L)  # NA only marks boundaries for splitting
  ranges <- .blocks_from_phase(length(ord), pv_common)

  blocks <- vector("list", length(ranges))
  block_id <- integer(length(ord))
  for (b in seq_along(ranges)) {
    idx <- ranges[[b]]; block_id[idx] <- b; mk <- ord[idx]
    if (length(mk) < min_block_markers) {
      blocks[[b]] <- list(block = b, markers = mk, contributing_dams = character(0),
                          fit = NULL, note = "block with < min_block_markers; r not estimated")
      next
    }
    ivars <- idx[-length(idx)]                        # interval indices inside the block
    # a dam contributes if its phase is fully resolved across this block
    contrib <- dams[vapply(seq_along(dams), function(d) !anyNA(pv_mat[ivars, d]), logical(1))]
    if (!length(contrib)) {
      blocks[[b]] <- list(block = b, markers = mk, contributing_dams = character(0),
                          fit = NULL, note = "no dam fully resolved across this block")
      next
    }
    phlist <- lapply(contrib, function(dn)
      structure(list(dam = dn, order = mk, phase_vec = as.integer(pv_mat[ivars, dn])),
                class = "HSMap.phased"))
    names(phlist) <- contrib
    class(phlist) <- "HSMap.phased.multi"
    fit <- hmm_map_joint(x, phased = phlist, dam = contrib, ...)
    blocks[[b]] <- list(block = b, markers = mk, contributing_dams = contrib, fit = fit)
  }

  itab <- .block_interval_table(ord, pv_common, block_id, blocks)
  out <- list(
    blocks                = blocks,
    block_id              = stats::setNames(block_id, ord),
    unresolved_boundaries = which(boundary),
    n_blocks              = length(ranges),
    order                 = ord,
    interval_table        = itab,
    dams                  = dams
  )
  class(out) <- "HSMap.map.blocks"
  out
}

# Contiguous index ranges split at NA (unresolved) intervals.
.blocks_from_phase <- function(n_markers, phase_vec) {
  Ti <- length(phase_vec)
  if (n_markers <= 1L) return(list(seq_len(n_markers)))
  bnd <- which(is.na(phase_vec))
  starts <- c(1L, bnd + 1L)
  ends   <- c(bnd, n_markers)
  Map(function(s, e) s:e, starts, ends)
}

# One row per adjacent interval: markers, block (NA across a boundary), status, r.
.block_interval_table <- function(ord, phase_vec, block_id, blocks) {
  Ti <- length(ord) - 1L
  if (Ti < 1L)
    return(data.frame(from = character(0), to = character(0), block = integer(0),
                      status = character(0), r = numeric(0), stringsAsFactors = FALSE))
  status <- ifelse(is.na(phase_vec), "between_blocks", "linked")
  blk <- ifelse(is.na(phase_vec), NA_integer_, block_id[seq_len(Ti)])
  rr <- rep(NA_real_, Ti)
  for (b in seq_along(blocks)) {
    fit <- blocks[[b]]$fit
    if (is.null(fit)) next
    rv <- as.numeric(fit$fit$r)
    mk <- blocks[[b]]$markers
    # map block-internal intervals back to global interval indices
    gi <- match(mk[-length(mk)], ord)
    if (length(rv) == length(gi)) rr[gi] <- rv
  }
  data.frame(from = ord[-length(ord)], to = ord[-1], block = blk,
             status = status, r = rr, stringsAsFactors = FALSE)
}
