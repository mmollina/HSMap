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
#' @param gap_r No-linkage threshold (default \code{0.499}). A within-block interval
#'   whose fitted \code{r} is at or above \code{gap_r} is classified as a
#'   \code{no_linkage_boundary} (a gap), not \code{linked}.
#'
#' @return An \code{HSMap.map.blocks}: a list with \code{blocks} (one entry per block,
#' each with \code{markers}, \code{block}, \code{contributing_dams}, and the block
#' \code{fit} or \code{NULL}), \code{block_id} (per marker), \code{unresolved_boundaries}
#' (interval indices that separate blocks), \code{n_blocks}, \code{gap_r}, and
#' \code{interval_table} (per adjacent interval: markers, block, fitted r, and a
#' \code{status} that is one of \code{linked}, \code{no_linkage_boundary},
#' \code{unresolved_phase}, \code{between_blocks}, \code{insufficient_information}).
#' Use [get_block_map()] for safe positions/lengths and [plot_block_map()] to plot.
#'
#' @seealso [hmm_map()], [get_map()], [get_block_map()], [plot_block_map()]
#' @export
hmm_map_blocks <- function(x, phased, ..., min_block_markers = 2L, gap_r = 0.499) {
  if (!inherits(x, "HSMap.data")) stop("`x` must be an HSMap.data object.")
  if (!is.numeric(gap_r) || length(gap_r) != 1L || !is.finite(gap_r) || gap_r <= 0 || gap_r > 0.5)
    stop("`gap_r` must be a single number in (0, 0.5].")
  if (inherits(phased, "HSMap.phased.multi"))
    return(.hmm_map_blocks_joint(x, phased, ..., min_block_markers = min_block_markers, gap_r = gap_r))
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

  itab <- .block_interval_table(ord, pv, block_id, blocks,
                                gap_r = gap_r, boundary_status = "unresolved_phase")
  out <- list(
    blocks                = blocks,
    block_id              = stats::setNames(block_id, ord),
    unresolved_boundaries = which(is.na(pv)),
    n_blocks              = length(ranges),
    order                 = ord,
    interval_table        = itab,
    gap_r                 = gap_r,
    dams                  = phased$dam
  )
  class(out) <- "HSMap.map.blocks"
  out
}

# Joint (multi-dam) blockwise fitting with a conservative common-block rule.
.hmm_map_blocks_joint <- function(x, phased_multi, ..., min_block_markers = 2L, gap_r = 0.499) {
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

  itab <- .block_interval_table(ord, pv_common, block_id, blocks,
                                gap_r = gap_r, boundary_status = "between_blocks")
  out <- list(
    blocks                = blocks,
    block_id              = stats::setNames(block_id, ord),
    unresolved_boundaries = which(boundary),
    n_blocks              = length(ranges),
    order                 = ord,
    interval_table        = itab,
    gap_r                 = gap_r,
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

# One row per adjacent interval, classified into one of:
#   linked                    within-block, fitted, r < gap_r
#   no_linkage_boundary       within-block, fitted, r at/above gap_r
#   insufficient_information  within-block but not fitted (block < 2 markers / no r)
#   between_blocks            boundary between blocks (joint conservative common rule)
#   unresolved_phase          boundary due to unresolved relative phase (single-dam)
# A within-block interval is NOT called linked merely because it was fitted; a fitted
# r at/near 0.5 is a no_linkage_boundary.
.block_interval_table <- function(ord, phase_vec, block_id, blocks,
                                  gap_r = 0.499, boundary_status = "unresolved_phase") {
  Ti <- length(ord) - 1L
  if (Ti < 1L)
    return(data.frame(from = character(0), to = character(0), block = integer(0),
                      status = character(0), r = numeric(0), stringsAsFactors = FALSE))
  rr  <- rep(NA_real_, Ti)
  for (b in seq_along(blocks)) {
    fit <- blocks[[b]]$fit
    if (is.null(fit)) next
    rv <- as.numeric(fit$fit$r); mk <- blocks[[b]]$markers
    gi <- match(mk[-length(mk)], ord)
    if (length(rv) == length(gi)) rr[gi] <- rv
  }
  blk <- ifelse(is.na(phase_vec), NA_integer_, block_id[seq_len(Ti)])
  status <- character(Ti)
  for (t in seq_len(Ti)) {
    if (is.na(phase_vec[t]))            status[t] <- boundary_status          # split boundary
    else if (is.na(rr[t]))             status[t] <- "insufficient_information"
    else if (rr[t] >= gap_r)           status[t] <- "no_linkage_boundary"    # fitted r ~ 0.5
    else                               status[t] <- "linked"
  }
  data.frame(from = ord[-length(ord)], to = ord[-1], block = blk,
             status = status, r = rr, stringsAsFactors = FALSE)
}


#' Safe map summary for a blockwise map
#'
#' @description
#' Extract marker positions, interval distances and statuses, per-block lengths, the
#' total finite linked length, and the gaps from an [hmm_map_blocks()] result --
#' without ever turning a no-linkage or unresolved interval into a large finite
#' distance. Positions are cumulative \emph{within} each block and reset between
#' blocks; gap intervals have \code{NA} distance.
#'
#' @param x An \code{HSMap.map.blocks} object.
#' @param map.function Distance function, one of \code{"haldane"}, \code{"kosambi"},
#'   \code{"morgan"}.
#' @param gap_r Optional override of the no-linkage threshold; defaults to the value
#'   stored on \code{x}. Intervals with fitted \code{r >= gap_r} are gaps.
#'
#' @return A list with: \code{positions} (data frame: \code{marker}, \code{block},
#'   within-block \code{pos}); \code{interval_table} (\code{from}, \code{to},
#'   \code{block}, \code{status}, \code{r}, \code{dist_cM}); \code{block_lengths};
#'   \code{total_linked_length}; \code{n_gaps}; \code{gap_intervals};
#'   \code{map.function}; \code{gap_r}.
#' @seealso [hmm_map_blocks()], [get_map()], [plot_block_map()]
#' @export
get_block_map <- function(x, map.function = c("haldane", "kosambi", "morgan"),
                          gap_r = NULL) {
  if (!inherits(x, "HSMap.map.blocks")) stop("`x` must be an HSMap.map.blocks object.")
  map.function <- match.arg(map.function)
  gr <- if (is.null(gap_r)) (x$gap_r %||% 0.499) else gap_r
  mf <- switch(map.function, haldane = inv_haldane, kosambi = inv_kosambi, morgan = inv_morgan)

  positions <- list(); block_lengths <- numeric(length(x$blocks))
  for (b in seq_along(x$blocks)) {
    blk <- x$blocks[[b]]
    if (is.null(blk$fit) || length(blk$markers) < 2L) {
      block_lengths[b] <- 0
      positions[[b]] <- data.frame(marker = blk$markers, block = b,
                                   pos = if (length(blk$markers)) rep(0, length(blk$markers)) else numeric(0),
                                   stringsAsFactors = FALSE)
      next
    }
    pos <- get_map(blk$fit, map.function, gap_r = gr)
    positions[[b]] <- data.frame(marker = names(pos), block = b, pos = as.numeric(pos),
                                 stringsAsFactors = FALSE)
    block_lengths[b] <- attr(pos, "total_linked_length")
  }
  pos_df <- if (length(positions)) do.call(rbind, positions) else
    data.frame(marker = character(0), block = integer(0), pos = numeric(0))

  itab <- x$interval_table
  itab$dist_cM <- ifelse(itab$status == "linked" & is.finite(itab$r), mf(itab$r), NA_real_)

  list(
    positions           = pos_df,
    interval_table      = itab,
    block_lengths       = block_lengths,
    total_linked_length = sum(block_lengths),
    n_gaps              = sum(itab$status != "linked"),
    gap_intervals       = which(itab$status != "linked"),
    map.function        = map.function,
    gap_r               = gr
  )
}


#' Block-aware map plot
#'
#' @description
#' Plot a blockwise map with each resolved phase block as its own panel, so
#' reset within-block positions are never drawn on top of one another as a single
#' continuous chromosome.
#'
#' @param x An \code{HSMap.map.blocks} object.
#' @param map.function Distance function; see [get_block_map()].
#' @return A \pkg{ggplot2} object (one panel per block).
#' @seealso [get_block_map()]
#' @export
plot_block_map <- function(x, map.function = c("haldane", "kosambi", "morgan")) {
  map.function <- match.arg(map.function)
  bm <- get_block_map(x, map.function)
  df <- bm$positions
  if (!nrow(df)) stop("No markers to plot.")
  df$block <- factor(df$block)
  ggplot2::ggplot(df, ggplot2::aes(x = pos, y = 0)) +
    ggplot2::geom_line(ggplot2::aes(group = block), colour = "grey60") +
    ggplot2::geom_point(size = 1.6, colour = "#00AFBB") +
    ggplot2::facet_wrap(~ block, scales = "free_x", ncol = 1,
                        labeller = ggplot2::labeller(.default = function(v) paste0("block ", v))) +
    ggplot2::labs(x = sprintf("Within-block position (cM, %s)", bm$map.function), y = NULL,
                  title = sprintf("Blockwise map: %d block(s), %d gap(s), total linked %.1f cM",
                                  length(bm$block_lengths), bm$n_gaps, bm$total_linked_length)) +
    ggplot2::theme_minimal(11) +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank())
}
