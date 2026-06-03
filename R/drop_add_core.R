#' Remove SNPs that create large adjacent gaps in an HSMap map
#'
#' Detect adjacent intervals whose Haldane distance exceeds a threshold
#' and remove either one marker per gap, or entire small clusters between
#' two large gaps. The procedure can iterate until no gap exceeds the
#' threshold. The remaining SNP names are returned in map order.
#'
#' When \code{cluster_size > 0}, any contiguous block of markers that lies
#' between two large gaps and has size less than or equal to \code{cluster_size}
#' is removed as a whole in a single pass. Edge blocks are not considered
#' clusters because they are not bounded by large gaps on both sides.
#'
#' The function accepts:
#' \itemize{
#'   \item a single-dam map of class \code{"HSMap.map"} with fields
#'         \code{$order} (marker names in map order) and \code{$fit$r}
#'         (adjacent recombination fractions), or
#'   \item a multi-dam map of class \code{"HSMap.map.multi"} with a
#'         list \code{$per_dam} containing one \code{"HSMap.map"} per dam.
#' }
#'
#' Distances in cM are computed from adjacent recombination fractions using
#' the Haldane mapping function: \code{d_cm = -50 * log(1 - 2r)} after clamping
#' \code{r} into \code{[1e-8, 0.499999]} to avoid numerical issues.
#'
#' @param x A single-dam \code{"HSMap.map"} or a multi-dam
#'   \code{"HSMap.map.multi"} with element \code{$per_dam}.
#' @param gap_cm Numeric threshold in centiMorgans. Any adjacent interval with
#'   distance greater than \code{gap_cm} is treated as a large gap. Default 10.
#' @param side Which endpoint to remove for each large gap when
#'   \code{cluster_size == 0}. One of:
#'   \describe{
#'     \item{\code{"auto"}}{remove the endpoint that belongs to the smaller
#'           neighboring block. Ties drop the right endpoint.}
#'     \item{\code{"right"}}{drop the right endpoint of each large gap.}
#'     \item{\code{"left"}}{drop the left endpoint of each large gap.}
#'     \item{\code{"both"}}{drop both endpoints of each large gap.}
#'   }
#'   Default \code{"auto"}.
#' @param cluster_size Integer. If greater than 0, remove entire blocks of
#'   consecutive markers that are bounded by two large gaps and whose size is
#'   less than or equal to \code{cluster_size}. If 0, remove markers one by one.
#'   Default 0.
#' @param iterate Logical. If \code{TRUE}, repeat passes until no large gap
#'   remains among the currently kept markers. Default \code{TRUE}.
#' @param dam For multi-dam maps, which dam(s) to process. Options:
#'   \itemize{
#'     \item character scalar: dam name (must be in \code{names(x$per_dam)}),
#'     \item integer scalar: 1-based index into \code{x$per_dam},
#'     \item \code{"all"}: process every dam and return a named list, or
#'     \item \code{NULL}: required to be non-NULL for multi-dam inputs; if omitted,
#'           an error is thrown to avoid ambiguity.
#'   }
#'
#' @return
#' \itemize{
#'   \item For a single-dam input: a character vector with the kept SNP names
#'         in map order. Attributes:
#'         \code{"removed"} (dropped markers),
#'         \code{"gap_threshold_cm"} (numeric threshold), and
#'         \code{"mode"} (\code{"cluster"} or \code{"per-gap"}).
#'   \item For a multi-dam input with \code{dam = "all"}: a named list of such
#'         character vectors (one per dam). The list has class
#'         \code{"HSMap.drop.multi"}.
#' }
#'
#' @examples
#' \dontrun{
#' # Single-dam example
#' kept1 <- drop_gap_markers(map, gap_cm = 12, cluster_size = 0, side = "auto")
#' attr(kept1, "removed")
#'
#' # Multi-dam example: remove small clusters up to size 3 for all dams
#' kept_all <- drop_gap_markers(multi_map, gap_cm = 10, cluster_size = 3, dam = "all")
#' kept_all$MOM2  # kept markers for dam MOM2
#' }
#' @export
drop_gap_markers <- function(x,
                             gap_cm = 10,
                             side = c("auto","right","left","both"),
                             cluster_size = 0L,
                             iterate = TRUE,
                             dam = NULL) {
  side <- match.arg(side)

  is_single <- function(obj) {
    is.list(obj) && !is.null(obj$order) && !is.null(obj$fit) && !is.null(obj$fit$r)
  }
  is_multi <- function(obj) {
    inherits(obj, "HSMap.map.multi") || (!is.null(obj$per_dam) && is.list(obj$per_dam))
  }

  if (is_single(x)) {
    return(.drop_gap_markers_one(x, gap_cm, side, as.integer(cluster_size), isTRUE(iterate)))
  }

  if (!is_multi(x)) {
    stop("`x` must be an HSMap.map (single-dam) or HSMap.map.multi (with $per_dam).")
  }

  per <- x$per_dam
  if (is.null(per) || !length(per)) stop("`x$per_dam` is empty.")

  # Require the user to disambiguate unless they ask for all.
  if (is.null(dam)) {
    stop("`x` contains multiple dams; please specify `dam` by name or index, or use dam = 'all'.")
  }

  dam_names <- names(per)

  if (identical(dam, "all")) {
    out <- lapply(seq_along(per), function(i) {
      .drop_gap_markers_one(per[[i]], gap_cm, side, as.integer(cluster_size), isTRUE(iterate))
    })
    names(out) <- if (!is.null(dam_names)) dam_names else paste0("Dam", seq_along(per))
    class(out) <- "HSMap.drop.multi"
    return(out)
  }

  if (is.character(dam) && length(dam) == 1L) {
    j <- match(dam, dam_names)
    if (is.na(j)) stop("Dam '", dam, "' not found. Available: ", paste(dam_names, collapse = ", "))
    return(.drop_gap_markers_one(per[[j]], gap_cm, side, as.integer(cluster_size), isTRUE(iterate)))
  }

  if (is.numeric(dam) && length(dam) == 1L) {
    j <- as.integer(dam)
    if (j < 1L || j > length(per)) stop("`dam` index out of range [1..", length(per), "].")
    return(.drop_gap_markers_one(per[[j]], gap_cm, side, as.integer(cluster_size), isTRUE(iterate)))
  }

  stop("`dam` must be 'all', a single name, or a single 1-based index.")
}

# --- internal worker for a single HSMap.map ----------------------------------
.drop_gap_markers_one <- function(x_one, gap_cm, side, cluster_size, iterate) {
  ord <- x_one$order
  r   <- x_one$fit$r
  if (is.null(ord) || is.null(r)) stop("Single-dam map must have fields `order` and `fit$r`.")
  if (length(r) != length(ord) - 1L)
    stop("`length(fit$r)` must equal `length(order) - 1`.")

  haldane_cm <- function(r) {
    r2 <- pmin(pmax(as.numeric(r), 1e-8), 0.499999)
    -50 * log(1 - 2 * r2)
  }
  d_cm_all <- haldane_cm(r)

  # indices (in original order) currently kept
  idx <- seq_along(ord)
  removed_idx <- integer(0)

  large_gap_flags <- function(idx_now) {
    if (length(idx_now) < 2L) return(rep(FALSE, 0L))
    left_pos  <- idx_now[-length(idx_now)]
    right_pos <- idx_now[-1L]
    is_adjacent <- (right_pos == left_pos + 1L)
    is_adjacent & (d_cm_all[left_pos] > gap_cm)
  }

  drop_clusters_once <- function(idx_now) {
    lg <- large_gap_flags(idx_now)
    if (!any(lg)) return(integer(0))
    block_id  <- 1L + cumsum(c(0L, as.integer(lg)))
    block_len <- as.integer(table(block_id))
    internal  <- seq_len(max(block_id))
    if (length(internal) <= 2L) return(integer(0))
    internal  <- internal[-c(1L, length(internal))]
    small     <- internal[block_len[internal] <= cluster_size]
    if (!length(small)) return(integer(0))
    to_drop_local <- which(block_id %in% small)
    idx_now[to_drop_local]
  }

  drop_per_gap_once <- function(idx_now) {
    lg <- large_gap_flags(idx_now)
    if (!any(lg)) return(integer(0))
    block_id  <- 1L + cumsum(c(0L, as.integer(lg)))
    block_len <- as.integer(table(block_id))
    gaps_t <- which(lg)

    if (identical(side, "both")) {
      return(idx_now[sort(unique(c(gaps_t, gaps_t + 1L)))])
    }
    if (identical(side, "left"))  return(idx_now[gaps_t])
    if (identical(side, "right")) return(idx_now[gaps_t + 1L])

    # side == "auto"
    to_drop_local <- integer(0)
    for (t in gaps_t) {
      L <- block_len[block_id[t]]
      R <- block_len[block_id[t + 1L]]
      if (L < R)        to_drop_local <- c(to_drop_local, t)
      else              to_drop_local <- c(to_drop_local, t + 1L)  # ties drop the right
    }
    idx_now[sort(unique(to_drop_local))]
  }

  mode_used <- if (cluster_size > 0L) "cluster" else "per-gap"

  repeat {
    drop_this <- integer(0)
    if (cluster_size > 0L) drop_this <- drop_clusters_once(idx)
    if (!length(drop_this)) drop_this <- drop_per_gap_once(idx)
    if (!length(drop_this)) break

    removed_idx <- c(removed_idx, drop_this)
    idx <- setdiff(idx, drop_this)

    if (!iterate) break
    if (length(idx) < 2L) break
  }

  kept <- ord[idx]
  attr(kept, "removed") <- ord[sort(unique(removed_idx))]
  attr(kept, "gap_threshold_cm") <- gap_cm
  attr(kept, "mode") <- mode_used
  kept
}
