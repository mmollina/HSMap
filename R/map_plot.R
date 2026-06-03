#' Plot one or more HSMap maps as linkage-group tracks
#'
#' Plot a single `HSMap.map` or a list of `HSMap.map` objects as simple tracks.
#' Each track is a rounded rectangle representing the linkage group with small
#' ticks at marker positions. The function returns a tidy data frame with the
#' plotted positions (invisibly), which can be useful for labeling or custom
#' annotations.
#'
#' @param x An object of class `HSMap.map` or a list of such objects.
#'   If a list is named, those names are used as linkage-group labels.
#' @param horiz Logical. Plot tracks horizontally if `TRUE`, vertically if `FALSE`.
#' @param map.function Character. One of `"haldane"`, `"kosambi"`, or `"morgan"`.
#'   Used to transform recombination fractions to map distance in centiMorgans.
#' @param col A vector of colors for the tracks. Recycled as needed.
#'   Defaults to `hs_pal(length(x))`.
#'
#' @return Invisibly returns a data frame with columns:
#'   `mrk` (marker id), `LG` (linkage-group label), and `pos` (position in cM).
#'
#' @examples
#' \dontrun{
#' # One map:
#' plot_map_list(my_map)
#'
#' # Multiple maps side by side:
#' plot_map_list(list(LG1 = map1, LG2 = map2), horiz = TRUE, map.function = "kosambi")
#' }
#'
#' @importFrom graphics plot axis rect lines
#' @export
plot_map_list <- function(x,
                          horiz = FALSE,
                          map.function = c("haldane", "kosambi", "morgan"),
                          col = hs_pal(length(x))) {
  map.function <- match.arg(map.function)

  # Normalize input to a list of HSMap.map
  if (inherits(x, "HSMap.map")) {
    obj_list <- list(x)
  } else if (is.list(x) && length(x) > 0) {
    ok <- vapply(x, function(y) inherits(y, "HSMap.map"), logical(1))
    if (!all(ok)) stop("All elements of 'x' must be HSMap.map objects.", call. = FALSE)
    obj_list <- x
  } else {
    stop("'x' must be an HSMap.map or a non-empty list of HSMap.map.", call. = FALSE)
  }

  # Map names
  lg_names <- names(obj_list)
  if (is.null(lg_names) || any(!nzchar(lg_names))) {
    lg_names <- paste0("LG", seq_along(obj_list))
  }

  # Positions for each map
  pos_list <- lapply(obj_list, get_map, map.function = map.function)

  # Basic sanity checks on positions
  lens <- vapply(pos_list, length, integer(1))
  if (any(lens == 0L)) stop("At least one map has zero markers after transformation.", call. = FALSE)

  # Determine range for axes
  max_dist <- max(vapply(pos_list, function(v) max(v, na.rm = TRUE), numeric(1)))
  if (!is.finite(max_dist) || max_dist < 0) stop("Invalid distances computed for maps.", call. = FALSE)

  # Colors recycled to number of maps
  col <- rep_len(col, length(pos_list))

  # Collect a tidy data frame while plotting
  out <- vector("list", length(pos_list))

  if (isTRUE(horiz)) {
    graphics::plot(0,
                   xlim = c(0, max_dist),
                   ylim = c(0, length(pos_list) + 1),
                   type = "n", axes = FALSE,
                   xlab = "Map position (cM)",
                   ylab = "Linkage groups"
    )
    graphics::axis(1)

    for (i in seq_along(pos_list)) {
      d <- pos_list[[i]]
      out[[i]] <- data.frame(mrk = names(d), LG = lg_names[i], pos = d, stringsAsFactors = FALSE)
      plot_one_map(d, i = i, horiz = TRUE, col = col[i])
    }
    graphics::axis(2, at = seq_along(pos_list), labels = lg_names, lwd = 0, las = 2)

  } else {
    graphics::plot(0,
                   ylim = c(-max_dist, 0),
                   xlim = c(0, length(pos_list) + 1),
                   type = "n", axes = FALSE,
                   ylab = "Map position (cM)",
                   xlab = "Linkage groups"
    )
    # Pretty y axis with positive labels for absolute distance
    at_y <- graphics::axis(2, labels = FALSE, lwd = 0)
    graphics::axis(2, at = at_y, labels = abs(at_y))

    for (i in seq_along(pos_list)) {
      d <- pos_list[[i]]
      out[[i]] <- data.frame(mrk = names(d), LG = lg_names[i], pos = d, stringsAsFactors = FALSE)
      plot_one_map(d, i = i, horiz = FALSE, col = col[i])
    }
    graphics::axis(3, at = seq_along(pos_list), labels = lg_names, lwd = 0, las = 2)
  }

  invisible(do.call(rbind, out))
}


#' Draw a single linkage-group track
#'
#' Low-level helper that draws one track given a vector of marker positions.
#' Used by [plot_map_list()]. Not exported.
#'
#' @param x Named numeric vector of cumulative positions in cM for a single map.
#'   Names are used for tick placement but are not drawn as text.
#' @param i Integer index of the track position in the panel.
#' @param horiz Logical. Horizontal if `TRUE`, vertical if `FALSE`.
#' @param col Fill color for the track rectangle.
#'
#' @return Invisibly returns `NULL`.
#' @keywords internal
plot_one_map <- function(x,
                         i = 0,
                         horiz = FALSE,
                         col = "lightgray") {
  if (!is.numeric(x) || !length(x)) return(invisible(NULL))

  half_w <- 0.25

  if (isTRUE(horiz)) {
    # Track body
    graphics::rect(xleft = x[1], ybottom = i - half_w,
                   xright = utils::tail(x, 1), ytop = i + half_w,
                   col = col, border = NA)
    # Marker ticks
    for (j in seq_along(x)) {
      graphics::lines(x = c(x[j], x[j]), y = c(i - half_w, i + half_w), lwd = 0.5)
    }
  } else {
    y <- -rev(x)  # flip so origin at top
    graphics::rect(xleft = i - half_w, ybottom = y[1],
                   xright = i + half_w, ytop = utils::tail(y, 1),
                   col = col, border = NA)
    for (j in seq_along(y)) {
      graphics::lines(y = c(y[j], y[j]), x = c(i - half_w, i + half_w), lwd = 0.5)
    }
  }

  invisible(NULL)
}


#' Compute map positions for an HSMap.map
#'
#' Convert recombination fractions to cumulative positions in cM according to
#' the chosen mapping function. Used internally by [plot_map_list()]. Not exported.
#'
#' @param map An object of class `HSMap.map` containing `map$fit$r` (vector of
#'   recombination fractions between consecutive markers) and `map$order`
#'   (character vector with marker names in map order).
#' @param map.function Character. One of `"haldane"`, `"kosambi"`, or `"morgan"`.
#'
#' @return A named numeric vector of cumulative positions in cM
#'   with names equal to `map$order`.
#' @keywords internal
get_map <- function(map, map.function = c("haldane", "kosambi", "morgan")) {
  map.function <- match.arg(map.function)

  if (!inherits(map, "HSMap.map")) {
    stop("'map' must be an HSMap.map object.", call. = FALSE)
  }
  if (is.null(map$fit) || is.null(map$fit$r)) {
    stop("'map$fit$r' is missing.", call. = FALSE)
  }
  if (is.null(map$order)) {
    stop("'map$order' is missing.", call. = FALSE)
  }

  mf <- switch(
    map.function,
    haldane = inv_haldane,
    kosambi = inv_kosambi,
    morgan  = inv_morgan
  )
  if (!is.function(mf)) {
    stop("Mapping function not found in the namespace: ", map.function, call. = FALSE)
  }

  r <- as.numeric(map$fit$r)
  if (any(!is.finite(r))) stop("Non-finite values in 'map$fit$r'.", call. = FALSE)
  if (any(r < 0 | r > 0.5, na.rm = TRUE)) {
    stop("Recombination fractions must be in [0, 0.5].", call. = FALSE)
  }

  pos <- cumsum(mf(c(0, r)))
  ord <- as.character(map$order)

  if (length(pos) != length(ord)) {
    stop("Length mismatch between positions and 'map$order'.", call. = FALSE)
  }

  stats::setNames(pos, ord)
}
