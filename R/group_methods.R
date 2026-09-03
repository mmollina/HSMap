
#' Plot an hsmap_group as a colored dendrogram with cluster rectangles
#'
#' @description
#' Colors branches by group (using `palette(k)`), draws `rect.hclust()` boxes,
#' and overlays group sizes mappoly2-style. Useful to visually assess the cut.
#'
#' @param x An object of class `"hsmap_group"`.
#' @param border Rectangle border color. Default `"black"`.
#' @param lwd_base Base line width for the dendrogram. Default `1.5`.
#' @param palette A function `f(n)` returning `n` colors (default [hs_pal()]).
#' @param ... Passed to `plot.dendrogram()`.
#' @return Invisibly returns a list of length `k`, each the marker names inside
#'   each rectangle from left to right.
#' @export
#' @importFrom stats as.dendrogram rect.hclust
plot.hsmap_group <- function(x,
                             border = "black",
                             lwd_base = 1.5,
                             palette = hs_pal,
                             ...) {
  stopifnot(inherits(x, "hsmap_group"))
  hc <- x$hc; k <- x$k
  if (is.null(hc) || is.null(k)) stop("Invalid hsmap_group: missing `hc` or `k`.")

  dend <- stats::as.dendrogram(hc)

  op <- graphics::par(lwd = lwd_base)
  on.exit(graphics::par(op), add = TRUE)

  if (requireNamespace("dendextend", quietly = TRUE)) {
    dend <- dendextend::color_branches(dend, k = k, col = palette(k))
  }
  plot(dend, leaflab = "none", ...)

  graphics::par(lwd = 4)
  rect.hclust(hc, k = k, border = border)
  graphics::par(lwd = lwd_base)

  # left-to-right membership "rectangles" (by run-length of cutree order)
  ord_names <- labels(dend)
  g_ord     <- as.integer(x$groups[ord_names])
  runs      <- rle(g_ord)
  ends      <- cumsum(runs$lengths)
  starts    <- ends - runs$lengths + 1L
  rects <- lapply(seq_along(starts), function(i) ord_names[starts[i]:ends[i]])

  # draw big dots with sizes
  lens <- vapply(rects, length, integer(1))
  xt   <- cumsum(lens) - ceiling(lens / 2)
  yt   <- 0.1
  graphics::points(x = xt, y = rep(yt, length(xt)), cex = 6, pch = 20, col = "lightgray")
  graphics::text(x = xt, y = yt, labels = lens, adj = 0.5, cex = .7)

  invisible(rects)
}

#' Print summary for hsmap_group
#'
#' @param x An `"hsmap_group"` object.
#' @param show_markers Integer; print up to this many marker names per group
#'   (in left-to-right order). Set 0 to suppress. Default `5`.
#' @param ... Unused.
#' @return Invisibly returns `x`.
#' @export
print.hsmap_group <- function(x, show_markers = 5, ...) {
  stopifnot(inherits(x, "hsmap_group"))
  k <- x$k; g <- x$groups
  if (is.null(k) || is.null(g)) stop("Invalid hsmap_group: missing `k` or `groups`.")
  n_markers <- length(x$markers %||% names(g))
  sizes <- table(factor(g, levels = seq_len(k)))

  cat("hsmap_group\n")
  cat("  Method  : ", x$method %||% "unknown", "\n", sep = "")
  cat("  Markers : ", n_markers, "\n", sep = "")
  cat("  Groups k: ", k, "\n", sep = "")
  cat("  Group sizes:\n")
  print(data.frame(Group = paste0("LG", seq_len(k)),
                   N = as.integer(sizes), row.names = NULL))
  if (is.numeric(show_markers) && show_markers > 0) {
    cat("  Example markers per LG (up to ", show_markers, "):\n", sep = "")
    sp <- split(names(g), g)
    for (i in seq_len(k)) {
      ids <- sp[[as.character(i)]]
      if (length(ids)) {
        show <- utils::head(ids, show_markers)
        more <- if (length(ids) > length(show)) " ..." else ""
        cat("   - LG", i, ": ", paste(show, collapse = ", "), more, "\n", sep = "")
      } else {
        cat("   - LG", i, ": (empty)\n", sep = "")
      }
    }
  }
  if (!is.null(x$comp)) {
    cat("  Comparison matrix available in `$comp` (LG x chromosome)\n")
  }
  cat("  Tip: plot(x) to inspect the dendrogram by groups\n")
  invisible(x)
}
