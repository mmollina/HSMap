#' Plot recombination or LOD matrix from an HSMap two-point result
#'
#' @description
#' Heatmap of either the recombination fraction matrix (\code{r}) or the LOD
#' matrix (\code{lod_r}) stored in an \code{HSMap.tpt} object. You can reorder
#' and/or remove markers, and optionally aggregate (downsample) the matrix to
#' speed up plotting for large datasets.
#'
#' @param tpt An \code{HSMap.tpt} object (as returned by \code{hsmap_twopoint()},
#'   containing two-point results under \code{tpt$fit}), or a plain list with
#'   square matrices \code{r} and \code{lod_r} that share identical dimnames.
#' @param type Either \code{"rf"} for recombination fraction or \code{"lod"} for LOD.
#' @param ord Optional ordering for markers. Either a character vector of marker
#'   IDs or a numeric index vector. The matrix is permuted as \code{w[ord, ord]}.
#'   Names not present are ignored with a warning.
#' @param rem Optional character vector of marker IDs to remove (dropped from rows
#'   and columns) before plotting.
#' @param main.text Optional plot title. If \code{NULL}, a default is chosen
#'   based on \code{type}.
#' @param index Logical. If \code{TRUE}, draw marker labels on heatmap cells.
#'   (For large matrices, the label size is auto-reduced.)
#' @param fact Integer aggregation factor. If \code{> 1}, the matrix is
#'   downsampled by averaging non-overlapping \code{fact x fact} blocks with
#'   \code{\link{aggregate_matrix}} (useful for very large matrices).
#' @param useRaster Logical passed to \code{fields::image.plot}. Setting
#'   \code{TRUE} is often faster for large matrices. Default \code{TRUE}.
#' @param ... Passed through to \code{fields::image.plot()}.
#'
#' @details
#' For \code{type = "lod"}, the matrix is \emph{log10}-transformed for plotting,
#' and the legend tick labels are shown back in the original LOD scale for
#' interpretability.
#'
#' @return Invisibly returns a list with:
#' \itemize{
#'   \item \code{w}: The matrix that was actually plotted (after any reordering,
#'         filtering, transformation, and aggregation).
#'   \item \code{markers}: The marker IDs in the displayed order.
#' }
#'
#' @seealso \code{\link{aggregate_matrix}}
#'
#' @examples
#' \dontrun{
#' # Basic RF heatmap
#' pairwise_heatmap(tpt, type = "rf")
#'
#' # LOD heatmap with a custom order and aggregation for speed
#' ord <- hsmap_mds_order(grp, tpt, lg = 1, plot_each = FALSE)[["LG1"]]
#' pairwise_heatmap(tpt, type = "lod", ord = ord, fact = 4)
#' }
#' @importFrom fields image.plot tim.colors
#' @export
pairwise_heatmap <- function(tpt,
                                  type = c("rf", "lod"),
                                  ord = NULL,
                                  rem = NULL,
                                  main.text = NULL,
                                  index = FALSE,
                                  fact = 1L,
                                  useRaster = TRUE,
                                  ...) {
  type <- match.arg(type)

  # extract matrices from HSMap.tpt or raw list
  src <- if (inherits(tpt, "HSMap.tpt")) tpt$fit else tpt
  if (is.null(src$r) || is.null(src$lod_r)) {
    stop("Input must contain matrices `r` and `lod_r` (for HSMap.tpt, in `tpt$fit`).")
  }
  R <- src$r
  L <- src$lod_r
  if (!is.matrix(R) || !is.matrix(L)) stop("`r` and `lod_r` must be matrices.")
  if (!identical(rownames(R), colnames(R)) ||
      !identical(dim(R), dim(L)) ||
      !identical(rownames(R), rownames(L)) ||
      !identical(colnames(R), colnames(L))) {
    stop("`r` and `lod_r` must be square with identical dimnames.")
  }

  # choose matrix to display
  if (type == "rf") {
    w <- R
    if (is.null(main.text)) main.text <- "Recombination fraction matrix"
  } else {
    w <- L
    if (is.null(main.text)) main.text <- "log10(LOD) score matrix"
  }
  if (inherits(ord, what = "hsmap_mds")){
    ord <- unlist(ord)
  }
  # optional ordering
  if (!is.null(ord)) {
    if (is.character(ord)) {
      miss <- setdiff(ord, colnames(w))
      if (length(miss)) {
        warning(length(miss), " markers in `ord` not found; ignoring them. Example: ",
                paste(utils::head(miss, 6), collapse = ", "),
                if (length(miss) > 6) " ..." else "")
      }
      ord <- intersect(ord, colnames(w))
      if (!length(ord)) stop("No markers from `ord` were found in the matrix.")
      w <- w[ord, ord, drop = FALSE]
    }
    else if (is.numeric(ord)) {
      ord <- as.integer(ord)
      if (any(ord < 1L | ord > ncol(w))) stop("Numeric `ord` has indices out of range.")
      w <- w[ord, ord, drop = FALSE]
    }
    else {
      stop("`ord` must be a character vector of marker IDs or a numeric index vector.")
    }
  }

  # optional removal
  if (!is.null(rem)) {
    drop_idx <- which(colnames(w) %in% rem)
    if (length(drop_idx)) w <- w[-drop_idx, -drop_idx, drop = FALSE]
  }

  # optional aggregation
  if (!is.numeric(fact) || length(fact) != 1L || is.na(fact) || fact < 1L) {
    stop("`fact` must be a positive integer.")
  }
  if (fact > 1L) {
    w <- aggregate_matrix(w, fact = fact)$R
  }

  # color scale and legend
  if (type == "rf") {
    # palette length tied to displayed RF range (but capped)
    max_rf <- max(w, na.rm = TRUE)
    ncols  <- max(2L, min(128L, ceiling(128 * max_rf) + 1L))
    col.range <- na.omit(rev(fields::tim.colors()))[seq_len(ncols)]
    lab_breaks <- NULL
  } else {
    # clip, log10-transform for plotting (labels back in LOD units)
    w[w < 1e-4] <- 1e-4
    w <- log10(w)
    max_lod <- max(10^w, na.rm = TRUE)
    ncols   <- max(10L, min(128L, ceiling(128 * max_lod) + 1L))
    col.range <- na.omit(fields::tim.colors()[seq_len(ncols)])
    if (all(is.finite(w))) {
      brks_log <- seq(min(w, na.rm = TRUE), max(w, na.rm = TRUE), length.out = 11)
      lab_breaks <- round(10^brks_log, 1)
    } else {
      lab_breaks <- NULL
    }
  }

  # draw heatmap
  fields::image.plot(
    w,
    col        = col.range,
    lab.breaks = lab_breaks,
    main       = main.text,
    useRaster  = useRaster,
    axes       = FALSE,
    ...
  )

  # optional labels
  if (isTRUE(index)) {
    ft <- if (ncol(w) < 100) 0.7 else 100 / ncol(w)
    graphics::text(
      x = seq(0, 1, length.out = ncol(w)),
      y = seq(0, 1, length.out = ncol(w)),
      labels = colnames(w),
      cex = ft
    )
  }

  invisible(list(w = w, markers = colnames(w)))
}
