#' Render HSMap maps with LinkageMapView (publication-quality PDF)
#'
#' @description
#' Bridge from fitted \code{HSMap.map} objects to
#' \code{LinkageMapView::lmv.linkage.plot()}, which draws the classic
#' publication figure: one vertical bar per linkage group with marker names (or,
#' with \code{denmap = TRUE}, a marker-density heat colouring) and a centimorgan
#' ruler, written to a PDF.
#'
#' Positions are computed with the same \strong{gap-aware} rule used by
#' \code{\link{plot_map_list}()}: an interval with no linkage (recombination
#' fraction at/above the no-linkage threshold), a non-finite estimate, or
#' unresolved phase is a \emph{gap}, and cumulative positions reset after it.
#' Because one continuous bar across a gap would be misleading, a map containing
#' gaps is drawn as \strong{one bar per resolved segment} (groups named
#' \code{"<LG>.1"}, \code{"<LG>.2"}, ...) when \code{split_gaps = TRUE} (the
#' default), or refused with an error when \code{split_gaps = FALSE}.
#'
#' \code{LinkageMapView} is used conditionally (listed in \code{Suggests}); the
#' function stops with an installation hint when the package is not available.
#' Note that \code{lmv.linkage.plot()} chooses its own PDF dimensions and prints
#' the sizes it used to the console.
#'
#' @param x An object of class \code{HSMap.map} or a non-empty (optionally
#'   named) list of them. List names label the linkage groups; unnamed input is
#'   labelled \code{LG1, LG2, ...}.
#' @param outfile Single file path for the PDF that
#'   \code{lmv.linkage.plot()} renders. Required (the file is created by
#'   LinkageMapView, not by this function).
#' @param map.function One of \code{"haldane"}, \code{"kosambi"},
#'   \code{"morgan"}; converts recombination fractions to centimorgans.
#' @param denmap Logical; if \code{TRUE} draw a marker-density heat map instead
#'   of marker names (LinkageMapView's \code{denmap} mode). Default \code{FALSE}.
#' @param split_gaps Logical; if \code{TRUE} (default) a map containing gap
#'   intervals is split at each gap and drawn as one bar per resolved segment.
#'   If \code{FALSE}, a map with gaps raises an error (mirroring
#'   \code{\link{plot_map_list}()}).
#' @param ... Further arguments passed to
#'   \code{LinkageMapView::lmv.linkage.plot()}, e.g. \code{mapthese},
#'   \code{autoconnadj}, \code{ruler}, \code{lgw}, \code{pdf.width},
#'   \code{pdf.height}, \code{qtldf}, \code{sectcoldf}, and the \code{cex}/font
#'   controls.
#'
#' @return Invisibly, the tidy data frame handed to LinkageMapView, with
#'   columns \code{group}, \code{position} (cM), and \code{locus} -- useful for
#'   custom annotation or for calling \code{lmv.linkage.plot()} directly.
#'
#' @seealso \code{\link{plot_map_list}} for a quick base-graphics track plot;
#'   \code{\link{get_block_map}} / \code{\link{plot_block_map}} for gap-aware
#'   blockwise summaries.
#'
#' @examplesIf requireNamespace("LinkageMapView", quietly = TRUE)
#' ## small simulated single-family map, rendered to a temporary PDF
#' sim <- sim_multi_pop(T_markers = 20, n_pops = 1, n_ind_per_pop = 80,
#'                      r_vec = rep(0.05, 19), phase_mode = "all_coupling",
#'                      maternal_geno_mode = "all_het",
#'                      paternal_pA_base = 0.4, error_rate = 0.01, seed = 1)
#' dat <- structure(list(G_list = sim$G_list, M_list = sim$M_list),
#'                  class = "HSMap.data")
#' tpt <- pairwise_rf(dat, threads = 1)
#' ph  <- phase_from_pairwise(tpt, order = sim$truth$markers_union)
#' map <- hmm_map(dat, phased = ph, epsilon = 0.01)
#' pdf_out <- tempfile(fileext = ".pdf")
#' plot_map_lmv(list(LG1 = map), outfile = pdf_out)
#' file.exists(pdf_out)
#' @export
plot_map_lmv <- function(x,
                         outfile,
                         map.function = c("haldane", "kosambi", "morgan"),
                         denmap = FALSE,
                         split_gaps = TRUE,
                         ...) {
  if (!requireNamespace("LinkageMapView", quietly = TRUE))
    stop("Package 'LinkageMapView' is required for plot_map_lmv(); install it ",
         "with install.packages(\"LinkageMapView\").", call. = FALSE)
  map.function <- match.arg(map.function)
  if (missing(outfile) || !is.character(outfile) || length(outfile) != 1L ||
      !nzchar(outfile))
    stop("`outfile` must be a single file path; LinkageMapView renders the ",
         "figure to this PDF.", call. = FALSE)

  ## normalize input exactly like plot_map_list()
  if (inherits(x, "HSMap.map")) {
    obj_list <- list(x)
  } else if (is.list(x) && length(x) > 0) {
    ok <- vapply(x, function(y) inherits(y, "HSMap.map"), logical(1))
    if (!all(ok)) stop("All elements of 'x' must be HSMap.map objects.", call. = FALSE)
    obj_list <- x
  } else {
    stop("'x' must be an HSMap.map or a non-empty list of HSMap.map.", call. = FALSE)
  }
  lg_names <- names(obj_list)
  if (is.null(lg_names) || any(!nzchar(lg_names))) {
    lg_names <- paste0("LG", seq_along(obj_list))
  }

  rows <- vector("list", length(obj_list))
  for (i in seq_along(obj_list)) {
    pos    <- get_map(obj_list[[i]], map.function = map.function)
    n_gaps <- attr(pos, "n_gaps") %||% 0L
    if (n_gaps > 0L && !isTRUE(split_gaps))
      stop("Map '", lg_names[i], "' contains ", n_gaps, " gap interval(s) ",
           "(no linkage or unresolved phase) and `split_gaps = FALSE`. ",
           "Positions reset at each gap, so one continuous bar would be ",
           "misleading. Use split_gaps = TRUE (one bar per resolved segment) ",
           "or a gap-aware summary via get_block_map().", call. = FALSE)
    grp <- if (n_gaps > 0L) paste0(lg_names[i], ".", attr(pos, "block"))
           else rep(lg_names[i], length(pos))
    rows[[i]] <- data.frame(group = grp,
                            position = as.numeric(pos),
                            locus = names(pos),
                            stringsAsFactors = FALSE)
  }
  df <- do.call(rbind, rows)
  rownames(df) <- NULL

  LinkageMapView::lmv.linkage.plot(mapthis = df, outfile = outfile,
                                   denmap = denmap, ...)
  invisible(df)
}
