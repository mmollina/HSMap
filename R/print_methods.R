# Concise print methods for the primary HSMap result classes, so interactive users get
# a readable one-screen summary instead of a raw list dump.

`%||%` <- function(a, b) if (is.null(a)) b else a

#' Print methods for HSMap result objects
#'
#' Concise one-screen summaries for the primary HSMap result classes
#' (two-point results, phase objects, single-dam and blockwise maps, and
#' per-offspring haplotype posteriors).
#'
#' @param x An HSMap result object.
#' @param ... Ignored.
#' @return \code{x}, invisibly.
#' @name hsmap-print
NULL

#' @rdname hsmap-print
#' @exportS3Method print HSMap.tpt
print.HSMap.tpt <- function(x, ...) {
  n_mk <- length(x$markers %||% x$fit$markers)
  n_dams <- x$fit$n_dams %||% length(x$fit$mom_phase_list %||% list())
  cat("HSMap.tpt (two-point result)\n")
  cat("  Markers      : ", n_mk, "\n", sep = "")
  cat("  Dams         : ", n_dams, "\n", sep = "")
  if (!is.null(x$time_sec)) cat("  Compute time : ", round(x$time_sec, 2), " s\n", sep = "")
  if (!is.null(x$filter_stats))
    cat("  Filtered     : yes (", length(x$filter_stats$kept_markers %||% character(0)),
        " markers kept)\n", sep = "")
  cat("  Next         : group_markers() -> mds_order() -> phase_from_pairwise()\n")
  invisible(x)
}

#' @rdname hsmap-print
#' @exportS3Method print HSMap.phased
print.HSMap.phased <- function(x, ...) {
  Ti <- length(x$phase_vec)
  cat("HSMap.phased (dam '", x$dam %||% "?", "')\n", sep = "")
  cat("  Markers            : ", length(x$order), "\n", sep = "")
  cat("  Resolved intervals : ", sum(!is.na(x$phase_vec)), " / ", Ti, "\n", sep = "")
  cat("  Phase components   : ", x$n_components %||% NA, "\n", sep = "")
  invisible(x)
}

#' @rdname hsmap-print
#' @exportS3Method print HSMap.phased.multi
print.HSMap.phased.multi <- function(x, ...) {
  cat("HSMap.phased.multi (", length(x), " dam(s))\n", sep = "")
  for (nm in names(x)) {
    p <- x[[nm]]; Ti <- length(p$phase_vec)
    cat("  ", nm, ": ", sum(!is.na(p$phase_vec)), "/", Ti, " intervals resolved, ",
        p$n_components %||% NA, " component(s)\n", sep = "")
  }
  invisible(x)
}

#' @rdname hsmap-print
#' @exportS3Method print HSMap.map
print.HSMap.map <- function(x, ...) {
  fit <- x$fit
  cat("HSMap.map (dam '", x$dam %||% "?", "')\n", sep = "")
  cat("  Markers   : ", length(x$order), "\n", sep = "")
  if (!is.null(fit$r)) cat("  Intervals : ", length(fit$r),
                           " (", sum(as.numeric(fit$r) >= 0.5 - 1e-9), " at no linkage)\n", sep = "")
  if (!is.null(fit$converged)) cat("  Converged : ", isTRUE(fit$converged),
                                   " (", fit$iters %||% NA, " iters)\n", sep = "")
  cat("  Summarise : get_map() / plot_map_list() (gap-safe cM positions)\n")
  invisible(x)
}

#' @rdname hsmap-print
#' @exportS3Method print HSMap.map.blocks
print.HSMap.map.blocks <- function(x, ...) {
  cat("HSMap.map.blocks (blockwise multipoint map)\n")
  cat("  Markers  : ", length(x$order), "\n", sep = "")
  cat("  Blocks   : ", x$n_blocks %||% length(x$blocks), " resolved phase block(s)\n", sep = "")
  cat("  Dams     : ", length(x$dams %||% character(0)), "\n", sep = "")
  cat("  gap_r    : ", x$gap_r %||% NA, "\n", sep = "")
  cat("  Summarise: get_block_map() / plot_block_map() (gaps -> NA, never huge cM)\n")
  invisible(x)
}

#' @rdname hsmap-print
#' @exportS3Method print HSMap.gamma
print.HSMap.gamma <- function(x, ...) {
  d <- dim(x$gamma)
  cat("HSMap.gamma (per-offspring haplotype posteriors, dam '", x$dam %||% "?", "')\n", sep = "")
  cat("  Array [homolog x marker x offspring]: ", paste(d, collapse = " x "), "\n", sep = "")
  invisible(x)
}
