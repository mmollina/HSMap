#' Print method for hsmap_mds
#'
#' @param x An \code{hsmap_mds} object.
#' @param show_bounds Integer; show this many markers from the start and end of each order.
#'   Default \code{3}.
#' @param max_groups Integer; print details for at most this many groups (the rest summarized).
#'   Default \code{12}.
#' @param ... Unused.
#' @return Invisibly returns \code{x}.
#' @export
print.hsmap_mds <- function(x, show_bounds = 3, max_groups = 12, ...) {
  stopifnot(inherits(x, "hsmap_mds") || (is.list(x) && !is.null(names(x))))

  `%||%` <- function(a, b) if (is.null(a)) b else a
  lg_names <- names(x) %||% paste0("LG", seq_along(x))
  n_groups <- length(x)
  sizes <- vapply(x, length, integer(1))
  total_markers <- sum(sizes)

  cat("hsmap_mds\n")
  cat("  Linkage groups: ", n_groups, "\n", sep = "")
  cat("  Total markers across groups: ", total_markers, "\n", sep = "")

  # summary table
  summ <- data.frame(LG = lg_names, N = sizes, stringsAsFactors = FALSE)
  ord <- seq_len(n_groups)
  suppressWarnings({
    lg_num <- as.integer(gsub("^[^0-9]*", "", lg_names))
    if (all(!is.na(lg_num))) ord <- order(lg_num)
  })
  summ <- summ[ord, , drop = FALSE]; rownames(summ) <- NULL
  cat("  Group sizes:\n")
  print(summ, row.names = FALSE)

  # details per group up to max_groups
  to_show <- min(n_groups, max_groups)
  if (to_show > 0) {
    cat("  Endpoints per group (first ", show_bounds, " and last ", show_bounds, " markers):\n", sep = "")
    for (i in seq_len(to_show)) {
      nm <- summ$LG[i]
      ord_vec <- x[[nm]] %||% x[[match(nm, lg_names)]]
      if (length(ord_vec) == 0) {
        cat("   - ", nm, ": (empty)\n", sep = ""); next
      }
      head_str <- paste(utils::head(ord_vec, show_bounds), collapse = ", ")
      tail_str <- paste(utils::tail(ord_vec, show_bounds), collapse = ", ")
      cat("   - ", nm, " [", length(ord_vec), "]: ",
          "[", head_str, "] ... [", tail_str, "]\n", sep = "")
    }
    if (n_groups > to_show) {
      cat("  ... ", n_groups - to_show, " more group(s) omitted\n", sep = "")
    }
  }

  invisible(x)
}
