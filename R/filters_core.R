#' Filter a two-point result (HSMap.tpt) by RF/LOD, with optional pre-ordering
#'
#' @description
#' This function filters the **two-point** output from HSMap by:
#' 1) optionally reordering markers by a user-supplied sequence `order`,
#' 2) keeping only pairs with LOD >= threshold and RF <= threshold
#'    (optionally within a diagonal band around the matrix diagonal), and
#' 3) retaining only those markers that have a sufficient number of surviving,
#'    non-`NA` pairs (controlled by a row-count quantile window).
#'
#' It returns the same \code{"HSMap.tpt"} object, **restricted to the retained
#' markers** and preserving the chosen order. All companion matrices
#' (`r`, `lod_r`, `lod_ph`, `logLik`) and each element of `mom_phase_list`
#' are subset consistently. The top-level fields `markers` and `fit$markers`
#' are updated. A `filter_stats` list is attached with thresholds and the kept IDs.
#'
#' @param x An object of class \code{"HSMap.tpt"} whose matrices live under
#'   \code{x$fit}: \code{r}, \code{lod_r}, \code{lod_ph}, \code{logLik},
#'   and \code{mom_phase_list}.
#' @param order Optional character vector of marker IDs to **reindex** the
#'   matrices before filtering (useful if you already have a map order).
#'   Markers not present are ignored (you'll get a warning); the filter then
#'   runs on the reindexed square submatrix.
#' @param thresh.LOD.ph Numeric. Kept for parity; we actually screen by
#'   \code{lod_r}. The effective LOD cutoff is
#'   \code{max(thresh.LOD.ph, thresh.LOD.rf)}. Default \code{5}.
#' @param thresh.LOD.rf Numeric LOD threshold for linkage (vs \eqn{r=0.5}).
#'   Combined with \code{thresh.LOD.ph} by \code{max()}. Default \code{5}.
#' @param thresh.rf Numeric RF cutoff; only pairs with \code{r <= thresh.rf}
#'   are considered informative. Default \code{0.15}.
#' @param probs Numeric length-2 vector of quantiles for the **row non-NA counts**
#'   after masking (see Details). Markers outside \code{[probs[1], probs[2]]}
#'   are dropped. Default \code{c(0.05, 1)} keeps the top 95% by information.
#' @param diag.markers Optional non-negative integer half-bandwidth. If supplied,
#'   we keep only entries with \code{|i - j| <= diag.markers} **after any reordering**.
#'   This is useful when you only want to keep near-diagonal (putative linked) pairs.
#' @param diagnostic.plot Logical. If \code{TRUE}, plot histograms of the per-marker
#'   non-NA counts before/after keeping markers. Default \code{TRUE}.
#' @param breaks Integer number of bins for the diagnostic histogram. Default \code{100}.
#'
#' @details
#' **Masking logic.** We build a logical mask \eqn{M_{ij}} that is \code{TRUE} when:
#' \deqn{ \text{LOD}_{ij} \ge \text{LOD\_thr} \;\;\text{and}\;\; r_{ij} \le \text{RF\_thr} }
#' where \eqn{\text{LOD\_thr} = \max(\text{thresh.LOD.ph}, \text{thresh.LOD.rf})}
#' and \eqn{\text{RF\_thr} = \text{thresh.rf}}. If \code{diag.markers} is given,
#' we further intersect with the diagonal band \eqn{|i-j|\le \text{diag.markers}}.
#'
#' For each marker \(i\), we count how many \code{TRUE} entries remain in row \(i\).
#' We then keep markers whose count lies within the quantile window
#' \code{[probs[1], probs[2]]} of this distribution. This drops outlier rows that are
#' **too sparse** (uninformative) or, if you set the upper quantile < 1, **too dense** (often
#' indicative of artefacts).
#'
#' @return The same object \code{x} (class \code{"HSMap.tpt"}) with:
#' \itemize{
#'   \item \code{x$fit$r}, \code{x$fit$lod_r}, \code{x$fit$lod_ph}, \code{x$fit$logLik}
#'         subset to the kept markers (square, same order).
#'   \item \code{x$fit$mom_phase_list[[g]]} subset similarly for each population.
#'   \item \code{x$fit$markers} and \code{x$markers} updated to the kept order.
#'   \item \code{x$filter_stats}: list with thresholds, counts, and \code{kept_markers}.
#' }
#'
#' @examples
#' \dontrun{
#' # Keep near-diagonal informative pairs and the best 90% markers by row counts
#' tpt2 <- tpt_filter(
#'   x = tpt,
#'   order = known_order,          # optional pre-order
#'   thresh.LOD.rf = 6,
#'   thresh.rf     = 0.2,
#'   diag.markers  = 200,
#'   probs         = c(0.10, 1),
#'   diagnostic.plot = TRUE
#' )
#' ncol(tpt2$fit$r)   # number of kept markers
#' }
#' @export
tpt_filter <- function(x,
                             order = NULL,
                             thresh.LOD.ph = 5,
                             thresh.LOD.rf = 5,
                             thresh.rf     = 0.15,
                             probs         = c(0.05, 1),
                             diag.markers  = NULL,
                             diagnostic.plot = TRUE,
                             breaks = 100) {
  if (!inherits(x, "HSMap.tpt"))
    stop("`x` must be an object of class 'HSMap.tpt'.")

  # --- pull matrices from the HSMap.tpt structure ----------------------------
  fit <- x$fit
  req <- c("r", "lod_r")
  if (!all(req %in% names(fit)))
    stop("`x$fit` must contain matrices `r` and `lod_r`.")
  R_all <- fit$r
  L_all <- fit$lod_r
  if (!is.matrix(R_all) || !is.matrix(L_all))
    stop("`x$fit$r` and `x$fit$lod_r` must be matrices.")
  if (is.null(rownames(R_all)) || !identical(rownames(R_all), colnames(R_all)))
    stop("`x$fit$r` must be square with matching dimnames.")
  if (!identical(dim(R_all), dim(L_all)) ||
      !identical(rownames(R_all), rownames(L_all)) ||
      !identical(colnames(R_all), colnames(L_all)))
    stop("`r` and `lod_r` must have identical dimensions and dimnames.")

  # --- optional pre-order (reindex) ------------------------------------------
  if (!is.null(order)) {
    order <- as.character(order)
    o <- intersect(order, rownames(R_all))
    if (length(o) < 2L) stop("Fewer than 2 markers from `order` are present in matrices.")
    if (length(o) < length(order)) {
      miss <- setdiff(order, o)
      warning(length(miss), " marker(s) from `order` not found; e.g., ",
              paste(utils::head(miss, 5), collapse = ", "))
    }
    R <- R_all[o, o, drop = FALSE]
    L <- L_all[o, o, drop = FALSE]
  } else {
    R <- R_all
    L <- L_all
  }

  # --- thresholds & mask ------------------------------------------------------
  lod_thr <- max(thresh.LOD.ph, thresh.LOD.rf)
  if (!is.numeric(thresh.rf) || thresh.rf < 0 || thresh.rf > 0.5)
    stop("`thresh.rf` must be in [0, 0.5].")

  ok <- !is.na(R) & !is.na(L) & (L >= lod_thr) & (R <= thresh.rf)

  # diagonal band (optional)
  if (!is.null(diag.markers)) {
    bw <- as.integer(diag.markers)
    if (!is.finite(bw) || bw < 0) stop("`diag.markers` must be a non-negative integer.")
    band <- abs(row(ok) - col(ok)) <= bw
    ok <- ok & band
  }

  # --- row non-NA counts & quantile window -----------------------------------
  z <- rowSums(ok, na.rm = TRUE)
  probs <- sort(probs)
  if (length(probs) != 2L || any(!is.finite(probs)) || probs[1] < 0 || probs[2] > 1)
    stop("`probs` must be a length-2 numeric vector within [0,1].")
  th <- stats::quantile(z, probs = probs, na.rm = TRUE)
  keep_ids <- names(z)[z >= th[1] & z <= th[2]]
  if (!length(keep_ids)) stop("No markers remain after filtering.")

  # keep order currently in R (which already reflects `order` if supplied)
  ids <- rownames(R)[rownames(R) %in% keep_ids]

  # --- diagnostics ------------------------------------------------------------
  if (isTRUE(diagnostic.plot)) {
    # compute a reasonable binwidth from the full histogram
    h <- graphics::hist(z, breaks = breaks, plot = FALSE)
    bw <- if (length(h$mids) > 1) diff(h$mids)[1] else 1
    dat <- rbind(
      data.frame(stage = "all",      value = z),
      data.frame(stage = "filtered", value = z[ids])
    )
    if (requireNamespace("ggplot2", quietly = TRUE)) {
      p <- ggplot2::ggplot(dat, ggplot2::aes(value, fill = stage)) +
        ggplot2::geom_histogram(alpha = 0.45, position = "identity", binwidth = bw) +
        ggplot2::scale_fill_manual(values = c(all = "#00AFBB", filtered = "#E7B800")) +
        ggplot2::labs(
          title = sprintf("Row non-NA counts after mask (LOD >= %.2f, r <= %.3f)", lod_thr, thresh.rf),
          x = "Count of non-NA in row", y = "Frequency"
        ) +
        ggplot2::theme_minimal(12)
      print(p)
    } else {
      graphics::hist(z, breaks = breaks, main = "Row non-NA counts (masked)",
                     xlab = "count")
      graphics::abline(v = th, col = "red", lty = 2)
    }
  }

  # --- helper to subset any square matrix to `ids` ----------------------------
  subset_sq <- function(mat) {
    if (!is.matrix(mat)) return(mat)
    rn <- rownames(mat); cn <- colnames(mat)
    if (is.null(rn) || is.null(cn)) return(mat)
    keep <- intersect(ids, rn)
    mat[keep, keep, drop = FALSE]
  }

  # --- write back into x (HSMap.tpt) ------------------------------------------
  x$fit$r      <- subset_sq(x$fit$r)
  x$fit$lod_r  <- subset_sq(x$fit$lod_r)
  if (!is.null(x$fit$lod_ph))  x$fit$lod_ph  <- subset_sq(x$fit$lod_ph)
  if (!is.null(x$fit$logLik))  x$fit$logLik  <- subset_sq(x$fit$logLik)

  if (!is.null(x$fit$mom_phase_list) && length(x$fit$mom_phase_list)) {
    for (nm in names(x$fit$mom_phase_list)) {
      ph <- x$fit$mom_phase_list[[nm]]
      if (is.matrix(ph) && !is.null(rownames(ph))) {
        keep <- intersect(ids, rownames(ph))
        x$fit$mom_phase_list[[nm]] <- ph[keep, keep, drop = FALSE]
      }
    }
  }

  # per-dam phase-LOD matrices (subset like mom_phase_list so all stay aligned)
  if (!is.null(x$fit$lod_ph_list) && length(x$fit$lod_ph_list)) {
    for (g in seq_along(x$fit$lod_ph_list)) {
      L <- x$fit$lod_ph_list[[g]]
      if (is.matrix(L) && !is.null(rownames(L))) {
        keep <- intersect(ids, rownames(L))
        x$fit$lod_ph_list[[g]] <- L[keep, keep, drop = FALSE]
      }
    }
  }

  # per-dam, per-marker q vectors (subset by marker name, same order as `ids`)
  if (!is.null(x$fit$q_list) && length(x$fit$q_list)) {
    for (g in seq_along(x$fit$q_list)) {
      qv <- x$fit$q_list[[g]]
      if (!is.null(names(qv))) {
        keep <- intersect(ids, names(qv))
        x$fit$q_list[[g]] <- qv[keep]
      }
    }
  }

  # update marker vectors
  x$fit$markers <- ids
  x$markers     <- ids

  # attach stats about the filtering
  x$filter_stats <- list(
    n_in_matrix  = nrow(R),
    n_kept       = length(ids),
    kept_markers = ids,
    thresholds   = list(
      LOD_used     = lod_thr,
      rf           = thresh.rf,
      probs        = probs,
      diag.markers = diag.markers
    ),
    order_supplied = !is.null(order)
  )

  class(x) <- "HSMap.tpt"
  x
}
