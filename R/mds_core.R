#' 1D MDS ordering of markers (HSMap), per group or for a given set of markers
#'
#' @description
#' Runs a 2D metric MDS (SMACOF) on a transformed RF distance, then fits a
#' principal curve and orders markers by the curve arc-length parameter. You can:
#' - Give a grouping object (\code{hsmap_group}) to order each linkage group, or
#' - Provide a character vector of marker IDs to order directly (single group).
#'
#' Optionally, markers with large perpendicular residuals (far from the curve)
#' can be dropped by an **envelope** rule whose width varies smoothly along the
#' curve (residual-based loess). Dropped marker names are stored in the
#' \code{"removed"} attribute of each returned order.
#'
#' @param x Either:
#'   \itemize{
#'     \item an \code{hsmap_group} object (with a \code{groups} or
#'           \code{groups.snp} integer vector named by marker IDs), or
#'     \item a character vector of marker IDs to order directly (then \code{lg} is ignored).
#'   }
#' @param tpt An \code{HSMap.tpt} object (preferred; matrices in \code{tpt$fit}) or
#'   a list with square matrices \code{r} and \code{lod_r} that share identical
#'   dimnames (marker IDs).
#' @param lg When \code{x} is an \code{hsmap_group}, which groups to process:
#'   \code{"all"} (default) or an integer vector of group IDs. Ignored if \code{x}
#'   is a character vector of markers.
#' @param plot_each Logical. If \code{TRUE}, plot each linkage group’s 2D MDS
#'   configuration, principal curve, and envelope while running. Default \code{TRUE}.
#' @param itmax Integer. Maximum SMACOF iterations (passed to \code{smacof::smacofSym}).
#'   Default \code{100000}.
#' @param drop_outside_envelope Logical. If \code{TRUE}, drop markers whose
#'   perpendicular residual to the principal curve exceeds
#'   \code{envelope_k * width(lambda)} (see Details). Default \code{FALSE}.
#' @param envelope_k Positive numeric multiplier for the local residual width;
#'   \code{2} is roughly a “~95%” rule under rough normality. Default \code{2}.
#' @param envelope_span Loess span used to smooth the residual width as a function
#'   of the curve parameter \code{lambda}. Default \code{0.25}.
#' @param min_keep Minimal number of markers to keep per group. If envelope
#'   filtering would leave fewer than this, the filter is **not** applied for that group.
#'   Default \code{3}.
#'
#' @details
#' \strong{Distance transform.} We transform RF into a metric dissimilarity via
#' \deqn{ d_{ij} \;=\; \frac{-50\,\log(1 - 2\,\min(\max(r_{ij}, 10^{-8}), 0.499999))}{100}, }
#' i.e. a monotone map similar to Haldane’s function, rescaled to roughly cM-like
#' units. NA entries are set to tiny positive numbers (for both distance and weights).
#'
#' \strong{Weights.} The weight matrix is \eqn{\mathrm{LOD}^2} with zero diagonal;
#' NA’s are replaced by a tiny positive number. This emphasizes pairs with strong
#' linkage signal during SMACOF.
#'
#' \strong{Envelope filter.} After fitting \code{princurve::principal_curve} on
#' the 2D configuration, we compute residual distances \eqn{r_i} from each point
#' to the curve, smooth \eqn{r_i} against \eqn{\lambda_i} (curve parameter) with
#' loess, and drop points with \eqn{r_i > k \cdot \widehat{w}(\lambda_i)}. If the
#' result would keep fewer than \code{min_keep} markers, we skip filtering for that group.
#'
#' @return A named list of class \code{"hsmap_mds"} where each element is a character
#' vector with the ordered marker names for one group (or for the one provided set).
#' Each element has attribute \code{"removed"} = names of markers dropped by the envelope
#' (possibly length zero).
#'
#' @examples
#' \dontrun{
#' # 1) By group:
#' ord <- mds_order(
#'   x   = grp,        # hsmap_group
#'   tpt = tpt,        # HSMap.tpt
#'   lg  = "all",
#'   drop_outside_envelope = TRUE,
#'   envelope_k   = 2,
#'   envelope_span= 0.3
#' )
#' # Get LG3 order and see dropped markers
#' ord$LG3
#' attr(ord$LG3, "removed")
#'
#' # 2) For a specific set of markers (ignores `lg`):
#' subset_ids <- sample(tpt$markers, 150)
#' ord2 <- mds_order(x = subset_ids, tpt = tpt, plot_each = FALSE)
#' }
#'
#' @importFrom smacof smacofSym
#' @importFrom princurve principal_curve
#' @importFrom graphics polygon
#' @export
mds_order <- function(x,
                            tpt,
                            lg = "all",
                            plot_each = TRUE,
                            itmax = 100000,
                            drop_outside_envelope = FALSE,
                            envelope_k = 2,
                            envelope_span = 0.25,
                            min_keep = 3) {
  # dependencies
  if (!requireNamespace("smacof", quietly = TRUE) ||
      !requireNamespace("princurve", quietly = TRUE)) {
    stop("Packages 'smacof' and 'princurve' are required.")
  }
  `%||%` <- function(a, b) if (is.null(a)) b else a

  # --- matrices from tpt (HSMap.tpt or raw list) -----------------------------
  mats <- if (inherits(tpt, "HSMap.tpt")) tpt$fit else tpt
  if (is.null(mats$r) || is.null(mats$lod_r))
    stop("`tpt` must provide matrices `r` and `lod_r` (in `tpt$fit` for HSMap.tpt).")

  R <- mats$r; L <- mats$lod_r
  if (!is.matrix(R) || !is.matrix(L)) stop("`r` and `lod_r` must be matrices.")
  if (!identical(rownames(R), colnames(R)) ||
      !identical(dim(R), dim(L)) ||
      !identical(rownames(R), rownames(L)) ||
      !identical(colnames(R), colnames(L))) {
    stop("`r` and `lod_r` must be square with identical dimnames.")
  }

  # --- figure out groups / ids to process ------------------------------------
  is_group_obj <- inherits(x, "hsmap_group")
  if (is_group_obj) {
    gvec <- x$groups.snp %||% x$groups
    if (is.null(gvec))
      stop("`hsmap_group` must contain a named integer vector `groups` (or `groups.snp`).")
    split_mrk <- split(names(gvec), gvec)
    all_ids <- sort(as.integer(names(split_mrk)))
    use_lg <- if (identical(lg, "all")) all_ids else intersect(as.integer(lg), all_ids)
    if (!length(use_lg)) stop("No valid linkage groups selected via `lg`.")
  } else {
    # x is expected to be a character vector of marker IDs
    if (!is.character(x)) stop("If `x` is not an 'hsmap_group', it must be a character vector of marker IDs.")
    miss <- setdiff(x, rownames(R))
    if (length(miss)) {
      stop("Some supplied marker IDs are not present in `tpt`: e.g., ",
           paste(utils::head(miss, 8), collapse = ", "), if (length(miss) > 8) " ..." else "")
    }
    # emulate a single LG
    split_mrk <- list(`1` = unique(x))
    use_lg <- 1L
  }

  # --- helper: draw envelope polygon -----------------------------------------
  add_pc_envelope <- function(X, pc, k = 2, span = 0.25,
                              col = grDevices::adjustcolor("black", 0.15),
                              border = NA) {
    lam <- pc$lambda
    ord <- order(lam)
    S <- pc$s[ord, , drop = FALSE]
    # residual per point
    r <- sqrt(rowSums((X - pc$s)^2))
    # smooth width along lambda
    fit_w <- stats::loess(r ~ lam, span = span,
                          control = stats::loess.control(surface = "direct"))
    wS <- as.numeric(stats::predict(fit_w, lam[ord]))
    if (anyNA(wS)) wS[is.na(wS)] <- stats::median(r, na.rm = TRUE)
    wS <- pmax(wS, .Machine$double.eps) * k
    # unit normal
    dS <- rbind(S[2, ] - S[1, ],
                S[3:nrow(S), ] - S[1:(nrow(S) - 2), ],
                S[nrow(S), ] - S[nrow(S) - 1, ])
    dS <- dS / sqrt(rowSums(dS^2))
    N  <- cbind(-dS[, 2], dS[, 1])
    upper <- S + N * wS
    lower <- S - N * wS
    poly  <- rbind(upper, lower[nrow(lower):1, , drop = FALSE])
    polygon(poly, col = col, border = border)
    invisible(NULL)
  }

  # --- core for one group -----------------------------------------------------
  mds_one <- function(ids, lg_label = NULL, do_plot = FALSE) {
    ids <- intersect(ids, rownames(R))
    ids <- ids[!duplicated(ids)]
    if (length(ids) < 3L) {
      out_names <- ids
      attr(out_names, "removed") <- character(0)
      return(out_names)
    }

    rf.mat  <- R[ids, ids, drop = FALSE]
    lod.mat <- L[ids, ids, drop = FALSE]

    na <- is.na(rf.mat) | is.na(lod.mat)
    rf.mat[na]  <- 1e-7
    lod.mat[na] <- 1e-7
    lod.mat     <- lod.mat^2
    diag(lod.mat) <- NA
    diag(rf.mat)  <- NA

    # distance transform
    imf_h <- function(r) -50 * log(1 - 2 * pmin(pmax(r, 1e-8), 0.499999))
    M <- imf_h(rf.mat) / 100

    fit <- smacof::smacofSym(delta = M, ndim = 2, weightmat = lod.mat, itmax = itmax)
    pc1 <- princurve::principal_curve(fit$conf, maxit = 150, smoother = "smooth_spline")

    ord <- pc1$ord
    locinames <- rownames(rf.mat)
    if (is.null(locinames)) locinames <- colnames(rf.mat)
    if (is.null(locinames)) locinames <- paste0("L", seq_len(nrow(rf.mat)))

    # envelope-based filtering (optional)
    removed <- character(0)
    keep_idx_ord <- seq_along(ord)
    if (isTRUE(drop_outside_envelope)) {
      lam <- pc1$lambda
      resid <- sqrt(rowSums((fit$conf - pc1$s)^2))
      fit_w <- stats::loess(resid ~ lam, span = envelope_span,
                            control = stats::loess.control(surface = "direct"))
      width_hat <- as.numeric(stats::predict(fit_w, lam))
      width_hat[is.na(width_hat)] <- stats::median(resid, na.rm = TRUE)
      thr <- envelope_k * width_hat
      outside <- resid > thr
      keep <- ord[!outside[ord]]
      if (length(keep) >= min_keep) {
        keep_idx_ord <- which(!outside[ord])
        removed <- locinames[ord[outside[ord]]]
      } else {
        warning(lg_label, ": envelope filtering would keep fewer than ", min_keep,
                " markers. Skipping envelope filter for this group.")
      }
    }

    ordered_names <- locinames[ord[keep_idx_ord]]

    if (isTRUE(do_plot)) {
      op <- graphics::par(mar = c(4, 4, 2, 1))
      on.exit(graphics::par(op), add = TRUE)
      graphics::plot(fit$conf, pch = 16, cex = 0.6,
                     xlab = "MDS-1", ylab = "MDS-2",
                     main = sprintf("%s: MDS + principal curve", lg_label %||% "Markers"))
      add_pc_envelope(fit$conf, pc1, k = envelope_k, span = envelope_span,
                      col = grDevices::adjustcolor("black", 0.12), border = NA)
      s_ord <- pc1$s[ord, , drop = FALSE]
      graphics::lines(s_ord, lwd = 2)
      if (length(removed)) {
        idx_removed <- match(removed, locinames)
        graphics::points(fit$conf[idx_removed, , drop = FALSE], pch = 16, cex = 0.7,
                         col = grDevices::adjustcolor("red", 0.6))
      }
    }

    attr(ordered_names, "removed") <- removed
    ordered_names
  }

  # --- run over groups (or single set) ---------------------------------------
  if (is_group_obj) {
    res <- vector("list", length(use_lg))
    names(res) <- paste0("LG", use_lg)
    for (i in seq_along(use_lg)) {
      g <- use_lg[i]
      ids <- split_mrk[[as.character(g)]]
      res[[i]] <- mds_one(ids, lg_label = paste0("LG", g), do_plot = isTRUE(plot_each))
    }
  } else {
    res <- list(LG1 = mds_one(split_mrk[[1]], lg_label = "LG1", do_plot = isTRUE(plot_each)))
  }

  class(res) <- "hsmap_mds"
  res
}
