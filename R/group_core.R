#' Group markers into linkage groups from two-point RF results (HSMap)
#'
#' @description
#' Clusters markers into **linkage groups (LGs)** using pairwise recombination
#' fractions from HSMap’s two-point analysis. You can pass either:
#' - the raw list returned by `"cpp_pairwise_rf()"` (must contain square matrices
#'   `r` and `lod_r` with matching dimnames), or
#' - an HSMap two-point wrapper result (class `"HSMap.tpt"`, whose `$fit` field
#'   is that same list).
#'
#' The function first symmetrizes and cleans the RF matrix, converts it to a
#' distance (`as.dist` on the symmetrized `r`), then performs hierarchical
#' clustering with `stats::hclust`. Optionally, it uses **LOD-derived weights**
#' to stabilize clustering if you have uneven information content across markers.
#'
#' If `inter = TRUE` in an interactive session and `k` is not specified,
#' you’ll get an iterative **dendrogram preview** to choose `k` (number of groups).
#'
#' @param w Either:
#'   - a **list** with components `r` and `lod_r`, both square numeric matrices
#'     with identical row/column names; or
#'   - an object of class `"HSMap.tpt"` whose `$fit` contains such a list.
#' @param k Integer. Desired number of linkage groups. If `NULL` and `inter = TRUE`,
#'   you can **choose interactively**; otherwise an error is thrown.
#' @param inter Logical. If `TRUE` and running interactively, show a dendrogram
#'   preview (using [plot.hsmap_group()]) and prompt for `k`. Default `TRUE`.
#' @param compare Logical. If `TRUE` and `chrom` is provided, attach a **comparison
#'   matrix** (LG × chromosome) to the result (`$comp`). Default `FALSE`.
#' @param weight_by_lod Logical. If `TRUE`, compute per-marker weights
#'   \eqn{w_i = \mathrm{round}(\mathrm{mean}_j[\mathrm{LOD}_{ij}^2])} and pass them
#'   as `members` to `hclust`. This tends to boost well-supported markers and can
#'   improve cluster stability. If you use weights and `method = "average"`,
#'   we switch to `method = "ward.D2"` (safer with `members`). Default `FALSE`.
#' @param chrom Optional **named** character vector giving chromosome labels per
#'   marker (names must match the matrix dimnames). Only used if `compare = TRUE`.
#' @param method Character linkage method for `stats::hclust`; one of
#'   `"average"`, `"complete"`, `"single"`, `"ward.D2"`, etc. Default `"average"`.
#' @param na_r Numeric in \[0, 0.5\]. Value used to fill `NA` in `r` before
#'   clustering. Default `0.5` (uninformative / unlinked).
#' @param palette A function `f(n)` returning `n` colors to color dendrogram branches
#'   (used by [plot.hsmap_group()]). Default: [hs_pal()].
#'
#' @details
#' ## What the algorithm does
#' Two-point analysis yields \eqn{r_{ij} \in [0, 0.5]} and typically a LOD score
#' for linkage vs \eqn{r = 0.5}. We:
#'
#' 1. **Symmetrize** \eqn{r} by \eqn{\tilde r = (r + r^\top)/2}, set the diagonal to 0,
#'    clamp to \[0, 0.5\], and replace `NA` by `na_r` (default 0.5).
#' 2. Treat \eqn{\tilde r} as a dissimilarity and feed `as.dist(\tilde r)` to
#'    `stats::hclust(method = ...)`.
#' 3. Optionally apply **LOD-based weights**: for marker \eqn{i}, compute
#'    \eqn{w_i = \mathrm{round}(\mathrm{mean}_j[\mathrm{LOD}_{ij}^2])} and pass
#'    `members = w` to `hclust`. With `members`, `"ward.D2"` is typically more
#'    appropriate than `"average"`; we auto-switch in that case.
#' 4. **Cut** the tree at `k` to obtain LG membership.
#'
#' ## Choosing `k`
#' - If `k` is known (e.g., number of chromosomes), provide it directly.
#' - If not, set `inter = TRUE` and pick *visually*: we’ll color the dendrogram
#'   by group and draw rectangles (mappoly2-style).
#'
#' ## Output object
#' The return is an `"hsmap_group"` object with the hclust tree, final `k`,
#' per-marker group calls, (optionally) the comparison matrix, and metadata. It
#' has friendly `print()` and `plot()` methods (see **Value**).
#'
#' @return An object of class `"hsmap_group"` with fields:
#' \describe{
#'   \item{hc}{An `hclust` tree built from the cleaned RF matrix.}
#'   \item{k}{Integer; number of groups used in the final cut.}
#'   \item{groups}{Integer vector, names = marker IDs, giving LG membership.}
#'   \item{markers}{Character; marker IDs used.}
#'   \item{method}{Clustering linkage method.}
#'   \item{comp}{If `compare = TRUE` and `chrom` was given, a LG × chromosome table.}
#' }
#'
#' @seealso `"cpp_pairwise_rf()"` for producing two-point RF/LOD matrices;
#'   [plot.hsmap_group()] and [print.hsmap_group()].
#'
#' @examples
#' \dontrun{
#' # Suppose `w` is the list returned by cpp_pairwise_rf()
#' g <- group_markers(w, k = 10)
#' print(g)
#' plot(g)                 # colored dendrogram w/ rectangles
#'
#' # Interactive choice of k
#' g2 <- group_markers(w, inter = TRUE)  # choose k in the prompt
#'
#' # With chromosome labels, get a comparison matrix
#' chr <- setNames(sim$truth$map$chrom, sim$truth$map$marker_id)
#' g3  <- group_markers(w, k = 10, compare = TRUE, chrom = chr)
#' g3$comp
#' }
#' @export
group_markers <- function(
    w,
    k = NULL,
    inter = TRUE,
    compare = FALSE,
    weight_by_lod = FALSE,
    chrom = NULL,
    method = "average",
    na_r = 0.5,
    palette = hs_pal
) {
  `%||%` <- function(a, b) if (is.null(a)) b else a

  # -- unwrap supported input types -------------------------------------------
  src <- if (inherits(w, "HSMap.tpt")) {
    if (is.null(w$fit) || !is.list(w$fit))
      stop("`HSMap.tpt` object lacks `$fit` (the C++ two-point result).")
    w$fit
  } else if (is.list(w)) {
    w
  } else {
    stop("`w` must be either a list with matrices `r` and `lod_r`, ",
         "or an object of class 'HSMap.tpt'.")
  }

  if (is.null(src$r) || is.null(src$lod_r))
    stop("Input must contain matrices `r` and `lod_r`.")

  Rmat <- src$r
  LOD  <- src$lod_r
  if (!is.matrix(Rmat) || !is.matrix(LOD))
    stop("`r` and `lod_r` must be matrices.")
  if (is.null(colnames(Rmat)) || is.null(rownames(Rmat)) ||
      !identical(colnames(Rmat), rownames(Rmat)))
    stop("`r` must be a square matrix with matching row/col names.")
  if (!identical(dim(Rmat), dim(LOD)) ||
      !identical(rownames(Rmat), rownames(LOD)) ||
      !identical(colnames(Rmat), colnames(LOD)))
    stop("`r` and `lod_r` must have identical dimensions and dimnames.")

  markers <- colnames(Rmat)

  # -- clean & symmetrize r ----------------------------------------------------
  Rsym <- 0.5 * (Rmat + t(Rmat))
  diag(Rsym) <- 0
  Rsym[is.na(Rsym)] <- na_r
  Rsym[Rsym < 0]   <- 0
  Rsym[Rsym > 0.5] <- 0.5

  # -- optional LOD weights ----------------------------------------------------
  members <- NULL
  hclust_method <- method
  if (isTRUE(weight_by_lod)) {
    wts <- rowMeans(LOD^2, na.rm = TRUE)
    wts[!is.finite(wts)] <- 0
    members <- pmax(1L, as.integer(round(wts)))
    if (identical(method, "average")) hclust_method <- "ward.D2"
  }

  # -- clustering --------------------------------------------------------------
  d  <- stats::as.dist(Rsym)
  hc <- stats::hclust(d, method = hclust_method, members = members)

  # -- interactive pick of k ---------------------------------------------------
  if (interactive() && isTRUE(inter) && is.null(k)) {
    message("Interactive dendrogram preview: enter a positive integer for k; ",
            "press Enter or 'Y' to accept.")
    repeat {
      # propose a default (rough): number of large branches at ~25% height
      k_proposed <- k %||% max(2L, length(stats::cutree(hc, h = 0.125)))
      k <- suppressWarnings(as.integer(
        readline(paste0("Enter k [", k_proposed, "]: "))
      ))
      if (is.na(k) || k < 1L) k <- k_proposed

      tmp <- list(hc = hc, k = k, groups = stats::cutree(hc, k = k),
                  markers = markers, method = hclust_method)
      class(tmp) <- "hsmap_group"
      plot(tmp, palette = palette)

      ans <- readline("Accept this k? [Y to accept, or type a new number]: ")
      if (!nzchar(ans) || tolower(substr(ans, 1, 1)) == "y") break
      newk <- suppressWarnings(as.integer(ans))
      k <- if (!is.na(newk) && newk > 0L) newk else NULL
    }
  }
  if (is.null(k)) stop("Provide `k`, or set `inter = TRUE` and choose interactively.")

  # -- final cut & optional comparison ----------------------------------------
  groups <- stats::cutree(hc, k = k)
  names(groups) <- markers

  comp <- NULL
  if (isTRUE(compare)) {
    if (is.null(chrom)) {
      warning("compare = TRUE but `chrom` was not provided; skipping comparison.")
    } else {
      if (is.null(names(chrom)))
        stop("`chrom` must be a *named* vector: names = marker IDs.")
      chr <- chrom[markers]
      chr[is.na(chr)] <- "NoChr"
      tab <- table(
        factor(groups, levels = seq_len(k)),
        factor(chr, levels = unique(chr))
      )
      tab <- cbind(tab, Total = rowSums(tab))
      tab <- rbind(tab, Total = colSums(tab))
      comp <- tab
    }
  }

  out <- list(
    hc      = hc,
    k       = as.integer(k),
    groups  = groups,
    markers = markers,
    method  = hclust_method,
    comp    = comp
  )
  class(out) <- "hsmap_group"
  out
}

