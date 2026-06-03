#' Palette for HSMap dendrograms
#' @param n Positive integer number of groups
#' @return Character vector of length `n` with hex colors
#' @importFrom grDevices colorRampPalette
#' @export
hs_pal <- function(n) {
  base <- c("#ffe119", "#f58231", "#e6194b", "#800000", "#911eb4",
            "#0202a1", "#4363d8", "#42d4f4", "#469990", "#3cb44b")
  if (n == 2)       return(base[c(3, 7)])
  else if (n == 3)  return(base[c(3, 7, 10)])
  else if (n == 4)  return(base[c(3, 7, 10, 5)])
  else if (n == 5)  return(base[c(5, 8, 10, 1, 3)])
  else if (n == 6)  return(base[c(3, 2, 1, 10, 9, 6)])
  colorRampPalette(base)(n)
}

#' Aggregate a square matrix by block-averaging
#'
#' @description
#' Downsample a square matrix by averaging values in non-overlapping
#' \code{fact x fact} blocks. This is useful for speeding up plotting of
#' very large RF/LOD matrices (e.g., in \code{\link{pairwise_heatmap}}).
#'
#' @param M A numeric \emph{square} matrix to be aggregated.
#' @param fact Positive integer factor by which to reduce the resolution.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{R}: Aggregated matrix of size \code{ceiling(n/fact)}.
#'   \item \code{markers}: A named list that maps each aggregated row/column
#'         group (e.g. \code{"Group_1"}) to the original marker IDs contained
#'         in that block.
#' }
#'
#' @examples
#' \dontrun{
#' agg <- aggregate_matrix(M = fit$r, fact = 4)
#' str(agg$markers$Group_3)  # which original markers were in block 3?
#' }
#' @export
aggregate_matrix <- function(M, fact) {
  if (!is.matrix(M)) stop("`M` must be a matrix.")
  if (!is.numeric(fact) || length(fact) != 1L || is.na(fact) || fact <= 0) {
    stop("`fact` must be a positive integer.")
  }
  n <- ncol(M)
  if (n != nrow(M)) stop("`M` must be square.")
  if (is.null(colnames(M)) || is.null(rownames(M))) {
    colnames(M) <- rownames(M) <- seq_len(n)
  }

  # define block index ranges
  id_starts <- seq(1L, n, by = fact)
  id_ends   <- c(id_starts[-1] - 1L, n)
  id        <- cbind(id_starts, id_ends)

  # build name mapping
  markers_list <- vector("list", nrow(id))
  group_names  <- sprintf("Group_%d", seq_len(nrow(id)))
  for (k in seq_len(nrow(id))) {
    idx <- id[k, 1]:id[k, 2]
    markers_list[[k]] <- colnames(M)[idx]
  }
  names(markers_list) <- group_names

  # aggregate by block means
  R <- matrix(NA_real_, nrow(id), nrow(id))
  for (i in seq_len(nrow(id))) {
    for (j in i:nrow(id)) {
      subM <- M[id[i, 1]:id[i, 2], id[j, 1]:id[j, 2], drop = FALSE]
      R[i, j] <- R[j, i] <- mean(subM, na.rm = TRUE)
    }
  }
  rownames(R) <- colnames(R) <- group_names

  list(R = R, markers = markers_list)
}
