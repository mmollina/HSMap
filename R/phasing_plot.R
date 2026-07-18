#' Tile plot of adjacent phase by dam (HSMap)
#'
#' @description
#' Visualize adjacent phases (coupling vs repulsion) along a marker order for
#' one dam or many dams. Each tile is an interval between consecutive markers.
#' Color encodes phase (e.g., coupling vs repulsion), and optional alpha encodes
#' LOD support for the adjacent pair from \code{tpt$fit$lod_ph}.
#'
#' @param phased An \code{HSMap.phased} or \code{HSMap.phased.multi} object
#'   returned by \code{phase_from_pairwise()}.
#' @param tpt An \code{HSMap.tpt} with \code{fit$lod_ph} to compute per-interval
#'   LOD weights (optional but recommended).
#' @param palette Named character vector of two colors:
#'   \code{c(repulsion="...", coupling="...")}. Defaults to \code{steelblue}
#'   (repulsion) and \code{tomato} (coupling).
#' @param alpha_by_lod Logical. If \code{TRUE}, tile alpha is proportional to
#'   \code{lod_ph\[ i, i+1 \]} in the given order (rescaled per dam). Default \code{TRUE}.
#' @param show_markers Logical. If \code{TRUE}, show x-axis labels for the interval
#'   endpoints (can be heavy for many markers). Default \code{FALSE}.
#'
#' @return A \code{ggplot} object.
#' @importFrom stats ave
#' @export
plot_phase <- function(
    phased,
    tpt = NULL,
    palette = c(repulsion = "steelblue", coupling = "tomato"),
    alpha_by_lod = TRUE,
    show_markers = FALSE
) {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Package 'ggplot2' is required for this plot.")

  as_list <- if (inherits(phased, "HSMap.phased.multi")) {
    phased
  } else if (inherits(phased, "HSMap.phased")) {
    setNames(list(phased), phased$dam %||% "Dam1")
  } else {
    stop("`phased` must be HSMap.phased or HSMap.phased.multi.")
  }

  # optional LOD matrix
  LOD <- NULL
  if (!is.null(tpt)) {
    if (!inherits(tpt, "HSMap.tpt") || !is.matrix(tpt$fit$lod_ph))
      stop("`tpt` must be HSMap.tpt with fit$lod_ph matrix.")
    LOD <- tpt$fit$lod_ph
  }

  # build long data
  rows <- list()
  for (nm in names(as_list)) {
    ph <- as_list[[nm]]
    ord <- ph$order
    pv  <- ph$phase_vec
    if (length(pv) != length(ord) - 1L) next

    # LOD for adjacent pairs
    adj_lod <- rep(NA_real_, length(pv))
    if (!is.null(LOD)) {
      ok <- ord[ord %in% rownames(LOD)]
      if (length(ok) == length(ord)) {
        for (i in seq_along(pv)) {
          adj_lod[i] <- LOD[ ord[i], ord[i+1] ]
        }
      }
    }
    rows[[nm]] <- data.frame(
      dam     = nm,
      interval= seq_along(pv),
      phase   = ifelse(pv == 1L, "coupling", "repulsion"),
      lod_adj = adj_lod,
      left_marker  = ord[-length(ord)],
      right_marker = ord[-1L],
      stringsAsFactors = FALSE
    )
  }
  df <- do.call(rbind, rows)
  if (!nrow(df)) stop("Nothing to plot.")

  # alpha scaling per dam (robust)
  if (isTRUE(alpha_by_lod) && "lod_adj" %in% names(df)) {
    df$alpha <- ave(df$lod_adj, df$dam, FUN = function(v) {
      v <- pmax(v, 0)
      if (all(is.na(v)) || max(v, na.rm = TRUE) == 0) return(rep(1, length(v)))
      v / max(v, na.rm = TRUE)
    })
  } else {
    df$alpha <- 1
  }

  gg <- ggplot2::ggplot(df, ggplot2::aes(x = interval, y = dam)) +
    ggplot2::geom_tile(ggplot2::aes(fill = phase, alpha = alpha), height = 0.9) +
    ggplot2::scale_fill_manual(values = palette) +
    ggplot2::scale_alpha(range = c(0.3, 1), guide = "none") +
    ggplot2::labs(x = "Interval (between consecutive markers)", y = NULL,
                  fill = "Phase") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(panel.grid = ggplot2::element_blank())

  if (isTRUE(show_markers)) {
    gg <- gg +
      ggplot2::scale_x_continuous(
        breaks = df$interval,
        labels = paste0(df$left_marker, "-", df$right_marker)
      ) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5))
  }

  gg
}

`%||%` <- function(a, b) if (is.null(a)) b else a
