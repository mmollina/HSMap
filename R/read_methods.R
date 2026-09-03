#' @export
print.HSMap.data <- function(x, ...) {
  stopifnot(inherits(x, "HSMap.data"))
  markers <- tryCatch(x$truth$markers_union, error = function(e) NULL)
  if (is.null(markers)) {
    # infer from first G
    if (length(x$G_list) && nrow(x$G_list[[1]]) >= 0) {
      markers <- colnames(x$G_list[[1]])
    }
  }
  Tm <- length(markers)
  Gn <- length(x$G_list)

  cat("HSMap.data\n")
  cat("  Markers     :", Tm, "\n")
  cat("  Populations :", Gn, "\n")

  if (!is.null(x$alleles)) {
    has_map <- all(c("marker_id","REF","ALT") %in% names(x$alleles))
    if (has_map) {
      cat("  Alleles     : ", nrow(x$alleles), " rows (",
          paste(intersect(c("marker_id","REF","ALT","chrom","position"), names(x$alleles)), collapse = ", "),
          ")\n", sep = "")
    }
  }

  # show brief per-pop stats if present
  if (!is.null(x$stats) && nrow(x$stats)) {
    cat("\nPer-population summary (first 6 rows):\n")
    print(utils::head(x$stats, 6))
  }
  invisible(x)
}

#' @export
plot.HSMap.data <- function(
    x,
    pops_order        = NULL,
    max_ind_per_pop   = Inf,
    show_maternal     = TRUE,
    palette           = c(`0` = "#1f78b4", `1` = "#fdbf6f", `2` = "#e31a1c"),
    na_color          = "grey90",
    x_tick_every      = 10,
    per_pop_panels    = FALSE,
    maternal_label    = "M",
    ...
) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Install 'ggplot2' to use plot.HSMap.data().")
  }
  if (!requireNamespace("tidyr", quietly = TRUE) ||
      !requireNamespace("dplyr", quietly = TRUE) ||
      !requireNamespace("tidyselect", quietly = TRUE)) {
    stop("Install 'tidyr', 'dplyr', and 'tidyselect' to use plot.HSMap.data().")
  }

  G_list <- x$G_list
  M_list <- x$M_list

  # union of markers
  all_markers <- tryCatch(x$truth$markers_union, error = function(e) NULL)
  if (is.null(all_markers) || !length(all_markers)) {
    all_markers <- unique(unlist(lapply(G_list, function(G) colnames(as.matrix(G)))))
  }
  if (is.null(all_markers) || any(is.na(all_markers)))
    stop("All genotype matrices must have column names (marker IDs).")

  # try numeric ordering m1..mT if applicable
  mk_num <- suppressWarnings(as.integer(gsub("[^0-9]", "", all_markers)))
  if (!any(is.na(mk_num))) all_markers <- all_markers[order(mk_num)]

  # population order
  pop_names <- names(G_list)
  if (is.null(pop_names)) pop_names <- paste0("Pop", seq_along(G_list))
  if (!is.null(pops_order)) {
    idx <- match(pops_order, pop_names)
    if (any(is.na(idx))) stop("pops_order contains names not in G_list.")
    G_list <- G_list[idx]
    M_list <- M_list[idx]
    pop_names <- pops_order
  }

  # alignment helpers
  align_G <- function(G) {
    G <- as.matrix(G)
    if (is.null(rownames(G))) rownames(G) <- paste0("O", seq_len(nrow(G)))
    out <- matrix(NA_integer_, nrow = nrow(G), ncol = length(all_markers),
                  dimnames = list(rownames(G), all_markers))
    common <- intersect(colnames(G), all_markers)
    if (length(common) > 0) out[, common] <- G[, common, drop = FALSE]
    out
  }
  align_M <- function(M) {
    out <- rep(NA_integer_, length(all_markers)); names(out) <- all_markers
    if (!is.null(M)) {
      common <- intersect(names(M), all_markers)
      if (length(common) > 0) out[common] <- as.integer(M[common])
    }
    out
  }

  # build long df
  dfs <- list()
  y_levels <- character(0)
  for (g in seq_along(G_list)) {
    pop <- pop_names[g]
    Gg  <- align_G(G_list[[g]])
    if (is.finite(max_ind_per_pop)) {
      keep <- seq_len(min(nrow(Gg), max_ind_per_pop))
      Gg <- Gg[keep, , drop = FALSE]
    }
    df_g <- as.data.frame(Gg, stringsAsFactors = FALSE)
    df_g$ind <- rownames(Gg)
    df_g$pop <- pop
    df_g <- tidyr::pivot_longer(df_g, cols = tidyselect::all_of(all_markers),
                                names_to = "marker", values_to = "dose")

    if (isTRUE(show_maternal) && !is.null(M_list) && length(M_list) >= g) {
      Mg <- align_M(M_list[[g]])
      df_m <- data.frame(ind = maternal_label, pop = pop,
                         marker = all_markers, dose = as.integer(Mg),
                         stringsAsFactors = FALSE)
      df  <- dplyr::bind_rows(df_m, df_g)
      y_levels <- c(y_levels, paste0(pop, ":", c(maternal_label, rownames(Gg))))
    } else {
      df  <- df_g
      y_levels <- c(y_levels, paste0(pop, ":", rownames(Gg)))
    }
    df$ykey <- paste0(df$pop, ":", df$ind)
    dfs[[g]] <- df
  }
  df_all <- dplyr::bind_rows(dfs)

  # factor levels for order
  y_levels_rev <- rev(unique(y_levels))
  df_all$ykey   <- factor(df_all$ykey,   levels = y_levels_rev)
  df_all$marker <- factor(df_all$marker, levels = all_markers)
  df_all$dose_f <- factor(df_all$dose,   levels = c(0,1,2), ordered = TRUE)

  # x-axis breaks
  mk_num2 <- suppressWarnings(as.integer(gsub("[^0-9]", "", all_markers)))
  if (!any(is.na(mk_num2))) {
    brk_idx  <- which(mk_num2 %% x_tick_every == 0)
    x_breaks <- all_markers[brk_idx]
    x_labels <- mk_num2[brk_idx]
  } else {
    k <- length(all_markers)
    brk_idx  <- unique(round(seq(1, k, length.out = min(10, k))))
    x_breaks <- all_markers[brk_idx]
    x_labels <- x_breaks
  }

  p <- ggplot2::ggplot(df_all, ggplot2::aes(x = marker, y = ykey, fill = dose_f)) +
    ggplot2::geom_tile(height = 0.95, width = 0.95) +
    ggplot2::scale_fill_manual(values = palette, na.value = na_color, name = "Dose") +
    ggplot2::scale_x_discrete(breaks = x_breaks, labels = x_labels, expand = c(0,0)) +
    ggplot2::labs(x = "Marker", y = NULL) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1),
      legend.position = "right"
    )

  if (isTRUE(per_pop_panels)) {
    p <- p + ggplot2::facet_grid(rows = ggplot2::vars(pop), scales = "free_y", space = "free_y")
  }

  print(p)
  invisible(p)
}
