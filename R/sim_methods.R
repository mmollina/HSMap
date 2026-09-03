
#' Print method for simulated multi-population datasets
#'
#' @description
#' Stylish, compact summary for objects of class \code{"sim_multi_pop"} produced by
#' \code{\link{sim_multi_pop}}.
#' The printer reports:
#' \itemize{
#'   \item union marker count and number of intervals,
#'   \item number of populations and offspring per population,
#'   \item brief summary of the recombination vector (if available),
#'   \item per-population marker presence, missingness, maternal heterozygosity rate,
#'         and whether a known sire was simulated.
#' }
#'
#' @param x An object of class \code{"sim_multi_pop"}.
#' @param ... Unused; for S3 compatibility.
#'
#' @return The input \code{x}, invisibly.
#'
#' @examples
#' \dontrun{
#' sim <- sim_multi_pop(100, n_pops = 2, n_ind_per_pop = c(40, 35))
#' print(sim)
#' }
#'
#' @importFrom cli rule cat_line cat_bullet symbol col_green col_cyan col_blue col_yellow col_magenta col_silver
#' @importFrom crayon bold italic %+%
#' @method print sim_multi_pop
#' @export
print.sim_multi_pop <- function(x, ...) {
  stopifnot(is.list(x), inherits(x, "sim_multi_pop"))

  # ---- small helpers --------------------------------------------------------
  nzchar1 <- function(z) !is.null(z) && length(z) > 0
  pct     <- function(x) if (is.na(x)) "NA" else sprintf("%.1f%%", 100 * x)
  pnum    <- function(x) if (is.na(x)) "NA" else sprintf("%.4f", x)
  span_str <- function(ids) {
    k <- length(ids)
    if (k == 0) return("(none)")
    if (k <= 8) paste(ids, collapse = ", ")
    else paste0(ids[1], ", ", ids[2], ", ", cli::symbol$ellipsis, ", ",
                ids[k-1], ", ", ids[k])
  }
  safe_colnames <- function(m) { cn <- tryCatch(colnames(m), error = function(e) NULL); cn %||% character() }
  `%||%` <- function(a, b) if (!is.null(a)) a else b

  # ---- union markers / intervals -------------------------------------------
  markers <- tryCatch(x$truth$markers_union, error = function(e) NULL)
  if (is.null(markers) || !length(markers)) {
    # fallback: infer from first G
    G1 <- x$G_list[[1]]
    markers <- safe_colnames(G1)
  }
  Tm <- length(markers)
  n_int <- max(0L, Tm - 1L)

  # ---- r_true summary (optional) -------------------------------------------
  r_true <- tryCatch(x$truth$r_true, error = function(e) NULL)
  has_r  <- !is.null(r_true) && length(r_true) == n_int

  # ---- populations ----------------------------------------------------------
  pop_names <- names(x$G_list)
  if (is.null(pop_names)) pop_names <- paste0("Pop", seq_along(x$G_list))
  Gn <- length(pop_names)

  # ---- header ---------------------------------------------------------------
  cli::rule(left = crayon::bold("sim_multi_pop"), right = cli::col_silver("synthetic dataset"))
  cli::cat_line(
    paste0(
      cli::col_cyan("Markers: "), crayon::bold(Tm),
      "  ", cli::symbol$line, "  ",
      cli::col_cyan("Populations: "), crayon::bold(Gn)
    )
  )

  # r summary if present
  if (has_r) {
    cli::cat_line(
      paste0(
        cli::col_silver("r_true"),
        ": min=", pnum(min(r_true)),
        "  med=", pnum(stats::median(r_true)),
        "  mean=", pnum(mean(r_true)),
        "  max=", pnum(max(r_true))
      )
    )
  }

  # Marker IDs (truncated)
  if (nzchar1(markers)) {
    cli::cat_line(paste0(cli::col_silver("Markers (IDs)"), ": ", span_str(markers)))
  }

  # ---- per-population summary ----------------------------------------------
  cli::rule(center = cli::col_blue("Per population"))
  for (g in seq_len(Gn)) {
    Gg <- as.matrix(x$G_list[[g]])
    Mg <- x$M_list[[g]]
    n_off <- nrow(Gg) %||% 0L

    # which markers present for this population? (any non-NA genotype)
    present_cols <- integer()
    if (!is.null(Gg) && length(Gg)) {
      present_cols <- which(colSums(!is.na(Gg)) > 0)
    }
    n_present <- length(present_cols)
    miss_rate <- if (n_off > 0 && n_present > 0) {
      mean(is.na(Gg[, present_cols, drop = FALSE]))
    } else NA_real_

    # maternal heterozygosity among present markers
    mhet <- if (!is.null(Mg) && n_present > 0) {
      mean(Mg[present_cols] == 1L, na.rm = TRUE)
    } else NA_real_

    # known sire? (either father_geno_list or pi_fixed_list has one-hot columns)
    known_sire <- FALSE
    known_cols <- 0L
    if (!is.null(x$father_geno_list) && !is.null(x$father_geno_list[[g]])) {
      fg <- x$father_geno_list[[g]]
      known_cols <- sum(!is.na(fg))
      known_sire <- known_cols > 0
    } else if (!is.null(x$pi_fixed_list) && !is.null(x$pi_fixed_list[[g]])) {
      pf <- x$pi_fixed_list[[g]]
      if (!is.null(dim(pf))) known_cols <- sum(colSums(!is.na(pf)) == 3L)
      known_sire <- known_cols > 0
    }

    # one bullet per population
    cli::cat_bullet(
      paste0(
        crayon::bold(pop_names[g]), "  ",
        cli::col_magenta("(offspring: ", n_off, ")"), "  ",
        cli::symbol$line, "  ",
        cli::col_green("present "), n_present, "/", Tm, " (", pct(n_present / max(1, Tm)), ")  ",
        cli::symbol$line, "  ",
        cli::col_yellow("missing "), pct(ifelse(is.na(miss_rate), 0, miss_rate)), "  ",
        cli::symbol$line, "  ",
        cli::col_cyan("M_het "), pct(ifelse(is.na(mhet), 0, mhet)), "  ",
        cli::symbol$line, "  ",
        (if (known_sire) cli::col_green("known sire: yes") else cli::col_silver("known sire: no")),
        if (known_cols > 0) paste0(" (", known_cols, " markers)") else ""
      ),
      bullet = "arrow_right"
    )
  }

  # hidden paths stored?
  if (!is.null(x$H_paths_list) && any(vapply(x$H_paths_list, function(M) !is.null(M), logical(1)))) {
    cli::cat_line(cli::col_blue("Maternal hidden paths stored: ") %+% crayon::bold("YES"))
  }

  invisible(x)
}


#' Plot a simulated multi-population dataset (offspring genotypes and maternal genotypes)
#'
#' @description
#' S3 plotting method for objects of class \code{"sim_multi_pop"} produced by
#' \code{\link{sim_multi_pop}}. It visualizes offspring genotype doses
#' (\code{0/1/2/NA}) as a heatmap across markers, stacked by population, and
#' (optionally) overlays a thin top row per population for the maternal genotype.
#'
#' @param x An object of class \code{"sim_multi_pop"}.
#' @param pops_order Optional character vector to set/override the display order
#'   of populations (must match \code{names(x$G_list)}).
#' @param max_ind_per_pop Positive integer; if finite, limit the number of
#'   offspring rows plotted per population (top rows kept). Default \code{Inf}.
#' @param show_maternal Logical; if \code{TRUE} (default), draws a thin top row
#'   per population with the maternal genotype (labelled \code{"M"}).
#' @param palette Named character vector mapping doses \code{0,1,2} to colors.
#'   Default: blue/orange/red (ColorBrewer-ish).
#' @param na_color Color used for \code{NA} entries (missing markers or calls).
#'   Default \code{"grey90"}.
#' @param x_tick_every Integer; place an x-axis tick/label every \code{x_tick_every}
#'   markers when marker IDs are of the form \code{m1..mT}. Default \code{10}.
#' @param per_pop_panels Logical; if \code{TRUE}, facet by population (one panel
#'   per pop). If \code{FALSE} (default), all populations are stacked in a single
#'   heatmap.
#' @param maternal_label Character label used for the maternal row. Default \code{"M"}.
#' @param ... Unused; for S3 compatibility.
#'
#' @details
#' The method first aligns each population's genotype matrix \code{G} and maternal
#' vector \code{M} to the **union of marker names** (column order preserved). Markers
#' absent in a population are padded with \code{NA}. Rows (offspring) are stacked by
#' population; if \code{show_maternal=TRUE}, a single maternal row is placed at the
#' top of each population block.
#'
#' @return A \code{ggplot} object (invisibly).
#'
#' @examples
#' \dontrun{
#' set.seed(42)
#' sim <- sim_multi_pop(
#'   T_markers = 100, n_pops = 2, n_ind_per_pop = c(60, 45),
#'   marker_intersection = 0.8, r_const = 0.02,
#'   phase_mode = "all_coupling", maternal_pA = 0.5,
#'   paternal_pA_base = 0.4, error_rate = 0.01
#' )
#' plot(sim, max_ind_per_pop = 40, x_tick_every = 10)
#' }
#'
#' @importFrom ggplot2 ggplot geom_tile scale_fill_manual scale_x_discrete labs
#' @importFrom ggplot2 theme_minimal theme element_blank element_text aes facet_grid
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr bind_rows
#' @importFrom tidyselect all_of
#' @method plot sim_multi_pop
#' @export
plot.sim_multi_pop <- function(
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
  G_list <- x$G_list
  M_list <- x$M_list

  # ---- union of markers, preserve numeric order if IDs look like m1..mT ----
  all_markers <- unique(unlist(lapply(G_list, function(G) colnames(as.matrix(G)))))
  if (is.null(all_markers) || any(is.na(all_markers)))
    stop("All G matrices must have column names (marker IDs).")

  mk_num <- suppressWarnings(as.integer(gsub("[^0-9]", "", all_markers)))
  if (!any(is.na(mk_num))) {
    all_markers <- all_markers[order(mk_num)]
  }

  # ---- population order ----
  pop_names <- names(G_list)
  if (is.null(pop_names)) pop_names <- paste0("Pop", seq_along(G_list))
  if (!is.null(pops_order)) {
    idx <- match(pops_order, pop_names)
    if (any(is.na(idx))) stop("pops_order contains names not in G_list.")
    G_list <- G_list[idx]
    M_list <- M_list[idx]
    pop_names <- pops_order
  }

  # ---- alignment helpers ----
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

  # ---- build long tibble for ggplot ----
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

    if (show_maternal && !is.null(M_list)) {
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

  # order y so first population appears at top (ggplot draws bottom->top)
  y_levels_rev <- rev(unique(y_levels))
  df_all$ykey   <- factor(df_all$ykey,   levels = y_levels_rev)
  df_all$marker <- factor(df_all$marker, levels = all_markers)
  df_all$dose_f <- factor(df_all$dose,   levels = c(0,1,2), ordered = TRUE)

  # x-axis breaks (every x_tick_every markers if numeric suffix)
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
    ggplot2::geom_tile(height = 0.95, width = 0.95, color = NA) +
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
