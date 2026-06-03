#' Simulate multi-population data across multiple chromosomes
#'
#' @description
#' Genome-aware wrapper around \code{sim_multi_pop()} that accepts a multi-chromosome
#' map (from \code{make_map()}) and simulates across all markers in order. Between
#' chromosomes, the recombination fraction is fixed to 0.5 so blocks are independent.
#'
#' @param map A data frame returned by \code{make_map()} (or equivalent) containing
#'   at least columns \code{marker_id}, \code{chrom}, \code{position}, and the
#'   attribute \code{"r_vec"} (a length-\code{T-1} numeric vector of recombination
#'   fractions across the genome).
#' @inheritParams sim_multi_pop
#'
#' @return An object of class \code{"sim_multi_pop"} identical in structure to
#'   \code{sim_multi_pop()}, with the truth bundle augmented by \code{chrom} and
#'   \code{position} per marker:
#'   \itemize{
#'     \item \code{truth$map}: the input \code{map} with the same row order as \code{markers_union}.
#'   }
#' @examples
#' \dontrun{
#' map <- make_map(2, c(10, 12), r_const_per_chrom = c(0.02, 0.03))
#' sim <- sim_multi_chrom(
#'   map, n_pops = 2, n_ind_per_pop = c(60, 50),
#'   maternal_geno_mode = "HWE", maternal_pA = 0.5,
#'   paternal_pA_base = 0.4, error_rate = 0.01, seed = 7
#' )
#' }
#' @export
sim_multi_chrom <- function(
    map,
    n_pops = 2,
    n_ind_per_pop = rep(50, n_pops),
    marker_intersection = 1,
    phase_mode = c("all_coupling","random","vector"),
    repulsion_rate = 0.2,
    phase_vector = NULL,
    maternal_geno_mode = c("HWE","all_het","vector"),
    maternal_M_given = NULL,
    maternal_pA = 0.5,
    paternal_pA_base = 0.5,
    paternal_pA_sd = 0,
    known_sire_prob_per_pop = 0.0,
    error_rate = 0.0,
    pedigree = NULL,
    seed = NULL,
    keep_paths = FALSE,
    miss_rate = miss_rate
){
  stopifnot(is.data.frame(map), all(c("marker_id","chrom","position") %in% names(map)))
  markers <- as.character(map$marker_id)
  Tm <- length(markers)
  if (Tm < 2) stop("Map must contain at least 2 markers in total.")
  r_vec <- attr(map, "r_vec")
  if (is.null(r_vec) || length(r_vec) != (Tm - 1L)) {
    stop("`map` must carry attribute 'r_vec' of length T-1 (see make_map()).")
  }

  # Temporarily override the default marker naming inside sim_multi_pop()
  # by setting the union markers to map$marker_id after the call.
  sim <- sim_multi_pop(
    T_markers = Tm,
    n_pops = n_pops,
    n_ind_per_pop = n_ind_per_pop,
    marker_intersection = marker_intersection,
    r_vec = r_vec,
    phase_mode = match.arg(phase_mode),
    repulsion_rate = repulsion_rate,
    phase_vector = phase_vector,
    maternal_geno_mode = match.arg(maternal_geno_mode),
    maternal_M_given = maternal_M_given,
    maternal_pA = maternal_pA,
    paternal_pA_base = paternal_pA_base,
    paternal_pA_sd = paternal_pA_sd,
    known_sire_prob_per_pop = known_sire_prob_per_pop,
    error_rate = error_rate,
    pedigree = pedigree,
    seed = seed,
    keep_paths = keep_paths,
    miss_rate = miss_rate
  )

  # Rename columns/entries to use the map's marker order
  # (sim_multi_pop already uses a union naming; we align to map order)
  for (g in seq_along(sim$G_list)) {
    colnames(sim$G_list[[g]]) <- markers
    sim$M_list[[g]] <- stats::setNames(sim$M_list[[g]], markers)
    colnames(sim$pi_true_list[[g]]) <- markers
    colnames(sim$pi_prior_list[[g]]) <- markers
    colnames(sim$pi_fixed_list[[g]]) <- markers
    names(sim$father_geno_list[[g]]) <- markers
    if (!is.null(sim$H_paths_list)) colnames(sim$H_paths_list[[g]]) <- markers
  }
  sim$truth$markers_union <- markers
  sim$truth$map <- map
  sim
}

#' Construct a simple multi-chromosome map
#'
#' @description
#' Convenience generator for a per-marker map that spans multiple chromosomes.
#' You give the number of chromosomes and the number of markers on each; the
#' function returns a data frame with marker IDs, chromosome labels, positions,
#' and the within-chromosome recombination fraction used between adjacent markers.
#' Recombination between chromosomes is **fixed to 0.5**.
#'
#' @param n_chrom Integer (>=1). Number of chromosomes.
#' @param markers_per_chrom Integer vector length \code{n_chrom}; markers on each chromosome.
#' @param chrom_prefix Character prefix used to label chromosomes. Default \code{"chr"}.
#' @param spacing_bp Numeric; base-pair spacing between adjacent markers on a chromosome
#'   (positions are \code{1, 1+spacing_bp, 1+2*spacing_bp, ...}). Default \code{1e6}.
#' @param r_const_per_chrom Numeric scalar or length-\code{n_chrom} vector; the within-chromosome
#'   recombination fraction used between adjacent markers for each chromosome (passed downstream).
#'
#' @return A data.frame with columns:
#' \itemize{
#'   \item \code{marker_id}: unique marker names (e.g., \code{"chr1_m1"}).
#'   \item \code{chrom}: chromosome label (e.g., \code{"chr1"}).
#'   \item \code{position}: integer base-pair position (monotone within chromosome).
#'   \item \code{r_within}: within-chromosome recombination fraction to the **next** marker
#'         on the same chromosome; \code{NA} for the last marker of each chromosome.
#' }
#' Additionally, the data frame carries an attribute \code{"r_vec"} which is
#' the concatenated \code{T-1} vector of recombination fractions across the genome:
#' within-chromosome \code{r_within}, and \code{0.5} at chromosome boundaries.
#'
#' @examples
#' map <- make_map(n_chrom = 3, markers_per_chrom = c(10, 8, 12), r_const_per_chrom = 0.02)
#' head(map); attr(map, "r_vec")[1:15]
#' @export
make_map <- function(n_chrom,
                     markers_per_chrom,
                     chrom_prefix = "chr",
                     spacing_bp = 1e6,
                     r_const_per_chrom = 0.02) {
  stopifnot(n_chrom >= 1L, length(markers_per_chrom) == n_chrom)
  if (length(r_const_per_chrom) == 1L) {
    r_const_per_chrom <- rep(r_const_per_chrom, n_chrom)
  }
  stopifnot(length(r_const_per_chrom) == n_chrom)

  out_list <- vector("list", n_chrom)
  r_vec_all <- numeric(0)

  for (c in seq_len(n_chrom)) {
    k <- as.integer(markers_per_chrom[c])
    if (k < 1L) stop("All chromosomes must have at least 1 marker.")
    chrom <- paste0(chrom_prefix, c)
    marker_id <- paste0(chrom, "_m", seq_len(k))
    position  <- 1L + spacing_bp * (seq_len(k) - 1L)
    r_within  <- rep(r_const_per_chrom[c], k)
    r_within[k] <- NA_real_  # last marker on the chrom has no within-chrom next

    dfc <- data.frame(
      marker_id = marker_id,
      chrom     = chrom,
      position  = as.integer(position),
      r_within  = r_within,
      stringsAsFactors = FALSE
    )
    out_list[[c]] <- dfc

    # build genome-wide r_vec: within-chrom r, and a 0.5 breakpoint between chroms
    if (k > 1L) r_vec_all <- c(r_vec_all, rep(r_const_per_chrom[c], k - 1L))
    if (c < n_chrom) r_vec_all <- c(r_vec_all, 0.5)  # chrom boundary
  }

  map <- do.call(rbind, out_list)
  attr(map, "r_vec") <- r_vec_all
  map
}
