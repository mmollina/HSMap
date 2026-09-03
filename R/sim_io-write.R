#' Write a wide genotype table (parents first, then offspring)
#'
#' @description
#' Produces (or writes) a marker-by-sample table with these columns:
#' \code{marker_id, REF, ALT, chrom, position, [one column per individual]}.
#' Parents (mothers, then sires if present) are placed first, followed by all offspring.
#' Genotypes are **doses** \code{0/1/2/NA}.
#'
#' @param sim Object returned by \code{sim_multi_pop()} or \code{sim_multi_chrom()}.
#' @param file Optional path. If provided, the table is written as CSV (UTF-8).
#' @param base_levels Character vector from \code{c("A","C","G","T")}.
#'   Used to randomly assign \code{REF} and \code{ALT}; ALT is sampled from
#'   \code{base_levels[base_levels != REF]}.
#'
#' @return A data.frame with one row per marker and one column per individual
#'   (parents first), in addition to the fixed header columns.
#' @importFrom utils write.csv
#' @export
write_sim_genotypes <- function(sim, file = NULL,
                                base_levels = c("A","C","G","T")) {
  stopifnot(inherits(sim, "sim_multi_pop"))
  if (!all(base_levels %in% c("A","C","G","T"))) {
    stop("`base_levels` must be subset of c('A','C','G','T').")
  }

  # --- marker metadata -------------------------------------------------------
  markers <- sim$truth$markers_union
  if (is.null(markers) || !length(markers)) stop("Simulation object lacks marker names.")
  chrom <- if (!is.null(sim$truth$map)) sim$truth$map$chrom else rep(NA_character_, length(markers))
  pos   <- if (!is.null(sim$truth$map)) sim$truth$map$position else rep(NA_integer_, length(markers))
  names(chrom) <- names(pos) <- markers

  # REF/ALT randomization
  set.seed(sample.int(.Machine$integer.max, 1L)) # harmless variability
  ref <- sample(base_levels, length(markers), replace = TRUE)
  alt <- vapply(ref, function(r) sample(base_levels[base_levels != r], 1L), character(1L))

  # --- collect IDs -----------------------------------------------------------
  pops <- sim$truth$pops_meta
  mother_ids    <- vapply(pops, `[[`, character(1), "mother_id")
  father_ids_all<- vapply(pops, function(p) if (is.na(p$father_id)) NA_character_ else p$father_id, character(1))
  offspring_ids <- unlist(lapply(pops, `[[`, "offspring_ids"), use.names = FALSE)

  # --- parents: each as a length-T vector (markers in rows) ------------------
  build_mom_vec <- function(g) {
    Mg <- sim$M_list[[g]]
    v  <- rep(NA_integer_, length(markers)); names(v) <- markers
    common <- intersect(names(Mg), markers)
    if (length(common)) v[common] <- as.integer(Mg[common])
    v
  }
  mom_cols <- lapply(seq_along(pops), build_mom_vec)
  mom_mat  <- do.call(cbind, mom_cols)
  rownames(mom_mat) <- markers
  colnames(mom_mat) <- mother_ids

  build_dad_vec <- function(g) {
    fg <- sim$father_geno_list[[g]]
    v  <- rep(NA_integer_, length(markers)); names(v) <- markers
    if (!is.null(fg)) {
      names(fg) <- if (is.null(names(fg))) markers else names(fg)
      common <- intersect(names(fg), markers)
      if (length(common)) v[common] <- as.integer(fg[common])
    }
    v
  }
  dad_cols <- lapply(seq_along(pops), build_dad_vec)
  have_dad <- !is.na(father_ids_all)
  dad_cols <- dad_cols[have_dad]
  father_ids <- father_ids_all[have_dad]
  dad_mat  <- if (length(dad_cols)) {
    out <- do.call(cbind, dad_cols)
    rownames(out) <- markers
    colnames(out) <- father_ids
    out
  } else NULL

  # --- offspring: return markers in rows, offspring in columns ---------------
  build_offspring_mat <- function(g) {
    Gg <- as.matrix(sim$G_list[[g]])
    Gout <- matrix(NA_integer_, nrow = length(markers), ncol = nrow(Gg),
                   dimnames = list(markers, rownames(Gg)))
    common <- intersect(colnames(Gg), markers)
    if (length(common)) {
      Gout[common, rownames(Gg)] <- t(Gg[, common, drop = FALSE])
    }
    Gout
  }
  kid_mats <- lapply(seq_along(pops), build_offspring_mat)
  kids_mat <- do.call(cbind, kid_mats)

  # --- bind: all have markers in rows ----------------------------------------
  X <- cbind(mom_mat, dad_mat, kids_mat)
  stopifnot(nrow(X) == length(markers))
  rownames(X) <- markers

  # --- finalize data.frame ---------------------------------------------------
  df <- data.frame(
    marker_id = markers,
    REF = ref,
    ALT = alt,
    chrom = chrom[markers],
    position = as.integer(pos[markers]),
    stringsAsFactors = FALSE
  )
  df <- cbind(df, as.data.frame(X, check.names = FALSE))

  if (!is.null(file)) {
    write.csv(df, file = file, row.names = FALSE, na = "", fileEncoding = "UTF-8")
  }
  df
}

#' Write a pedigree table from a simulation
#'
#' @description
#' Creates a simple pedigree \code{data.frame} (and optionally writes it to disk as CSV)
#' from a \code{"sim_multi_pop"} object: one row per individual with IDs and parents.
#'
#' @param sim Object returned by \code{sim_multi_pop()} or \code{sim_multi_chrom()}.
#' @param file Optional path. If provided, the table is written as CSV (UTF-8).
#'
#' @return A data.frame with columns:
#' \itemize{
#'   \item \code{id}: individual ID,
#'   \item \code{mother}: mother ID (\code{NA} for mothers and sires),
#'   \item \code{father}: father ID (\code{NA} if unknown),
#'   \item \code{generation}: \code{1} for parents, \code{2} for offspring,
#'   \item \code{family_id}: population/family identifier (the mother ID).
#' }
#' @examples
#' \dontrun{
#' ped <- write_sim_pedigree(sim, file = "pedigree.csv")
#' }
#' @export
write_sim_pedigree <- function(sim, file = NULL) {
  stopifnot(inherits(sim, "sim_multi_pop"))
  pops <- sim$truth$pops_meta
  rows <- list()
  for (g in seq_along(pops)) {
    fam <- pops[[g]]
    mom <- fam$mother_id
    dad <- if (is.na(fam$father_id)) NA_character_ else fam$father_id

    # parents
    rows[[length(rows)+1L]] <- data.frame(
      id = mom, mother = NA_character_, father = NA_character_,
      generation = 1L, family_id = mom, stringsAsFactors = FALSE
    )
    if (!is.na(dad)) {
      rows[[length(rows)+1L]] <- data.frame(
        id = dad, mother = NA_character_, father = NA_character_,
        generation = 1L, family_id = mom, stringsAsFactors = FALSE
      )
    }

    # offspring
    if (length(fam$offspring_ids)) {
      rows[[length(rows)+1L]] <- data.frame(
        id = fam$offspring_ids,
        mother = mom,
        father = if (is.na(dad)) NA_character_ else dad,
        generation = 2L,
        family_id = mom,
        stringsAsFactors = FALSE
      )
    }
  }
  ped <- do.call(rbind, rows)
  if (!is.null(file)) {
    write.csv(ped, file = file, row.names = FALSE, na = "", fileEncoding = "UTF-8")
  }
  ped
}
