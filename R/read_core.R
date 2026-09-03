#' Read HSMap pedigree and genotype CSVs
#'
#' @param pedigree Path to pedigree CSV (columns: id, mother, father, generation, family_id).
#' @param genotypes Path to genotype CSV (columns: marker_id, REF, ALT, chrom, position, + samples).
#' @param na_strings Character vector of values treated as NA.
#'
#' @return An object of class HSMap.data with:
#'   - G_list: list of offspring genotype matrices per family (offspring x markers).
#'   - M_list: list of maternal genotype vectors per family.
#'   - alleles: data.frame (marker_id, REF, ALT, chrom, position).
#'   - pedigree: the pedigree data frame.
#'   - stats: per-family summary statistics.
#'
#' @examples
#' # A small simulated open-pollinated example dataset ships with the package.
#' ped  <- system.file("extdata", "example_pedigree.csv",  package = "HSMap")
#' geno <- system.file("extdata", "example_genotypes.csv", package = "HSMap")
#' dat  <- read_HSMap_data(ped, geno)
#' dat
#' @seealso \code{\link{pairwise_rf}}, \code{\link{hmm_map_blocks}}; the
#'   \code{vignette("getting-started", package = "HSMap")} for the full workflow.
#' @importFrom stats na.omit setNames
#' @export
read_HSMap_data <- function(pedigree, genotypes, na_strings = c("NA",".","")) {
  ped <- utils::read.csv(pedigree, na.strings = na_strings, stringsAsFactors = FALSE)
  geno_df <- utils::read.csv(genotypes, na.strings = na_strings, stringsAsFactors = FALSE)

  # minimal checks
  reqp <- c("id","mother","father","generation","family_id")
  if (!all(reqp %in% names(ped))) {
    stop("`pedigree` must contain columns: ", paste(reqp, collapse = ", "))
  }
  reqg <- c("marker_id","REF","ALT")
  if (!all(reqg %in% names(geno_df))) {
    stop("`genotypes` must contain columns: marker_id, REF, ALT (plus samples).")
  }

  # allele metadata
  alleles <- geno_df[, c("marker_id","REF","ALT",
                         intersect(c("chrom","position"), names(geno_df))), drop = FALSE]
  markers <- alleles$marker_id

  # genotype matrix markers x samples
  sample_cols <- setdiff(names(geno_df), names(alleles))
  X <- as.matrix(geno_df[, sample_cols, drop = FALSE])
  storage.mode(X) <- "integer"
  rownames(X) <- markers
  colnames(X) <- sample_cols

  # build per-family G_list and M_list
  fam_ids <- unique(ped$family_id)
  G_list <- M_list <- setNames(vector("list", length(fam_ids)), fam_ids)
  stats_rows <- vector("list", length(fam_ids))

  for (k in seq_along(fam_ids)) {
    fam <- fam_ids[k]
    fam_ped <- ped[ped$family_id == fam, ]

    mom_id <- unique(na.omit(fam_ped$mother))
    mom_id <- if (length(mom_id) == 1L) mom_id else fam

    # maternal vector
    Mvec <- rep(NA_integer_, length(markers)); names(Mvec) <- markers
    if (mom_id %in% colnames(X)) {
      Mvec <- X[, mom_id]
    }

    # offspring
    kids <- fam_ped$id[fam_ped$generation == 2 & fam_ped$id %in% colnames(X)]
    Gmat <- if (length(kids)) t(X[, kids, drop = FALSE]) else
      matrix(NA_integer_, 0, length(markers), dimnames = list(character(0), markers))

    G_list[[fam]] <- Gmat
    M_list[[fam]] <- Mvec

    # quick stats
    miss_rate <- if (nrow(Gmat)) mean(is.na(Gmat)) else NA_real_
    mhet <- mean(Mvec == 1L, na.rm = TRUE)
    stats_rows[[k]] <- data.frame(
      family_id = fam,
      mother_id = mom_id,
      n_offspring = nrow(Gmat),
      n_markers = length(markers),
      missing_rate = miss_rate,
      maternal_het_rate = if (is.finite(mhet)) mhet else NA_real_,
      stringsAsFactors = FALSE
    )
  }
  stats_df <- do.call(rbind, stats_rows)

  out <- list(
    G_list = G_list,
    M_list = M_list,
    alleles = alleles,
    pedigree = ped,
    stats = stats_df
  )
  class(out) <- "HSMap.data"
  out
}
