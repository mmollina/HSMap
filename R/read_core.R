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
#'   Cross-aware fields (added in the known-sire extension; safe to ignore for
#'   open-pollinated workflows):
#'   - cross_table: one row per cross with \code{cross_id}, \code{mother_id},
#'     \code{father_id}, \code{family_type} (\code{"open_pollinated"},
#'     \code{"known_sire_genotyped"}, or \code{"known_sire_untyped"}),
#'     \code{n_offspring}, and genotyped flags.
#'   - crosses: named list keyed by \code{cross_id}, each holding the cross's
#'     \code{mother_id}, \code{father_id}, \code{family_type}, \code{offspring} IDs,
#'     maternal genotype \code{M}, paternal genotype \code{F} (\code{NA} when the sire
#'     is unknown or untyped), and offspring genotype matrix \code{G}.
#'   - parent_genotypes: named list keyed by parent ID storing each parent's genotype
#'     vector once (\code{NA} for a named-but-untyped sire); shared references so the
#'     same mother/sire across crosses carries one genotype.
#'   - F_list: paternal genotype vector per \code{cross_id} (\code{NA} vector when the
#'     sire is unknown/untyped).
#'
#'   The unknown sire is represented by the token \code{"__unknown_sire__"} in
#'   \code{father_id}; an open-pollinated cross keeps \code{cross_id == mother_id} so
#'   legacy \code{G_list}/\code{M_list} names are unchanged.
#'
#' @section Migration:
#'   Existing open-pollinated CSVs (mother + offspring, empty/NA \code{father}) read
#'   exactly as before: \code{G_list}/\code{M_list} keep their family_id names and
#'   values, so any code using them is unaffected. An old \code{HSMap.data} object
#'   created before this extension simply lacks the cross-aware fields; the full-sib
#'   API requires those fields, so re-read such data with the current
#'   \code{read_HSMap_data()} to populate them. Known father IDs in the pedigree are
#'   never dropped: they are recorded in \code{cross_table$father_id} and drive
#'   \code{family_type}.
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

  # ---- cross-aware structures (backward-compatible add-on) ------------------
  # The legacy family_id-keyed G_list/M_list above are left UNCHANGED so existing
  # open-pollinated workflows are byte-identical. The cross-aware structures below
  # group offspring by (mother, father) and retain sire identity/genotypes.
  cross <- .build_crosses(ped, X, markers)

  out <- list(
    G_list          = G_list,
    M_list          = M_list,
    alleles         = alleles,
    pedigree        = ped,
    stats           = stats_df,
    # cross-aware fields
    cross_table     = cross$cross_table,
    crosses         = cross$crosses,
    parent_genotypes = cross$parent_genotypes,
    F_list          = cross$F_list
  )
  class(out) <- "HSMap.data"
  out
}

# Token for an unknown / open-pollinated sire.
HSMAP_UNKNOWN_SIRE <- "__unknown_sire__"

# Build cross-aware structures from a pedigree and the marker x sample matrix X.
# Groups offspring by (mother_id, father_id); stores each parent's genotype ONCE by
# parent ID; assigns a family_type per cross; validates parent identity/coding.
.build_crosses <- function(ped, X, markers) {
  samp <- colnames(X)
  fa_raw <- ped$father
  # normalize an unknown father to the token
  fa <- ifelse(is.na(fa_raw) | !nzchar(as.character(fa_raw)), HSMAP_UNKNOWN_SIRE,
               as.character(fa_raw))
  mo <- as.character(ped$mother)

  # offspring rows = generation 2 individuals that are genotyped
  is_off <- ped$generation == 2L & ped$id %in% samp
  off_ped <- ped[is_off, , drop = FALSE]
  off_fa  <- fa[is_off]
  off_mo  <- mo[is_off]

  # consistency: an offspring id must not appear with conflicting parents
  dup_off <- off_ped$id[duplicated(off_ped$id)]
  if (length(dup_off)) {
    bad <- vapply(unique(dup_off), function(idi) {
      r <- off_ped[off_ped$id == idi, , drop = FALSE]
      length(unique(paste(r$mother, r$father))) > 1L
    }, logical(1))
    if (any(bad))
      stop("read_HSMap_data(): offspring with conflicting parents in pedigree: ",
           paste(unique(dup_off)[bad], collapse = ", "), call. = FALSE)
  }

  # cross key: mother_id for an unknown sire (backward-compatible with the legacy
  # family_id == mother naming); mother__x__father for a known sire.
  ckey <- ifelse(off_fa == HSMAP_UNKNOWN_SIRE, off_mo,
                 paste(off_mo, off_fa, sep = "__x__"))
  cross_ids <- unique(ckey)

  # parent genotypes stored once by parent id (mothers and fathers)
  all_parents <- unique(c(off_mo, off_fa[off_fa != HSMAP_UNKNOWN_SIRE]))
  parent_genotypes <- stats::setNames(vector("list", length(all_parents)), all_parents)
  for (p in all_parents) {
    if (p %in% samp) {
      g <- X[, p]; names(g) <- markers
      parent_genotypes[[p]] <- as.integer(g)
    } else {
      parent_genotypes[[p]] <- NA   # known id but not genotyped
    }
  }

  crosses <- stats::setNames(vector("list", length(cross_ids)), cross_ids)
  F_list  <- stats::setNames(vector("list", length(cross_ids)), cross_ids)
  ct_rows <- vector("list", length(cross_ids))
  for (i in seq_along(cross_ids)) {
    cid <- cross_ids[i]
    sel <- ckey == cid
    mom_id <- unique(off_mo[sel]); mom_id <- mom_id[1]
    fat_id <- unique(off_fa[sel]); fat_id <- fat_id[1]
    kids   <- off_ped$id[sel]

    Mvec <- rep(NA_integer_, length(markers)); names(Mvec) <- markers
    if (mom_id %in% samp) Mvec <- as.integer(X[, mom_id])
    names(Mvec) <- markers

    Fvec <- rep(NA_integer_, length(markers)); names(Fvec) <- markers
    father_known <- fat_id != HSMAP_UNKNOWN_SIRE
    father_typed <- father_known && (fat_id %in% samp)
    if (father_typed) { Fvec <- as.integer(X[, fat_id]); names(Fvec) <- markers }

    ftype <- if (!father_known) "open_pollinated"
             else if (father_typed) "known_sire_genotyped"
             else "known_sire_untyped"

    Gmat <- if (length(kids)) t(X[, kids, drop = FALSE]) else
      matrix(NA_integer_, 0, length(markers), dimnames = list(character(0), markers))
    storage.mode(Gmat) <- "integer"; colnames(Gmat) <- markers

    crosses[[cid]] <- list(
      cross_id = cid, mother_id = mom_id, father_id = fat_id,
      family_type = ftype, offspring = kids,
      M = Mvec, F = Fvec, G = Gmat
    )
    F_list[[cid]] <- Fvec
    ct_rows[[i]] <- data.frame(
      cross_id = cid, mother_id = mom_id, father_id = fat_id,
      family_type = ftype, n_offspring = length(kids),
      mother_genotyped = mom_id %in% samp, father_genotyped = father_typed,
      stringsAsFactors = FALSE)
  }
  cross_table <- do.call(rbind, ct_rows)

  # validation: a genotyped parent must actually be in the genotype file (it is, by
  # construction); a mother required for a usable cross must be genotyped.
  no_mom <- cross_table$cross_id[!cross_table$mother_genotyped]
  if (length(no_mom))
    warning("read_HSMap_data(): mother genotype missing for cross(es): ",
            paste(utils::head(no_mom, 5), collapse = ", "),
            if (length(no_mom) > 5) " ..." else "", call. = FALSE)
  untyped <- cross_table$cross_id[cross_table$family_type == "known_sire_untyped"]
  if (length(untyped))
    warning("read_HSMap_data(): ", length(untyped), " cross(es) name a sire that is ",
            "NOT genotyped (family_type = 'known_sire_untyped'); these are not treated ",
            "as full-sib unless you opt into the open-pollinated fallback. Example: ",
            paste(utils::head(untyped, 3), collapse = ", "), call. = FALSE)

  list(cross_table = cross_table, crosses = crosses,
       parent_genotypes = parent_genotypes, F_list = F_list)
}
