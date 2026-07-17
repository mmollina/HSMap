# Shared utilities for the known-sire / full-sib extension.
# See dev/known_sire_design.md (Section 3) for the phase representation.

#' Build phased haplotype allele labels from a genotype and an adjacent phase vector
#'
#' @description
#' Construct a parent's two homolog allele sequences (the oracle "phased haplotypes")
#' from its genotype and an adjacent phase vector, using the package's standard
#' anchored convention: homolog labelling is carried forward unchanged across a
#' coupling interval and swapped across a repulsion interval, anchored so that the
#' first heterozygous marker places the A allele on homolog 1. Phase is thus applied
#' exactly once, in this construction; the resulting labels feed the full-sib emission.
#'
#' @param genotype Integer vector (length z) of parent genotypes coded 0/1/2 = aa/Aa/AA
#'   (\code{NA} allowed for a missing genotype at a marker).
#' @param phase_vec Integer vector (length z-1) of adjacent phases, 1 = coupling,
#'   0 = repulsion. \code{NA} entries are treated as coupling for the orientation
#'   chain (they carry the labelling forward unchanged).
#'
#' @return An integer matrix \code{2 x z} with rows = homolog 1 and 2 and entries in
#'   \code{{0,1}} (allele a / A); \code{NA} at markers with a missing genotype.
#' @export
phase_to_haplotypes <- function(genotype, phase_vec) {
  g <- as.integer(genotype)
  z <- length(g)
  if (z < 1L) stop("`genotype` must have length >= 1.")
  if (length(phase_vec) != z - 1L)
    stop("`phase_vec` must have length length(genotype) - 1.")
  pv <- as.integer(phase_vec)
  # orientation parity: flip at every repulsion interval (NA treated as coupling)
  orient <- integer(z)
  for (k in 2:z) {
    flip <- if (z >= 2 && !is.na(pv[k - 1L]) && pv[k - 1L] == 0L) 1L else 0L
    orient[k] <- bitwXor(orient[k - 1L], flip)
  }
  if (z == 1L) orient <- 0L
  H <- matrix(NA_integer_, 2L, z)
  for (k in seq_len(z)) {
    gk <- g[k]
    if (is.na(gk)) next
    if (gk == 2L)      { H[1L, k] <- 1L; H[2L, k] <- 1L }
    else if (gk == 0L) { H[1L, k] <- 0L; H[2L, k] <- 0L }
    else {                                   # heterozygous
      a1 <- if (orient[k] == 0L) 1L else 0L
      H[1L, k] <- a1; H[2L, k] <- 1L - a1
    }
  }
  H
}

#' Recover an adjacent phase vector from phased haplotypes (inverse of
#' \code{phase_to_haplotypes} at heterozygous markers)
#'
#' @param H Integer matrix \code{2 x z} of homolog alleles (0/1).
#' @return Integer vector length z-1: 1 = coupling, 0 = repulsion; \code{NA} for an
#'   interval whose two markers are not both heterozygous (phase undefined there).
#' @export
haplotypes_to_phase <- function(H) {
  z <- ncol(H)
  het <- H[1L, ] != H[2L, ]
  pv <- rep(NA_integer_, max(z - 1L, 0L))
  for (k in seq_len(z - 1L)) {
    if (isTRUE(het[k]) && isTRUE(het[k + 1L]) &&
        !anyNA(H[, k]) && !anyNA(H[, k + 1L])) {
      pv[k] <- as.integer(H[1L, k] == H[1L, k + 1L])  # same allele on homolog 1 -> coupling
    }
  }
  pv
}
