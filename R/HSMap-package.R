#' HSMap: maternal linkage mapping in open-pollinated diploid families
#'
#' HSMap builds maternal linkage maps from open-pollinated / unknown-sire diploid
#' half-sib families: one known, genotyped dam, many offspring, unobserved fathers.
#' The production chromosome-wide estimator, \code{\link{hmm_map_source_aware}}, is a
#' hidden Markov model whose state pairs the transmitted maternal homolog with a
#' paternal source state, so persistent paternal haplotype-sharing tracts (as arise
#' when some pollen parents are related to the dam) are modelled explicitly rather
#' than absorbed into the maternal map; setting the sharing entry rate to zero
#' recovers the population-mode estimator \code{\link{hmm_map}} exactly. Several dams
#' can be combined by a joint EM that estimates one shared recombination map while
#' keeping phase and the paternal model dam-specific (companion-paper extension).
#' Maps are fitted within resolved phase blocks and no-linkage intervals are reported
#' as gaps. See the "Getting started" vignette and \code{\link{hmm_map_blocks}}.
#'
#' @keywords internal
#' @import RcppParallel
#' @importFrom Rcpp evalCpp sourceCpp
#' @useDynLib HSMap, .registration = TRUE
"_PACKAGE"
