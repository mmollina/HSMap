#' HSMap: maternal linkage mapping in open-pollinated diploid families
#'
#' HSMap builds maternal linkage maps from open-pollinated / unknown-sire diploid
#' half-sib families using a fast C++ (Rcpp) maternal hidden Markov model. The unknown
#' paternal contribution is integrated out through a per-marker, dam-specific paternal
#' gametic frequency; several dams are combined by a joint EM that estimates one shared
#' recombination map while keeping phase and the paternal model dam-specific. Maps are
#' fitted within resolved phase blocks and no-linkage intervals are reported as gaps.
#' See the "Getting started" vignette and \code{\link{hmm_map_blocks}}.
#'
#' @keywords internal
#' @import RcppParallel
#' @importFrom Rcpp evalCpp sourceCpp
#' @useDynLib HSMap, .registration = TRUE
"_PACKAGE"
