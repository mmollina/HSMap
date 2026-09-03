#' Simulate offspring genotypes for one population
#'
#' @description
#' Low-level wrapper around the C++ core `simulate_offspring_cpp()` that
#' simulates all offspring in a *single* population (nuclear family with a
#' known mother and unknown or mixed paternal contribution).
#'
#' @details
#' - Maternal genotypes are supplied as an integer vector `M` with entries
#'   `0,1,2,NA` for `aa, Aa, AA, missing`. Missing maternal markers produce
#'   `NA` in the corresponding offspring genotypes.
#' - The paternal contribution is modeled per marker by a 3 x T matrix
#'   `pi_true` (rows: `AA, Aa, aa`) whose **columns** represent the paternal
#'   genotype mixture for each marker. Non-finite values are disallowed.
#'   Columns are normalized here if they do not sum to 1.
#' - Recombination follows a 2-state hidden path `H` with
#'   `Pr(H_{t+1} != H_t) = r_true[t]`. The optional `z_phase` (0/1 per marker)
#'   flips the maternal allele used at heterozygous maternal loci, which lets
#'   you impose coupling or repulsion along the chromosome.
#' - Genotyping error (if `error_rate > 0`) randomly replaces a generated
#'   call `Y_t in {0,1,2}` with a *different* value, uniformly over the other
#'   two values.
#'
#' @param M Integer vector (length T) with entries in `0,1,2,NA`,
#'   the maternal genotype per marker (0=aa, 1=Aa, 2=AA).
#' @param pi_true Numeric matrix `3 x T` with rows `AA, Aa, aa`;
#'   paternal mixture per marker. Columns must be positive and are normalized
#'   here if needed.
#' @param r_true Numeric vector length `T-1`; recombination fractions between
#'   adjacent markers.
#' @param n_offspring Integer, number of offspring to simulate.
#' @param error_rate Numeric in `[0,1]`; per-marker probability to replace a
#'   simulated `Y_t` with a different value in `{0,1,2}`.
#' @param keep_paths Logical; if `TRUE`, also return the hidden 0/1 paths `H`.
#' @param z_phase Optional integer vector length `T` with values in `{0,1}`.
#'   When the mother is heterozygous at marker `t`, the transmitted allele is
#'   `(H_t XOR z_phase[t])`. Use this to encode a desired coupling/repulsion
#'   pattern when simulating.
#' @param miss_rate Numeric in [0,1]; additional random missingness rate applied
#'   to offspring genotype calls *after* simulation (MCAR). Default 0 means no
#'   extra missingness.
#'
#' @return
#' - If `keep_paths = FALSE`: list with `G` (integer matrix `n_offspring x T`
#'   with entries `0/1/2/NA`).
#' - If `keep_paths = TRUE`: list with `G` and `H` (hidden 0/1 paths).
#'
#' @seealso [sim_multi_pop()] for a higher-level multi-population
#'   generator; `simulate_offspring_cpp()` for the C++ core.
#' @keywords internal
#' @noRd
sim_family_genotypes <- function(M,
                                 pi_true,
                                 r_true,
                                 n_offspring,
                                 error_rate = 0,
                                 keep_paths = FALSE,
                                 z_phase = NULL,
                                 miss_rate = 0) {
  # Basic coercions
  M        <- as.integer(M)
  pi_true  <- as.matrix(pi_true)
  r_true   <- as.numeric(r_true)

  # Shape checks
  stopifnot(
    nrow(pi_true) == 3L,
    length(M)     == ncol(pi_true),
    length(r_true) == length(M) - 1L
  )

  # Optional z_phase
  if (!is.null(z_phase)) {
    z_phase <- as.integer(z_phase)
    stopifnot(length(z_phase) == length(M))
  }

  if (!is.numeric(miss_rate) || length(miss_rate) != 1L ||
      !is.finite(miss_rate) || miss_rate < 0 || miss_rate > 1) {
    stop("`miss_rate` must be a single number in [0,1].")
  }

  # Validate paternal mixture at present markers
  ok <- which(!is.na(M))
  if (length(ok)) {
    if (any(!is.finite(pi_true[, ok, drop = FALSE])))
      stop("`pi_true` has NA/Inf at present markers")
    s <- colSums(pi_true[, ok, drop = FALSE])
    if (any(!is.finite(s) | s <= 0))
      stop("`pi_true` columns at present markers must sum to positive values")
    if (any(abs(s - 1) > 1e-6)) {
      # Normalize columns that differ numerically from 1
      pi_true[, ok] <- sweep(pi_true[, ok, drop = FALSE], 2, s, "/")
    }
  }

  # Delegate to C++ core (see src/cpp)
  out <- simulate_offspring_cpp(
    M               = M,
    pi_true         = pi_true,
    r_true          = r_true,
    n_offspring     = as.integer(n_offspring),
    error_rate      = as.numeric(error_rate),
    keep_paths      = isTRUE(keep_paths),
    z_phase_in      = if (is.null(z_phase)) NULL else z_phase
  )

  # Optional MCAR masking on G
  if (miss_rate > 0) {
    G <- out$G
    mask <- matrix(stats::runif(length(G)) < miss_rate, nrow = nrow(G))
    # only overwrite non-NA entries
    G[mask & !is.na(G)] <- NA_integer_
    out$G <- G
  }

  out
}
