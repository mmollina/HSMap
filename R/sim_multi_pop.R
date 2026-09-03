#' Simulate multi-population (or pedigree-derived) datasets for the maternal-only HMM
#'
#' @description
#' Generate synthetic genotype data for one or more populations, each treated
#' as a nuclear family with a known mother and multiple offspring. Paternal
#' contribution can be unknown (mixture per marker) or known (single sire).
#' Optionally, input a pedigree to form families; otherwise families are
#' created as independent populations. This function orchestrates per-population
#' simulation by calling an internal single-population core.
#'
#' @details
#' **Markers and recombination**
#'
#' - The union of markers has length `T_markers = T`. If a pedigree is not
#'   provided, each population is assigned a subset of these markers according
#'   to `marker_intersection`. Shared markers retain the same order across
#'   populations; non-shared markers are scattered (not blocky).
#' - The recombination vector `r_true` is common to all populations, either
#'   supplied via `r_vec` (length `T-1`) or filled as a constant `r_const`.
#'
#' **Maternal genotypes**
#'
#' Controlled by `maternal_geno_mode`:
#' - `"HWE"` (default): for markers present in a population, simulate maternal
#'   genotype under HWE at allele-A frequency `maternal_pA` (scalar or vector
#'   of length `T`). Absent markers remain `NA`.
#' - `"all_het"`: force all present maternal markers to be heterozygous (Aa).
#' - `"vector"`: use the provided integer vector `maternal_M_given` (0/1/2/NA),
#'   aligned to the union marker set. Absent markers in a population are kept
#'   as `NA` regardless of values in `maternal_M_given`.
#'
#' **Paternal composition per marker**
#'
#' - If a population has an *unknown sire*, the paternal genotype mixture at
#'   present markers follows HWE with allele-A frequency
#'   `paternal_pA_base` (scalar or length-T vector). To induce population
#'   heterogeneity, a zero-mean Gaussian jitter with SD `paternal_pA_sd` is
#'   added **per marker** and truncated into `(1e-4, 1-1e-4)`.
#' - If a population has a *known sire* (probability
#'   `known_sire_prob_per_pop` when no pedigree is given), the sire genotype is
#'   drawn **per marker** under HWE at the population's perturbed allele frequency
#'   and encoded as a one-hot column in the 3 x T `pi_true` matrix for that
#'   population. This is a **simulation-only** device that fixes a paternal
#'   *genotype* composition; the paternal genotypes are generated **independently
#'   across markers** (there is no linked paternal meiosis and no paternal
#'   haplotype). It does **not** correspond to, and must not be read as, support for
#'   known-sire / full-sib linkage mapping: the estimation functions in this package
#'   model the open-pollinated / unknown-sire case only (the paternal contribution is
#'   integrated out per marker through a dam-specific gametic frequency).
#'
#' **Phase along the chromosome**
#'
#' - Maternal phase per population is drawn by `generate_phase()` using
#'   `phase_mode`:
#'   - `"all_coupling"`: all intervals coupling;
#'   - `"random"`: each interval is coupling with probability `1 - repulsion_rate`
#'     and repulsion with probability `repulsion_rate`;
#'   - `"vector"`: supplied `phase_vector` of length `T-1` with entries `0/1`.
#' - Internally we maintain both a `{+1, -1}` representation (`z`) and the
#'   standard `{0,1}` interval representation (`v`). We pass `{0,1}` to the
#'   C++ simulator by converting `z` to a per-marker flip vector.
#'
#' **Genotyping error**
#'
#' - Each offspring/marker call is independently replaced with a *different*
#'   genotype with probability `error_rate`. This matches the C++ core behavior.
#'
#' **Reproducibility**
#'
#' - If `seed` is provided, `set.seed(seed)` is called at the start. Note that
#'   downstream C++ draws share R's RNG state via `RNGScope`, so the simulation
#'   is fully reproducible under the same R and package versions.
#'
#' @param T_markers Integer. Total number of markers in the union (chromosome length `T`).
#' @param n_pops Integer. Number of populations (nuclear families) to simulate.
#' @param n_ind_per_pop Integer vector of length `n_pops`. Number of offspring per population.
#' @param marker_intersection Numeric in `[0,1]`. Fraction of markers shared by all populations.
#'   Non-shared markers are scattered across the non-intersection while preserving order.
#' @param r_vec Optional numeric of length `T_markers - 1`. If supplied, used as
#'   the recombination vector for all populations.
#' @param r_const Numeric. Constant recombination fraction per interval when `r_vec` is `NULL`.
#' @param phase_mode Character, one of `"all_coupling"`, `"random"`, `"vector"`.
#' @param repulsion_rate Numeric in `[0,1]`. Used only when `phase_mode = "random"`, the
#'   fraction of intervals simulated in repulsion.
#' @param phase_vector Integer vector of length `T_markers - 1` with entries `0/1`,
#'   used only when `phase_mode = "vector"`. Ignored otherwise.
#' @param maternal_geno_mode Character, one of `"HWE"`, `"all_het"`, `"vector"`.
#' @param maternal_M_given Optional integer vector (0/1/2/NA) of length `T_markers`
#'   or named by union markers; used only when `maternal_geno_mode = "vector"`.
#' @param maternal_pA Numeric scalar or length-`T_markers` vector in `(0,1)` for HWE
#'   maternal genotypes when `maternal_geno_mode = "HWE"`. Values are clipped to `(1e-6, 1-1e-6)`.
#' @param paternal_pA_base Numeric scalar or length-`T_markers` vector in `(0,1)` for HWE
#'   paternal mixture per marker. Values are clipped to `(1e-6, 1-1e-6)` before population jitter.
#' @param paternal_pA_sd Non-negative numeric. Per-population Gaussian perturbation SD
#'   added to `paternal_pA_base` *per marker*. The perturbed values are then
#'   truncated into `(1e-4, 1-1e-4)`.
#' @param known_sire_prob_per_pop Numeric in `[0,1]`. When no pedigree is passed, each
#'   population independently has a known sire with this probability. \strong{This is a
#'   simulation-only, fixed paternal-genotype mechanism:} the sire's genotype is drawn
#'   independently per marker (no linked paternal meiosis, no paternal haplotype), and it
#'   does \strong{not} imply or enable known-sire / full-sib mapping --- the estimators
#'   here model the open-pollinated / unknown-sire case only.
#' @param error_rate Numeric in `[0,1]`. Per offspring/marker error rate for replacing
#'   a generated genotype with a different value.
#' @param pedigree Optional `data.frame` with columns `id, mother, father, generation, family_id`.
#'   If provided, families are derived from the pedigree; all families are assigned the full
#'   union marker set in order.
#' @param seed Optional integer for reproducibility.
#' @param keep_paths Logical; if `TRUE`, the hidden paths `H` per population are returned.
#' @param miss_rate Numeric scalar in `[0,1]` or length-`n_pops` vector. Additional
#'   random missingness (MCAR) applied per population to offspring genotypes.
#'   Default 0 (no extra missingness).
#'
#' @return A list of class `"sim_multi_pop"` with elements:
#' \describe{
#'   \item{G_list}{List of length `n_pops`. Each element is an integer matrix
#'     `n_offspring(g) x T` with entries `0/1/2/NA`, the simulated offspring genotypes.}
#'   \item{M_list}{List of length `n_pops`. Each element is an integer vector
#'     length `T` with entries `0/1/2/NA`, the maternal genotype per marker
#'     (missing where the marker is absent in that population if applicable).}
#'   \item{pi_true_list}{List of length `n_pops`. Each element is a numeric matrix
#'     `3 x T` (rows `AA,Aa,aa`) with the paternal mixture actually used to
#'     generate the data in that population. Columns sum to 1 at present markers.}
#'   \item{pi_prior_list}{Same shape as `pi_true_list`. Default equals the truth here,
#'     provided for downstream workflows that expect priors.}
#'   \item{pi_fixed_list}{List of `3 x T` matrices with one-hot columns at markers where
#'     the sire is known, `NA` elsewhere.}
#'   \item{father_geno_list}{List of integer vectors length `T` with entries `0/1/2/NA`,
#'     the drawn sire genotype per marker when a known sire is present; `NA` otherwise.}
#'   \item{truth}{A list with ground truth for reproducibility:
#'     \describe{
#'       \item{markers_union}{Character vector of union marker names, length `T`.}
#'       \item{r_true}{Numeric vector length `T-1`, recombination fractions.}
#'       \item{z_true}{List of length `n_pops`, each a `{+1,-1}` per-marker phase vector.}
#'       \item{v_true}{List of length `n_pops`, each a `{0,1}` per-interval phase vector.}
#'       \item{maternal_pA}{Named numeric vector length `T`, maternal A-allele frequency used for HWE.}
#'       \item{paternal_pA_base}{Named numeric vector length `T`, the base paternal A-allele frequency before per-population jitter.}
#'       \item{pops_meta}{List with per-population metadata: IDs, assigned marker sets, and offspring IDs.}
#'     }}
#'   \item{H_paths_list}{If `keep_paths = TRUE`, a list of hidden 0/1 path matrices
#'     `n_offspring(g) x T`, otherwise `NULL`.}
#' }
#'
#' @section Column and row naming:
#' Offspring matrices in `G_list` have row names equal to `offspring_ids` and
#' column names equal to the union marker names. The same column naming applies
#' to `H_paths_list` and all per-marker vectors/matrices in the return.
#'
#' @section Invariants and checks:
#' - `nrow(pi_true) == 3`, `ncol(pi_true) == T`, columns at present markers sum to 1
#'   (normalized internally if small numerical deviations occur).
#' - `length(r_true) == T-1`.
#' - Maternal genotype vectors are integer in `0/1/2/NA`. When `maternal_geno_mode = "vector"`,
#'   values outside this set trigger an error.
#'
#' @examples
#' \dontrun{
#' # 1) Single population, all markers shared, all moms Aa, unknown sire HWE at p=0.4:
#' sim1 <- sim_multi_pop(
#'   T_markers = 50, n_pops = 1, n_ind_per_pop = 200,
#'   marker_intersection = 1,
#'   r_const = 0.02,
#'   phase_mode = "all_coupling",
#'   maternal_geno_mode = "all_het",
#'   paternal_pA_base = 0.4,
#'   error_rate = 0.01,
#'   seed = 11
#' )
#' str(sim1$G_list[[1]])
#'
#' # 2) Two populations with different paternal allele frequencies per marker:
#' T <- 30
#' pA_base <- runif(T, 0.2, 0.8)   # per-marker paternal pA
#' sim2 <- sim_multi_pop(
#'   T_markers = T, n_pops = 2, n_ind_per_pop = c(80, 60),
#'   marker_intersection = 1,
#'   r_const = 0.03,
#'   phase_mode = "random", repulsion_rate = 0.25,
#'   maternal_geno_mode = "HWE", maternal_pA = 0.55,
#'   paternal_pA_base = pA_base, paternal_pA_sd = 0.05,
#'   error_rate = 0.01, seed = 7
#' )
#' # Check per-marker paternal mixture used in pop 1:
#' head(t(sim2$pi_true_list[[1]]))
#'
#' # 3) Vector maternal genotypes (fixed) and known sires in some populations:
#' T <- 20
#' M_vec <- sample(c(0L,1L,2L), T, replace = TRUE)  # fixed mom genotypes
#' sim3 <- sim_multi_pop(
#'   T_markers = T, n_pops = 3, n_ind_per_pop = c(50, 50, 50),
#'   marker_intersection = 0.8,
#'   r_const = 0.02,
#'   phase_mode = "vector", phase_vector = sample(0:1, T-1, TRUE),
#'   maternal_geno_mode = "vector", maternal_M_given = M_vec,
#'   paternal_pA_base = 0.5, known_sire_prob_per_pop = 0.3,
#'   error_rate = 0.005, seed = 123
#' )
#' }
#'
#' @export
sim_multi_pop <- function(
    T_markers,
    n_pops = 2,
    n_ind_per_pop = rep(50, n_pops),
    marker_intersection = 1,
    r_vec = NULL,
    r_const = 0.05,
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
    miss_rate = 0
){
  if (!is.null(seed)) set.seed(seed)
  Tm <- as.integer(T_markers)
  if (Tm < 2) stop("Need at least 2 markers")

  maternal_geno_mode <- match.arg(maternal_geno_mode)
  phase_mode         <- match.arg(phase_mode)

  # --- Map ---
  r_true <- if (is.null(r_vec)) {
    rep(r_const, Tm - 1)
  } else {
    if (length(r_vec) != (Tm - 1)) stop("`r_vec` must have length T_markers - 1")
    as.numeric(r_vec)
  }

  # normalize miss_rate to length n_pops
  if (length(miss_rate) == 1L) miss_rate <- rep(miss_rate, n_pops)
  if (length(miss_rate) != n_pops || any(!is.finite(miss_rate)) ||
      any(miss_rate < 0 | miss_rate > 1)) {
    stop("`miss_rate` must be a scalar or length-n_pops vector with values in [0,1].")
  }

  # --- Phase per population ---
  z_true_list <- vector("list", n_pops)
  v_true_list <- vector("list", n_pops)
  for (g in seq_len(n_pops)) {
    ph_g <- generate_phase(
      Tm,
      mode         = phase_mode,
      repulsion_rate = repulsion_rate,
      v_vector     = if (!is.null(phase_vector) && length(phase_vector) == Tm - 1) phase_vector else NULL
    )
    z_true_list[[g]] <- ph_g$z
    v_true_list[[g]] <- ph_g$v
  }

  # --- Marker sets (define union_markers) ---
  if (is.null(pedigree)) {
    if (length(n_ind_per_pop) != n_pops) stop("`n_ind_per_pop` must have length `n_pops`")
    marker_assign   <- assign_marker_sets(Tm, n_pops, marker_intersection)
    union_markers   <- marker_assign$markers_union
    per_pop_markers <- marker_assign$per_pop
    pops <- lapply(seq_len(n_pops), function(g) {
      list(
        pop_id       = paste0("P", g),
        mother_id    = paste0("MOM", g),
        # Simulation-only: with prob. known_sire_prob_per_pop this population gets a
        # fixed sire GENOTYPE (drawn independently per marker, no linked paternal
        # meiosis). This does NOT enable known-sire/full-sib mapping; the estimators
        # model the open-pollinated / unknown-sire case only.
        father_id    = ifelse(stats::runif(1) < known_sire_prob_per_pop, paste0("DAD", g), NA_character_),
        offspring_ids= paste0("P", g, "_O", seq_len(n_ind_per_pop[g])),
        markers      = per_pop_markers[[g]]
      )
    })
  } else {
    fams <- families_from_pedigree(pedigree)
    n_pops <- length(fams)
    pops <- fams
    union_markers <- paste0("m", seq_len(Tm))
    for (k in seq_along(pops)) pops[[k]]$markers <- union_markers
  }

  # --- Allele frequencies (AFTER union_markers known) ---
  to_lenT <- function(x) if (length(x) == 1L) rep(x, Tm) else x
  pA_m    <- stats::setNames(pmin(pmax(to_lenT(maternal_pA),      1e-6), 1 - 1e-6), union_markers)
  pA_base <- stats::setNames(pmin(pmax(to_lenT(paternal_pA_base), 1e-6), 1 - 1e-6), union_markers)

  # --- If user provides a maternal genotype vector, align it ---
  if (identical(maternal_geno_mode, "vector")) {
    if (is.null(maternal_M_given)) stop("`maternal_M_given` must be provided when maternal_geno_mode='vector'")
    if (!is.null(names(maternal_M_given))) {
      maternal_M_given <- as.integer(maternal_M_given[union_markers])
    } else {
      if (length(maternal_M_given) != Tm) stop("`maternal_M_given` must have length T_markers")
      maternal_M_given <- as.integer(maternal_M_given)
    }
    if (!all(maternal_M_given %in% c(0L,1L,2L,NA_integer_)))
      stop("`maternal_M_given` must contain only 0/1/2/NA")
  }

  # --- Outputs ---
  G_list <- vector("list", n_pops)
  M_list <- vector("list", n_pops)
  pi_true_list     <- vector("list", n_pops)
  pi_prior_list    <- vector("list", n_pops)
  pi_fixed_list    <- vector("list", n_pops)
  father_geno_list <- vector("list", n_pops)
  H_paths_list     <- if (keep_paths) vector("list", n_pops) else NULL
  names(G_list) <- names(M_list) <- names(pi_true_list) <- names(pi_prior_list) <-
    names(pi_fixed_list) <- names(father_geno_list) <- vapply(pops, `[[`, character(1L), "pop_id")

  # --- Population loop ---
  for (g in seq_len(n_pops)) {
    pop  <- pops[[g]]
    cols <- pop$markers  # marker names present in this population

    # Per-pop paternal AF jitter (if any), then truncate
    if (paternal_pA_sd > 0) {
      delta  <- stats::setNames(stats::rnorm(Tm, 0, paternal_pA_sd), union_markers)
      pA_pat <- pmin(pmax(pA_base + delta, 1e-4), 1 - 1e-4)
    } else {
      pA_pat <- pA_base
    }

    # Maternal genotypes on the union (NA where absent)
    M_g <- rep(NA_integer_, Tm); names(M_g) <- union_markers
    if (identical(maternal_geno_mode, "HWE")) {
      M_g[cols] <- simulate_maternal_geno(p_m = pA_m[cols])
    } else if (identical(maternal_geno_mode, "all_het")) {
      M_g[cols] <- 1L
    } else { # "vector"
      M_g <- maternal_M_given
      # ensure absent markers in this population stay NA
      M_g[setdiff(union_markers, cols)] <- NA_integer_
    }

    # Paternal composition (true) and optional known sire
    pi_true     <- matrix(NA_real_, 3, Tm, dimnames = list(c("AA","Aa","aa"), union_markers))
    father_geno <- rep(NA_integer_, Tm); names(father_geno) <- union_markers

    if (is.na(pop$father_id)) {
      # Unknown sire: HWE mixture at pA_pat (rows are "AA","Aa","aa")
      p  <- pA_pat[cols]
      AA <- p^2
      Aa <- 2 * p * (1 - p)
      aa <- (1 - p)^2
      pi_true[c("AA","Aa","aa"), cols] <- rbind(AA, Aa, aa)
      cs <- colSums(pi_true[, cols, drop = FALSE])
      pi_true[, cols] <- sweep(pi_true[, cols, drop = FALSE], 2, cs, "/")
    } else {
      # Known sire: draw genotype per marker; one-hot pi
      father_geno[cols] <- rHWE_geno012(pA_pat[cols])
      pi_true[, cols]   <- geno012_to_pi(father_geno[cols])
    }

    # Priors (default to truth) and fixed-pi (one-hot where sire known)
    pi_prior <- pi_true
    pi_fixed <- matrix(NA_real_, 3, Tm, dimnames = list(c("AA","Aa","aa"), union_markers))
    if (any(!is.na(father_geno))) {
      pi_fixed[, !is.na(father_geno)] <- geno012_to_pi(father_geno[!is.na(father_geno)])
    }

    # Convert phase from {+1, -1} to {0, 1} per marker for the C++ call.
    # Convention: Coupling (z=+1) -> 0; Repulsion (z=-1) -> 1
    n_ind <- length(pop$offspring_ids)
    z_phase_01 <- as.integer(z_true_list[[g]] == -1L)

    res_cpp <- sim_family_genotypes(
      M           = M_g,
      pi_true     = pi_true,
      r_true      = r_true,
      n_offspring = n_ind,
      error_rate  = error_rate,
      keep_paths  = keep_paths,
      z_phase     = z_phase_01,
      miss_rate   = miss_rate[g]
    )

    Gg <- res_cpp$G
    rownames(Gg) <- pop$offspring_ids
    colnames(Gg) <- union_markers
    if (keep_paths) {
      H_paths_list[[g]] <- res_cpp$H
      rownames(H_paths_list[[g]]) <- pop$offspring_ids
      colnames(H_paths_list[[g]]) <- union_markers
    }

    # Write outputs
    G_list[[g]]           <- Gg
    M_list[[g]]           <- M_g
    pi_true_list[[g]]     <- pi_true
    pi_prior_list[[g]]    <- pi_prior
    pi_fixed_list[[g]]    <- pi_fixed
    father_geno_list[[g]] <- father_geno
  }

  # Truth bundle
  truth <- list(
    markers_union    = union_markers,
    r_true           = r_true,
    z_true           = z_true_list,
    v_true           = v_true_list,
    maternal_pA      = pA_m,
    paternal_pA_base = pA_base,
    pops_meta        = pops
  )

  out <- list(
    G_list           = G_list,
    M_list           = M_list,
    pi_true_list     = pi_true_list,
    pi_prior_list    = pi_prior_list,
    pi_fixed_list    = pi_fixed_list,
    father_geno_list = father_geno_list,
    truth            = truth,
    H_paths_list     = H_paths_list
  )
  class(out) <- "sim_multi_pop"
  out
}
