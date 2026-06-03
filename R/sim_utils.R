#' Generate maternal phase along markers
#'
#' @description
#' Produces a phase sign vector \eqn{z \in \{+1,-1\}^T} and a flip indicator
#' vector \eqn{v \in \{0,1\}^{T-1}} with \eqn{v_t = 1\{z_{t+1} \neq z_t\}}.
#' Modes:
#' \itemize{
#'   \item \code{"all_coupling"}: \eqn{z_t \equiv +1}, \eqn{v \equiv 0}.
#'   \item \code{"random"}: independent flips with probability \code{repulsion_rate}.
#'   \item \code{"vector"}: use a user-supplied 0/1 vector \code{v_vector}.
#' }
#'
#' @param T integer (>= 2). Number of markers.
#' @param mode character; one of \code{"all_coupling"}, \code{"random"}, \code{"vector"}.
#' @param repulsion_rate numeric in \[0,1\]. Flip probability per interval for
#'   \code{mode = "random"}.
#' @param v_vector integer vector of length \code{T-1} with entries \code{0/1}
#'   indicating flips; used when \code{mode = "vector"}.
#'
#' @return A list with \code{z} (length \code{T}, values \code{+1/-1}) and
#'   \code{v} (length \code{T-1}, values \code{0/1}).
#'
#' @examples
#' \dontrun{
#' generate_phase(10, mode = "all_coupling")
#' generate_phase(10, mode = "random", repulsion_rate = 0.2)
#' }
#' @keywords internal
#' @importFrom stats rbinom
#' @noRd
generate_phase <- function(T,
                           mode = c("all_coupling", "random", "vector"),
                           repulsion_rate = 0.2,
                           v_vector = NULL) {
  if (!is.numeric(T) || length(T) != 1L || T < 2 || !is.finite(T)) {
    stop("`T` must be a single finite integer >= 2.")
  }
  T <- as.integer(T)
  mode <- match.arg(mode)

  if (mode == "vector") {
    if (is.null(v_vector) || length(v_vector) != (T - 1L)) {
      stop("`v_vector` must be provided and have length T-1 when mode = 'vector'.")
    }
    v <- as.integer(v_vector != 0L)
    z <- integer(T); z[1] <- 1L
    for (t in seq_len(T - 1L)) z[t + 1L] <- if (v[t] == 1L) -z[t] else z[t]
    return(list(z = z, v = v))
  }

  if (mode == "all_coupling") {
    return(list(z = rep.int(1L, T), v = rep.int(0L, T - 1L)))
  }

  # random flips
  if (!is.numeric(repulsion_rate) || length(repulsion_rate) != 1L ||
      !is.finite(repulsion_rate) || repulsion_rate < 0 || repulsion_rate > 1) {
    stop("`repulsion_rate` must be a single number in [0, 1].")
  }
  v <- rbinom(T - 1L, size = 1L, prob = repulsion_rate)
  z <- integer(T); z[1] <- 1L
  for (t in seq_len(T - 1L)) z[t + 1L] <- if (v[t] == 1L) -z[t] else z[t]
  list(z = z, v = v)
}

#' Assign marker sets to populations with a scattered intersection
#'
#' @description
#' Builds an ordered union of \eqn{T} marker IDs (`m1..mT`) and assigns to each of
#' `n_pops` a subset that preserves genomic order. A fraction `intersection` of
#' markers is shared by **all** populations (chosen at random positions and kept
#' in order). The remaining markers are split **disjointly** across populations
#' and placed at **random positions** (scattered) while preserving order.
#'
#' @param T integer (>= 1). Total markers in the union.
#' @param n_pops integer (>= 1). Number of populations.
#' @param intersection numeric in \[0,1\]. Fraction of markers shared by all populations.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{markers_union}: character vector \code{m1..mT}.
#'   \item \code{per_pop}: list of length \code{n_pops}; each is a character vector
#'         of marker IDs for that population, in genomic order.
#' }
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' str(assign_marker_sets(T = 20, n_pops = 3, intersection = 0.6))
#' }
#' @keywords internal
#' @noRd
assign_marker_sets <- function(T, n_pops, intersection = 1) {
  # Validate
  if (!is.numeric(T) || length(T) != 1L || T < 1 || !is.finite(T)) {
    stop("`T` must be a single finite integer >= 1.")
  }
  if (!is.numeric(n_pops) || length(n_pops) != 1L || n_pops < 1 || !is.finite(n_pops)) {
    stop("`n_pops` must be a single finite integer >= 1.")
  }
  if (!is.numeric(intersection) || length(intersection) != 1L ||
      !is.finite(intersection) || intersection < 0 || intersection > 1) {
    stop("`intersection` must be a single number in [0, 1].")
  }

  T <- as.integer(T); n_pops <- as.integer(n_pops)
  all_markers <- paste0("m", seq_len(T))

  # Trivial case: one population gets all markers
  if (n_pops == 1L) {
    return(list(markers_union = all_markers, per_pop = list(all_markers)))
  }

  # 1) Random core (present in all pops), positions sampled then sorted
  core_sz  <- max(0L, min(T, round(intersection * T)))
  core_idx <- if (core_sz > 0L) sort(sample.int(T, core_sz, replace = FALSE)) else integer(0)

  # 2) Remaining indices split disjointly & balanced across populations
  rem_idx <- setdiff(seq_len(T), core_idx)
  k <- length(rem_idx)
  if (k == 0L) {
    per_pop <- replicate(n_pops, all_markers, simplify = FALSE)
    return(list(markers_union = all_markers, per_pop = per_pop))
  }
  base  <- k %/% n_pops
  extra <- k - base * n_pops
  counts <- rep.int(base, n_pops)
  if (extra > 0L) counts[seq_len(extra)] <- counts[seq_len(extra)] + 1L

  rem_perm  <- if (k > 1L) sample(rem_idx, k, replace = FALSE) else rem_idx
  assign_id <- rep.int(seq_len(n_pops), times = counts)  # length k
  stopifnot(length(assign_id) == length(rem_perm))

  # 3) Build per-pop marker sets in ORIGINAL ORDER
  per_pop <- vector("list", n_pops)
  for (g in seq_len(n_pops)) {
    idx_g <- c(core_idx, rem_perm[assign_id == g])
    idx_g <- sort(idx_g)
    per_pop[[g]] <- all_markers[idx_g]
  }

  list(markers_union = all_markers, per_pop = per_pop)
}

#' Sample genotypes under Hardy-Weinberg equilibrium (0/1/2 coding)
#'
#' @description
#' Draws genotypes coded as A-allele counts \code{0,1,2} under HWE at allele
#' frequency \eqn{p}. Supports a scalar \eqn{p} with \code{n} draws, or a vector
#' of per-marker allele frequencies.
#'
#' @param p numeric scalar or vector in (0,1). A-allele frequency.
#' @param n optional integer. Number of draws when \code{p} is scalar. Ignored when
#'   \code{p} is a vector.
#'
#' @return Integer vector of length \code{n} (if \code{p} scalar) or
#'   \code{length(p)} (if \code{p} vector) with values in \code{0,1,2}.
#'
#' @examples
#' \dontrun{
#' set.seed(7)
#' rHWE_geno012(0.3, n = 5)
#' rHWE_geno012(c(0.2, 0.5, 0.8))
#' }
#' @keywords internal
#' @noRd
rHWE_geno012 <- function(p, n = NULL) {
  if (length(p) == 1L && !is.null(n)) {
    if (!is.numeric(n) || length(n) != 1L || n < 1 || !is.finite(n)) {
      stop("`n` must be a single finite integer >= 1 when `p` is scalar.")
    }
    p <- rep(p, as.integer(n))
  } else if (length(p) >= 1L && is.null(n)) {
    # ok, vector p
  } else {
    stop("Provide either scalar `p` with `n`, or a vector `p` with `n = NULL`.")
  }

  p <- pmin(pmax(as.numeric(p), 1e-6), 1 - 1e-6)
  aa <- (1 - p)^2
  Aa <- 2 * p * (1 - p)
  u  <- stats::runif(length(p))

  # 0 with prob aa, 1 with prob Aa, 2 otherwise
  out <- ifelse(u < aa, 0L, ifelse(u < aa + Aa, 1L, 2L))
  as.integer(out)
}

#' Simulate maternal genotypes per marker under HWE
#'
#' @description
#' Convenience wrapper around \code{rHWE_geno012()} to draw the known maternal
#' genotype vector (0=aa, 1=Aa, 2=AA) for a population at per-marker allele
#' frequencies \code{p_m}. Returns a vector with names preserved (if any).
#'
#' @param p_m numeric vector in (0,1). Maternal A-allele frequencies per marker.
#'
#' @return Integer vector of the same length (and names) as \code{p_m}, with entries
#'   \code{0/1/2}.
#'
#' @examples
#' \dontrun{
#' set.seed(3)
#' simulate_maternal_geno(p_m = rep(0.5, 5))
#' }
#' @keywords internal
#' @noRd
simulate_maternal_geno <- function(p_m) {
  g <- rHWE_geno012(p = p_m)
  names(g) <- names(p_m)
  g
}

#' Convert 0/1/2 genotype to one-hot paternal mixture columns (AA/Aa/aa)
#'
#' @description
#' Converts a genotype vector coded as A-allele counts \code{2,1,0} (i.e.,
#' \code{2 = AA}, \code{1 = Aa}, \code{0 = aa}) into a \code{3 x T} matrix of
#' one-hot columns with row names \code{c("AA","Aa","aa")}. \code{NA} inputs yield
#' \code{NA}-filled columns.
#'
#' @param g integer vector with values in \code{2,1,0,NA}.
#'
#' @return Numeric matrix \code{3 x length(g)} with one-hot columns (or \code{NA}s),
#'   rownames \code{AA,Aa,aa}.
#'
#' @examples
#' \dontrun{
#' geno012_to_pi(c(2, 1, 0, NA, 2))
#' }
#' @keywords internal
#' @noRd
geno012_to_pi <- function(g) {
  g <- as.integer(g)
  Tm <- length(g)
  out <- matrix(NA_real_, 3L, Tm, dimnames = list(c("AA","Aa","aa"), NULL))
  if (Tm == 0L) return(out)

  # Fill columns
  for (t in seq_len(Tm)) {
    gt <- g[t]
    if (is.na(gt)) next
    if (gt == 2L)      out[, t] <- c(1, 0, 0)
    else if (gt == 1L) out[, t] <- c(0, 1, 0)
    else if (gt == 0L) out[, t] <- c(0, 0, 1)
    else stop("`g` must contain only {0,1,2,NA}.")
  }
  out
}

#' Collapse a pedigree into nuclear-family populations (by mother)
#'
#' @description
#' Given a pedigree data frame with columns \code{id}, \code{mother}, \code{father},
#' \code{generation}, and \code{family_id} (equal to the mother ID for each family),
#' returns a list of populations, each describing one nuclear family (one known
#' mother, optional known sire, and the list of offspring IDs).
#'
#' If multiple distinct father IDs appear within a family, \code{father_id} is set
#' to \code{NA_character_}.
#'
#' @param ped A data frame with at least the columns
#'   \code{c("id","mother","father","generation","family_id")}.
#'
#' @return A named list of populations; each element is a \code{list} with fields:
#' \itemize{
#'   \item \code{pop_id}: character, same as \code{family_id}.
#'   \item \code{mother_id}: character.
#'   \item \code{father_id}: character or \code{NA_character_}.
#'   \item \code{offspring_ids}: character vector.
#' }
#'
#' @examples
#' \dontrun{
#' ped <- data.frame(
#'   id = c("M1","D1","K1","K2"),
#'   mother = c(NA, NA, "M1", "M1"),
#'   father = c(NA, NA, "D1", "D1"),
#'   generation = c(1,1,2,2),
#'   family_id = c("M1","D1","M1","M1"),
#'   stringsAsFactors = FALSE
#' )
#' str(families_from_pedigree(ped))
#' }
#' @keywords internal
#' @noRd
families_from_pedigree <- function(ped) {
  if (!is.data.frame(ped)) stop("`ped` must be a data.frame.")
  req <- c("id", "mother", "father", "generation", "family_id")
  if (!all(req %in% names(ped))) {
    stop("`ped` must contain columns: ", paste(req, collapse = ", "))
  }

  # Normalize types
  ped$id         <- as.character(ped$id)
  ped$mother     <- as.character(ped$mother)
  ped$father     <- as.character(ped$father)
  ped$family_id  <- as.character(ped$family_id)

  mothers <- unique(stats::na.omit(ped$family_id))
  pops <- lapply(mothers, function(m) {
    is_child <- !is.na(ped$mother) & ped$mother == m
    kids <- ped$id[is_child]
    dads <- unique(stats::na.omit(ped$father[is_child]))
    father_id <- if (length(dads) == 1L) dads else NA_character_
    list(pop_id = m, mother_id = m, father_id = father_id, offspring_ids = kids)
  })
  names(pops) <- mothers
  pops
}

