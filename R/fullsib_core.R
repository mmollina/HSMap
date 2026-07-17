# Known-sire full-sib mapping driver (oracle phase). See dev/known_sire_design.md.

# Per-interval genetic distance and status for a sex-specific consensus map, mirroring
# the OP block reporter: distance is inv_haldane(r) ONLY for a linked interval; an
# interval with r at/above gap_r is a no_linkage_boundary (gap, distance NA), and one
# with no informative meioses is insufficient_information (distance NA). r = 0.5 (no
# linkage) is therefore never turned into a huge finite map distance.
.fs_interval_report <- function(r, tot, gap_r) {
  r <- as.numeric(r); Ti <- length(r)
  status <- rep("linked", Ti)
  bad_tot <- is.na(tot) | tot <= 0                 # Inf (dispatch sentinel) counts as ample data
  status[bad_tot | is.na(r)] <- "insufficient_information"
  status[status == "linked" & r >= gap_r] <- "no_linkage_boundary"
  dist <- ifelse(status == "linked", inv_haldane(pmin(r, 0.5 - 1e-12)), NA_real_)
  list(dist = dist, status = status)
}

# Per-interval structural informativeness of the maternal and paternal chains, and a
# global label-exchangeability flag. A parent's transmission is informative at interval k
# only if that parent is HETEROZYGOUS at both markers k and k+1 (a homozygous flanking
# marker makes the transmitted homolog indistinguishable). Maternal and paternal labels
# are non-identifiable (exchangeable) when every cross has identical maternal and paternal
# haplotypes, so swapping (r_m, r_p) leaves the likelihood invariant.
.fs_informativeness <- function(mat_builts, pat_builts, Ti) {
  mat_inf <- logical(Ti); pat_inf <- logical(Ti); n_off <- integer(Ti)
  het <- function(A) A[1, ] != A[2, ]
  for (b in mat_builts) {
    mh <- het(b$Am); nk <- nrow(b$G)
    for (k in seq_len(Ti)) { n_off[k] <- n_off[k] + nk; if (mh[k] && mh[k + 1L]) mat_inf[k] <- TRUE }
  }
  for (b in pat_builts) {
    ph <- het(b$Ap)
    for (k in seq_len(Ti)) if (ph[k] && ph[k + 1L]) pat_inf[k] <- TRUE
  }
  exch <- length(pat_builts) > 0L && length(mat_builts) == length(pat_builts) &&
    all(vapply(pat_builts, function(b) !is.null(b$Ap) && identical(b$Am, b$Ap), logical(1)))
  list(n_off = n_off, mat_inf = mat_inf, pat_inf = pat_inf, globally_exchangeable = exch)
}

# Overlay structural informativeness / exchangeability onto an r/gap interval report:
# a structurally uninformative interval becomes insufficient_information (NA distance);
# under global exchangeability the doubly-informative intervals become
# nonidentifiable_exchangeable (NA distance) rather than forcing a per-parent estimate.
.fs_apply_informativeness <- function(rep, informative, exch, both_inf) {
  rep$status[!informative] <- "insufficient_information"
  rep$dist[!informative] <- NA_real_
  if (isTRUE(exch)) { rep$status[both_inf] <- "nonidentifiable_exchangeable"; rep$dist[both_inf] <- NA_real_ }
  rep
}

# Shared numerical-argument validation for the full-sib / mixed fitters.
.fs_validate_inputs <- function(epsilon, tol, maxit, r_start, lambda = NULL, q0 = NULL) {
  num1 <- function(v) is.numeric(v) && length(v) == 1L && !is.na(v)
  if (!num1(epsilon) || epsilon < 0 || epsilon >= 1)
    stop("`epsilon` must be a single number in [0, 1).", call. = FALSE)
  if (!num1(tol) || tol <= 0)
    stop("`tol` must be a single positive number.", call. = FALSE)
  if (!num1(maxit) || maxit < 1 || maxit != round(maxit))
    stop("`maxit` must be a single positive integer.", call. = FALSE)
  if (!num1(r_start) || r_start < 0 || r_start > 0.5)
    stop("`r_start` must be a single number in [0, 0.5].", call. = FALSE)
  if (!is.null(lambda) && (!num1(lambda) || lambda < 0))
    stop("`lambda` must be a single number >= 0.", call. = FALSE)
  if (!is.null(q0) && (!num1(q0) || q0 < 0 || q0 > 1))
    stop("`q0` must be a single number in [0, 1].", call. = FALSE)
  invisible(TRUE)
}

# Material (scale-aware) decrease of the active objective. Used both inside the loop and
# for the final-step check (objective at the final returned params vs the preceding one).
.fs_obj_decreased <- function(new, old, rel = 1e-8) {
  is.finite(old) && is.finite(new) && new < old - rel * (1 + abs(old))
}

# Build a cross's maternal/paternal 2xT allele matrices from parent genotypes and
# Resolve a parent's oracle 2 x length(order) allele matrix (entries 0/1) from EITHER
# an explicit haplotype matrix (the canonical input) OR a fully-resolved legacy phase
# vector. Unresolved (NA) phase or NA haplotype entries are rejected -- never coerced to
# coupling. The same parent id always resolves to the same matrix (shared object).
.resolve_parent_alleles <- function(id, role, haplotypes, phased, parent_geno, order, fn) {
  z <- length(order)
  H <- if (!is.null(haplotypes)) haplotypes[[id]] else NULL
  if (!is.null(H)) {                                   # canonical: explicit haplotype matrix
    H <- as.matrix(H)
    if (nrow(H) != 2L)
      stop(fn, "(): haplotypes for ", role, " '", id, "' must be a 2 x T matrix.", call. = FALSE)
    if (!is.null(colnames(H))) {
      if (!all(order %in% colnames(H)))
        stop(fn, "(): haplotypes for ", role, " '", id, "' are missing markers in `order`.", call. = FALSE)
      H <- H[, order, drop = FALSE]
    } else if (ncol(H) != z) {
      stop(fn, "(): haplotypes for ", role, " '", id, "' have ", ncol(H),
           " columns but `order` has ", z, " and the matrix has no column names.", call. = FALSE)
    }
    if (anyNA(H))
      stop(fn, "(): haplotypes for ", role, " '", id, "' contain NA (unresolved) at fitted ",
           "marker(s); the core requires fully resolved oracle haplotypes.", call. = FALSE)
    if (!all(H %in% c(0L, 1L)))
      stop(fn, "(): haplotype entries for ", role, " '", id, "' must be 0/1.", call. = FALSE)
    return(matrix(as.integer(H), 2))
  }
  pv <- if (!is.null(phased)) phased[[id]] else NULL   # legacy: adjacent phase vector
  if (is.null(pv))
    stop(fn, "(): no oracle phase for ", role, " '", id, "'. Supply haplotypes or a phase vector.", call. = FALSE)
  if (length(pv) != z - 1L)
    stop(fn, "(): phase vector for ", role, " '", id, "' must have length length(order) - 1.", call. = FALSE)
  if (anyNA(pv))
    stop(fn, "(): unresolved (NA) phase supplied for ", role, " '", id, "' via the legacy ",
         "phase-vector interface; NA phase is never treated as coupling. Provide an explicit ",
         "haplotype matrix (haplotypes_", substr(role, 1, 1), ") instead.", call. = FALSE)
  g <- parent_geno[[id]]
  if (is.null(g) || all(is.na(g)))
    stop(fn, "(): ", role, " '", id, "' has no genotype.", call. = FALSE)
  names(g) <- names(g) %||% order; g <- g[order]
  H <- phase_to_haplotypes(g, pv)
  if (anyNA(H))
    stop(fn, "(): missing required ", role, " genotype at fitted marker(s) for '", id,
         "'. The core does not invent a parental allele.", call. = FALSE)
  matrix(as.integer(H), 2)
}

# Build a cross's maternal + paternal allele matrices (from haplotypes or phase).
.fs_build_cross <- function(cr, parent_geno, haplotypes_m, haplotypes_p, phased_m, phased_p, order, on_missing) {
  Am <- .resolve_parent_alleles(cr$mother_id, "mother", haplotypes_m, phased_m, parent_geno, order, "hmm_map_fullsib")
  Ap <- .resolve_parent_alleles(cr$father_id, "sire",   haplotypes_p, phased_p, parent_geno, order, "hmm_map_fullsib")
  G <- cr$G[, order, drop = FALSE]; storage.mode(G) <- "integer"
  list(G = G, Am = Am, Ap = Ap)
}

#' Fit a known-sire full-sib recombination map (oracle phase)
#'
#' @description
#' Estimate a shared maternal recombination map \code{r_m} and a shared paternal
#' recombination map \code{r_p} from genotyped known-sire full-sib crosses, using a
#' four-state hidden inheritance process (maternal homolog, paternal homolog) with the
#' parental phase held fixed (oracle phase). Maternal and paternal recombination
#' fractions are estimated separately and are never assumed equal. See
#' \code{dev/known_sire_design.md}.
#'
#' @param x An \code{HSMap.data} object with cross-aware fields (\code{crosses},
#'   \code{parent_genotypes}).
#' @param haplotypes_m,haplotypes_p \strong{Canonical oracle input.} Named lists
#'   mapping each mother / sire id to a fully resolved \code{2 x T} haplotype allele
#'   matrix (entries 0/1, rows = homolog 1/2), with column names giving the markers
#'   (or exactly \code{length(order)} columns in \code{order}). This is preferred over
#'   \code{phased_*} because it is unambiguous when heterozygous markers are separated
#'   by homozygous markers. \code{NA} (unresolved) entries are rejected. The same
#'   parent id resolves to the same matrix across crosses.
#' @param phased_m,phased_p \emph{Legacy compatibility only.} Named lists mapping each
#'   mother / sire id to a fully resolved adjacent phase vector (length \code{z-1},
#'   1 = coupling, 0 = repulsion). \strong{Any \code{NA} (unresolved) phase is rejected}
#'   -- unresolved phase is never treated as coupling. Used only when the corresponding
#'   \code{haplotypes_*} entry is absent. Prefer \code{haplotypes_*}.
#' @param crosses Optional character vector of \code{cross_id}s to fit; default all
#'   \code{known_sire_genotyped} crosses.
#' @param order Optional marker order (character); default the common column order of
#'   the selected crosses' genotype matrices.
#' @param epsilon Symmetric genotyping error rate.
#' @param tol Convergence tolerance on the relative active objective and on \code{r}.
#' @param maxit Maximum EM iterations.
#' @param r_start Initial recombination fraction for both maps.
#' @param gap_r Recombination fraction at/above which an interval is reported as a
#'   no-linkage boundary (a gap): its map distance is \code{NA}, never a large finite
#'   value. Default \code{0.499}.
#' @param on_missing One of \code{"error"} (default): stop on a missing required parent
#'   genotype at a fitted marker.
#' @param return_posterior Logical; if \code{TRUE}, include per-cross posterior
#'   inheritance probabilities.
#'
#' @return An object of class \code{HSMap.fullsib}. It reports a \strong{sex-specific
#'   consensus} map (\code{map_scope = "sex_specific_consensus"}): one maternal map
#'   \code{r_m} pooling all mothers and one paternal map \code{r_p} pooling all genotyped
#'   sires -- NOT parent-specific maps. \code{$fit$maternal_meioses_by_mother} and
#'   \code{$fit$paternal_meioses_by_sire} give the per-interval meiosis counts pooled
#'   into each map. See also \code{$fit} for r/d/interval status and traces.
#' @export
hmm_map_fullsib <- function(x, phased_m = NULL, phased_p = NULL,
                            haplotypes_m = NULL, haplotypes_p = NULL,
                            crosses = NULL, order = NULL,
                            epsilon = 0.05, tol = 1e-6, maxit = 1000L, r_start = 0.1,
                            gap_r = 0.499, on_missing = c("error"), return_posterior = FALSE) {
  if (!inherits(x, "HSMap.data")) stop("`x` must be an HSMap.data object.")
  if (is.null(x$crosses)) stop("`x` lacks cross-aware fields; re-read with read_HSMap_data().")
  on_missing <- match.arg(on_missing)
  .fs_validate_inputs(epsilon, tol, maxit, r_start)
  all_cross <- x$crosses
  fs_ids <- names(all_cross)[vapply(all_cross, function(c) identical(c$family_type, "known_sire_genotyped"), logical(1))]
  if (!is.null(crosses)) fs_ids <- intersect(crosses, fs_ids)
  if (!length(fs_ids)) stop("hmm_map_fullsib(): no known_sire_genotyped crosses to fit.")

  if (is.null(order)) order <- colnames(all_cross[[fs_ids[1]]]$G)
  order <- as.character(order)
  for (cid in fs_ids) {
    oc <- colnames(all_cross[[cid]]$G)
    if (!all(order %in% oc))
      stop("hmm_map_fullsib(): cross '", cid, "' is missing markers from `order`.")
  }
  z <- length(order); Ti <- z - 1L
  if (z < 2L) stop("Need >= 2 markers.")

  # precompute per-cross aligned data + allele matrices (phase fixed)
  built <- lapply(fs_ids, function(cid)
    .fs_build_cross(all_cross[[cid]], x$parent_genotypes, haplotypes_m, haplotypes_p, phased_m, phased_p, order, on_missing))
  names(built) <- fs_ids
  moms <- unique(vapply(fs_ids, function(cid) all_cross[[cid]]$mother_id, character(1)))
  sires <- unique(vapply(fs_ids, function(cid) all_cross[[cid]]$father_id, character(1)))

  r_m <- rep(r_start, Ti); r_p <- rep(r_start, Ti)
  clampr <- function(r) pmin(pmax(r, 0.0), 0.5 - 1e-12)

  ll_trace <- numeric(0); dm_trace <- numeric(0); dp_trace <- numeric(0)
  prev_ll <- -Inf; obj_dec <- FALSE
  converged <- FALSE; conv_reason <- "maxit"; it <- 0L
  OBJ_DEC_REL <- 1e-8
  for (it in seq_len(maxit)) {
    m_sw <- numeric(Ti); p_sw <- numeric(Ti); tot <- numeric(Ti); ll <- 0
    for (cid in fs_ids) {
      b <- built[[cid]]
      es <- fs_estep_cpp(b$G, b$Am, b$Ap, r_m, r_p, epsilon, FALSE)
      m_sw <- m_sw + es$m_switch; p_sw <- p_sw + es$p_switch; tot <- tot + es$total
      ll <- ll + es$loglik
    }
    ll_trace <- c(ll_trace, ll)                      # active objective at the current params
    if (it > 1L && .fs_obj_decreased(ll, prev_ll)) obj_dec <- TRUE
    r_m_new <- clampr(ifelse(tot > 0, m_sw / tot, r_m))
    r_p_new <- clampr(ifelse(tot > 0, p_sw / tot, r_p))
    dm <- max(abs(r_m_new - r_m)); dp <- max(abs(r_p_new - r_p))
    dm_trace <- c(dm_trace, dm); dp_trace <- c(dp_trace, dp)
    rel <- if (is.finite(prev_ll)) abs(ll - prev_ll) / (1 + abs(prev_ll)) else Inf
    r_m <- r_m_new; r_p <- r_p_new
    if (it > 1L && rel < tol && max(dm, dp) < tol) {
      converged <- TRUE; conv_reason <- "relative_objective_and_params_stable"; prev_ll <- ll; break }
    prev_ll <- ll
  }
  # objective recomputed at the FINAL returned parameters (after the last M-step); the
  # traces end exactly at this value, and a decrease at the final step is caught too.
  # Per-interval meiosis counts are also attributed to each mother / sire to make the
  # consensus pooling explicit.
  inm <- paste0(order[-z], "-", order[-1])
  final_ll <- 0; final_tot <- numeric(Ti)
  mat_by_mom <- matrix(0, Ti, length(moms), dimnames = list(inm, moms))
  pat_by_sire <- matrix(0, Ti, length(sires), dimnames = list(inm, sires))
  for (cid in fs_ids) { b <- built[[cid]]
    es <- fs_estep_cpp(b$G, b$Am, b$Ap, r_m, r_p, epsilon, FALSE)
    final_ll <- final_ll + es$loglik; final_tot <- final_tot + es$total
    mat_by_mom[, all_cross[[cid]]$mother_id] <- mat_by_mom[, all_cross[[cid]]$mother_id] + es$total
    pat_by_sire[, all_cross[[cid]]$father_id] <- pat_by_sire[, all_cross[[cid]]$father_id] + es$total }
  if (.fs_obj_decreased(final_ll, prev_ll)) obj_dec <- TRUE     # decrease at the final M-step
  ll_trace <- c(ll_trace, final_ll); dm_trace <- c(dm_trace, 0); dp_trace <- c(dp_trace, 0)
  if (obj_dec) { converged <- FALSE; conv_reason <- "objective_decreased" }
  if (!converged) {
    if (identical(conv_reason, "objective_decreased"))
      warning("hmm_map_fullsib(): the active EM objective decreased materially; ",
              "results may be unreliable.", call. = FALSE)
    else
      warning("hmm_map_fullsib(): EM did not converge in ", it, " iterations; the ",
              "returned estimates are the last iterate. Increase `maxit` or relax `tol`.",
              call. = FALSE)
  }

  post <- NULL
  if (isTRUE(return_posterior)) {
    post <- lapply(fs_ids, function(cid) { b <- built[[cid]]
      es <- fs_estep_cpp(b$G, b$Am, b$Ap, r_m, r_p, epsilon, TRUE); es$gamma })
    names(post) <- fs_ids
  }

  info <- .fs_informativeness(built, built, Ti)
  both_inf <- info$mat_inf & info$pat_inf
  rep_m <- .fs_apply_informativeness(.fs_interval_report(r_m, final_tot, gap_r), info$mat_inf, info$globally_exchangeable, both_inf)
  rep_p <- .fs_apply_informativeness(.fs_interval_report(r_p, final_tot, gap_r), info$pat_inf, info$globally_exchangeable, both_inf)
  ident_tab <- data.frame(
    interval = inm, n_offspring = info$n_off,
    maternal_informative = info$mat_inf, paternal_informative = info$pat_inf,
    exchangeable = info$globally_exchangeable & both_inf,
    status_m = rep_m$status, status_p = rep_p$status, stringsAsFactors = FALSE)
  fit <- list(
    r_m = stats::setNames(r_m, inm),
    r_p = stats::setNames(r_p, inm),
    d_m = stats::setNames(rep_m$dist, inm), d_p = stats::setNames(rep_p$dist, inm),
    interval_status_m = stats::setNames(rep_m$status, inm),
    interval_status_p = stats::setNames(rep_p$status, inm),
    identifiability = ident_tab, identifiable_labels = !info$globally_exchangeable,
    gap_r = gap_r, meiosis_count = stats::setNames(final_tot, inm),
    # sex-specific CONSENSUS maps: r_m pools all maternal meioses across mothers, r_p
    # pools all paternal meioses across genotyped sires. These are NOT parent-specific
    # maps; a repeated mother/sire contributes to one shared map (see design note).
    map_scope = "sex_specific_consensus",
    maternal_meioses_by_mother = mat_by_mom,
    paternal_meioses_by_sire = pat_by_sire,
    logLik = final_ll, objective = final_ll,
    converged = converged, conv_reason = conv_reason, iters = it,
    objective_decreased = obj_dec,
    loglik_trace = ll_trace, objective_trace = ll_trace,   # no penalty for pure full-sib
    max_dr_m_trace = dm_trace, max_dr_p_trace = dp_trace,
    epsilon = epsilon, posterior = post
  )
  out <- list(
    order = order, fit = fit, map_scope = "sex_specific_consensus",
    contributing_crosses = fs_ids,
    contributing_mothers = moms, contributing_sires = sires,
    maternal_crosses = fs_ids, paternal_crosses = fs_ids,
    family_type = stats::setNames(rep("known_sire_genotyped", length(fs_ids)), fs_ids),
    parent_phase = list(maternal = phased_m[moms], paternal = phased_p[sires])
  )
  class(out) <- "HSMap.fullsib"
  out
}

# Build maternal 2xT allele matrix for an OP cross (phase in emission).
.op_build_cross <- function(cr, parent_geno, haplotypes_m, phased_m, order, on_missing) {
  Am <- .resolve_parent_alleles(cr$mother_id, "mother", haplotypes_m, phased_m, parent_geno, order, "hmm_map_mixed")
  G <- cr$G[, order, drop = FALSE]; storage.mode(G) <- "integer"
  list(G = G, Am = Am, mother = cr$mother_id)
}

#' Fit a mixed open-pollinated + known-sire recombination map
#'
#' @description
#' Estimate a shared maternal recombination map from a dataset containing both
#' open-pollinated (OP) and known-sire full-sib crosses, plus a separate shared
#' paternal map from the full-sib crosses. OP crosses use the two-state maternal model
#' with a dam-specific paternal gametic frequency \code{q}; full-sib crosses use the
#' four-state model. With only OP crosses present, dispatches to the existing engine so
#' OP-only results are unchanged. See \code{dev/known_sire_design.md}.
#'
#' @param x An \code{HSMap.data} object with cross-aware fields.
#' @param haplotypes_m,haplotypes_p Canonical oracle input: named lists of fully
#'   resolved \code{2 x T} parental haplotype allele matrices (see
#'   \code{\link{hmm_map_fullsib}}). \code{haplotypes_p} (or \code{phased_p}) is
#'   required when any full-sib cross is fitted.
#' @param phased_m,phased_p Legacy compatibility: fully resolved adjacent phase vectors
#'   (mother/sire id -> length z-1). \code{NA} phase is rejected. For the open-pollinated
#'   dispatch a resolved \code{phased_m} is used directly; explicit haplotypes are
#'   converted where possible.
#' @param order Optional marker order; default the common genotype column order.
#' @param epsilon,tol,maxit,r_start As in \code{hmm_map_fullsib}.
#' @param lambda,q0 OP paternal-\code{q} pseudocount total and shrinkage target
#'   (MAP update), matching \code{hmm_map}.
#' @param gap_r Recombination fraction at/above which an interval is a no-linkage
#'   boundary (gap): map distance \code{NA}. Default \code{0.499}.
#' @param untyped_sire How to treat \code{known_sire_untyped} crosses: \code{"error"}
#'   (default) stop; \code{"open_pollinated"} model them with the OP model (reported).
#' @param on_missing \code{"error"} (default) on a missing required parent genotype.
#' @param return_posterior Logical; include full-sib posterior inheritance probabilities.
#'
#' @return An object of class \code{HSMap.mixed} reporting a \strong{sex-specific
#'   consensus} map (\code{map_scope = "sex_specific_consensus"}): one maternal map
#'   pooling ALL crosses (open-pollinated + full-sib) and one paternal map pooling the
#'   genotyped sires; \code{maternal_crosses} / \code{paternal_crosses} and the
#'   per-interval \code{maternal_meioses_by_mother} / \code{paternal_meioses_by_sire}
#'   show exactly what is pooled. Not parent-specific.
#' @export
hmm_map_mixed <- function(x, phased_m = NULL, phased_p = NULL,
                          haplotypes_m = NULL, haplotypes_p = NULL, order = NULL,
                          epsilon = 0.05, lambda = 20, q0 = 0.5, tol = 1e-6,
                          maxit = 1000L, r_start = 0.1, gap_r = 0.499,
                          untyped_sire = c("error", "open_pollinated"),
                          on_missing = c("error"), return_posterior = FALSE) {
  if (!inherits(x, "HSMap.data")) stop("`x` must be an HSMap.data object.")
  if (is.null(x$crosses)) stop("`x` lacks cross-aware fields; re-read with read_HSMap_data().")
  untyped_sire <- match.arg(untyped_sire); on_missing <- match.arg(on_missing)
  .fs_validate_inputs(epsilon, tol, maxit, r_start, lambda, q0)
  cx <- x$crosses
  ftype <- vapply(cx, function(c) c$family_type, character(1))

  untyped <- names(cx)[ftype == "known_sire_untyped"]
  if (length(untyped)) {
    if (identical(untyped_sire, "error"))
      stop("hmm_map_mixed(): ", length(untyped), " cross(es) name an untyped sire ",
           "(e.g. ", paste(utils::head(untyped, 3), collapse = ", "), "). Set ",
           "untyped_sire = 'open_pollinated' to model them with the OP model, or ",
           "supply the sire genotypes.", call. = FALSE)
  }
  op_ids <- names(cx)[ftype == "open_pollinated" |
                        (ftype == "known_sire_untyped" & untyped_sire == "open_pollinated")]
  fs_ids <- names(cx)[ftype == "known_sire_genotyped"]
  if (!length(op_ids) && !length(fs_ids)) stop("hmm_map_mixed(): no fittable crosses.")

  if (is.null(order)) order <- colnames(cx[[c(fs_ids, op_ids)[1]]]$G)
  order <- as.character(order); z <- length(order); Ti <- z - 1L

  # ---- OP-only: dispatch to the existing engine (unchanged results) ---------
  if (!length(fs_ids)) {
    # the legacy allele-state engine consumes an adjacent phase vector; derive one from
    # explicit haplotypes if needed and reject an unresolved (NA) result.
    op_phase <- function(mom) {
      pv <- if (!is.null(phased_m)) phased_m[[mom]] else NULL
      if (is.null(pv) && !is.null(haplotypes_m) && !is.null(haplotypes_m[[mom]]))
        pv <- haplotypes_to_phase(as.matrix(haplotypes_m[[mom]]))
      if (is.null(pv)) stop("hmm_map_mixed(): no maternal phase for OP mother '", mom, "'.", call. = FALSE)
      if (anyNA(pv))
        stop("hmm_map_mixed(): unresolved (NA) phase for OP mother '", mom, "'; the legacy ",
             "open-pollinated dispatch needs a fully resolved phase vector (homozygous markers ",
             "between heterozygous ones can leave phase undefined). Supply a resolved phase_m.",
             call. = FALSE)
      as.integer(pv)
    }
    phased_list <- lapply(op_ids, function(cid)
      structure(list(dam = cx[[cid]]$mother_id, order = order,
                     phase_vec = op_phase(cx[[cid]]$mother_id)),
                class = "HSMap.phased"))
    names(phased_list) <- vapply(op_ids, function(cid) cx[[cid]]$mother_id, character(1))
    if (length(phased_list) == 1L) {
      res <- hmm_map(x, phased = phased_list[[1]], dam = names(phased_list)[1],
                     epsilon = epsilon, tol = tol, lambda = lambda, maxit = maxit, r_start = r_start)
    } else {
      ph <- phased_list; class(ph) <- "HSMap.phased.multi"
      res <- hmm_map(x, phased = ph, dam = "all", epsilon = epsilon, tol = tol,
                     lambda = lambda, maxit = maxit, r_start = r_start, method = "joint")
    }
    rfit <- res$fit
    inm <- paste0(order[-z], "-", order[-1])
    rep_m <- .fs_interval_report(as.numeric(rfit$r), rep(Inf, Ti), gap_r)  # OP r already fitted
    out <- list(order = order, fit = list(
      r_m = rfit$r, r_p = NULL,
      d_m = stats::setNames(rep_m$dist, inm), d_p = NULL,
      interval_status_m = stats::setNames(rep_m$status, inm), interval_status_p = NULL,
      gap_r = gap_r, map_scope = "maternal_consensus",   # OP-only: one maternal map, no paternal
      logLik = rfit$logLik, objective = rfit$objective %||% rfit$penalized_obj,
      converged = rfit$converged, conv_reason = rfit$conv_reason, iters = rfit$iters,
      objective_decreased = rfit$objective_decreased %||% FALSE,
      loglik_trace = rfit$loglik_trace, epsilon = epsilon, q = rfit$q),
      contributing_crosses = op_ids,
      contributing_mothers = unique(names(phased_list)), contributing_sires = character(0),
      family_type = stats::setNames(rep("open_pollinated", length(op_ids)), op_ids),
      dispatched = TRUE, op_result = res)
    class(out) <- "HSMap.mixed"
    return(out)
  }

  # ---- mixed EM (OP + full-sib) ---------------------------------------------
  if (length(fs_ids) && is.null(phased_p) && is.null(haplotypes_p))
    stop("hmm_map_mixed(): full-sib crosses present but no sire phase supplied ",
         "(provide haplotypes_p or phased_p).", call. = FALSE)
  fs_built <- lapply(fs_ids, function(cid)
    .fs_build_cross(cx[[cid]], x$parent_genotypes, haplotypes_m, haplotypes_p, phased_m, phased_p, order, on_missing))
  names(fs_built) <- fs_ids
  op_built <- lapply(op_ids, function(cid)
    .op_build_cross(cx[[cid]], x$parent_genotypes, haplotypes_m, phased_m, order, on_missing))
  names(op_built) <- op_ids

  # OP crosses grouped by dam (one q per dam)
  op_dam <- vapply(op_ids, function(cid) cx[[cid]]$mother_id, character(1))
  dams_op <- unique(op_dam)
  q_list <- stats::setNames(lapply(dams_op, function(d) rep(q0, z)), dams_op)
  alpha_p <- lambda * q0; beta_p <- lambda * (1 - q0)

  r_m <- rep(r_start, Ti); r_p <- rep(r_start, Ti)
  clampr <- function(r) pmin(pmax(r, 0), 0.5 - 1e-12)
  q_pen <- function(qv) if (lambda > 0) sum(alpha_p * log(pmin(pmax(qv,1e-12),1-1e-12)) +
                                            beta_p * log(pmin(pmax(1-qv,1e-12),1-1e-12))) else 0

  obj_trace <- numeric(0); ll_trace <- numeric(0)
  dm_trace <- numeric(0); dp_trace <- numeric(0); dq_trace <- numeric(0)
  prev_obj <- -Inf; obj_dec <- FALSE
  converged <- FALSE; conv_reason <- "maxit"; it <- 0L; OBJ_DEC_REL <- 1e-8
  for (it in seq_len(maxit)) {
    m_sw <- numeric(Ti); m_tot <- numeric(Ti); p_sw <- numeric(Ti); p_tot <- numeric(Ti)
    ll <- 0
    NA_d <- stats::setNames(lapply(dams_op, function(d) numeric(z)), dams_op)
    Na_d <- stats::setNames(lapply(dams_op, function(d) numeric(z)), dams_op)
    for (cid in fs_ids) { b <- fs_built[[cid]]
      es <- fs_estep_cpp(b$G, b$Am, b$Ap, r_m, r_p, epsilon, FALSE)
      m_sw <- m_sw + es$m_switch; m_tot <- m_tot + es$total
      p_sw <- p_sw + es$p_switch; p_tot <- p_tot + es$total; ll <- ll + es$loglik }
    for (cid in op_ids) { b <- op_built[[cid]]; d <- b$mother
      es <- op_estep_cpp(b$G, b$Am, r_m, q_list[[d]], epsilon)
      m_sw <- m_sw + es$m_switch; m_tot <- m_tot + es$total
      NA_d[[d]] <- NA_d[[d]] + es$N_A; Na_d[[d]] <- Na_d[[d]] + es$N_a; ll <- ll + es$loglik }

    pen <- sum(vapply(q_list, q_pen, numeric(1)))
    obj <- ll + pen; obj_trace <- c(obj_trace, obj); ll_trace <- c(ll_trace, ll)
    if (it > 1L && .fs_obj_decreased(obj, prev_obj)) obj_dec <- TRUE

    r_m_new <- clampr(ifelse(m_tot > 0, m_sw / m_tot, r_m))
    r_p_new <- clampr(ifelse(p_tot > 0, p_sw / p_tot, r_p))
    q_new <- q_list; dq <- 0
    for (d in dams_op) {
      qd <- (NA_d[[d]] + alpha_p) / (NA_d[[d]] + Na_d[[d]] + alpha_p + beta_p)
      qd <- pmin(pmax(qd, 1e-9), 1 - 1e-9); dq <- max(dq, max(abs(qd - q_list[[d]]))); q_new[[d]] <- qd
    }
    dm <- max(abs(r_m_new - r_m)); dp <- max(abs(r_p_new - r_p))
    dm_trace <- c(dm_trace, dm); dp_trace <- c(dp_trace, dp); dq_trace <- c(dq_trace, dq)
    rel <- if (is.finite(prev_obj)) abs(obj - prev_obj) / (1 + abs(prev_obj)) else Inf
    r_m <- r_m_new; r_p <- r_p_new; q_list <- q_new
    if (it > 1L && rel < tol && max(dm, dp) < tol && dq < tol) {
      converged <- TRUE; conv_reason <- "relative_objective_and_params_stable"; prev_obj <- obj; break }
    prev_obj <- obj
  }
  # objective recomputed at the FINAL returned parameters:
  # final_objective = final_logLik + final_q_penalty. Traces end exactly here. Per-
  # interval meiosis counts are attributed to each mother / sire to make the consensus
  # pooling explicit (maternal pools OP + full-sib; paternal is full-sib sires only).
  inm <- paste0(order[-z], "-", order[-1])
  moms  <- unique(c(vapply(fs_ids, function(cid) cx[[cid]]$mother_id, character(1)), dams_op))
  sires <- unique(vapply(fs_ids, function(cid) cx[[cid]]$father_id, character(1)))
  final_ll <- 0; final_m_tot <- numeric(Ti); final_p_tot <- numeric(Ti)
  mat_by_mom <- matrix(0, Ti, length(moms), dimnames = list(inm, moms))
  pat_by_sire <- if (length(sires)) matrix(0, Ti, length(sires), dimnames = list(inm, sires)) else NULL
  for (cid in fs_ids) { b <- fs_built[[cid]]
    es <- fs_estep_cpp(b$G, b$Am, b$Ap, r_m, r_p, epsilon, FALSE)
    final_ll <- final_ll + es$loglik; final_m_tot <- final_m_tot + es$total; final_p_tot <- final_p_tot + es$total
    mat_by_mom[, cx[[cid]]$mother_id] <- mat_by_mom[, cx[[cid]]$mother_id] + es$total
    pat_by_sire[, cx[[cid]]$father_id] <- pat_by_sire[, cx[[cid]]$father_id] + es$total }
  for (cid in op_ids) { b <- op_built[[cid]]
    es <- op_estep_cpp(b$G, b$Am, r_m, q_list[[b$mother]], epsilon)
    final_ll <- final_ll + es$loglik; final_m_tot <- final_m_tot + es$total
    mat_by_mom[, b$mother] <- mat_by_mom[, b$mother] + es$total }
  final_pen <- sum(vapply(q_list, q_pen, numeric(1)))
  final_obj <- final_ll + final_pen
  if (.fs_obj_decreased(final_obj, prev_obj)) obj_dec <- TRUE   # decrease at the final M-step
  obj_trace <- c(obj_trace, final_obj); ll_trace <- c(ll_trace, final_ll)
  dm_trace <- c(dm_trace, 0); dp_trace <- c(dp_trace, 0); dq_trace <- c(dq_trace, 0)
  if (obj_dec) { converged <- FALSE; conv_reason <- "objective_decreased" }
  if (!converged) {
    if (identical(conv_reason, "objective_decreased"))
      warning("hmm_map_mixed(): the active EM objective decreased materially; ",
              "results may be unreliable.", call. = FALSE)
    else
      warning("hmm_map_mixed(): EM did not converge in ", it, " iterations; the ",
              "returned estimates are the last iterate. Increase `maxit` or relax `tol`.",
              call. = FALSE)
  }

  post <- NULL
  if (isTRUE(return_posterior)) {
    post <- lapply(fs_ids, function(cid) { b <- fs_built[[cid]]
      fs_estep_cpp(b$G, b$Am, b$Ap, r_m, r_p, epsilon, TRUE)$gamma }); names(post) <- fs_ids
  }
  info <- .fs_informativeness(c(fs_built, op_built), fs_built, Ti)
  both_inf <- info$mat_inf & info$pat_inf
  rep_m <- .fs_apply_informativeness(.fs_interval_report(r_m, final_m_tot, gap_r), info$mat_inf, info$globally_exchangeable, both_inf)
  rep_p <- .fs_apply_informativeness(.fs_interval_report(r_p, final_p_tot, gap_r), info$pat_inf, info$globally_exchangeable, both_inf)
  ident_tab <- data.frame(
    interval = inm, n_offspring = info$n_off,
    maternal_informative = info$mat_inf, paternal_informative = info$pat_inf,
    exchangeable = info$globally_exchangeable & both_inf,
    status_m = rep_m$status, status_p = rep_p$status, stringsAsFactors = FALSE)
  fit <- list(
    r_m = stats::setNames(r_m, inm),
    r_p = stats::setNames(r_p, inm),
    d_m = stats::setNames(rep_m$dist, inm), d_p = stats::setNames(rep_p$dist, inm),
    interval_status_m = stats::setNames(rep_m$status, inm),
    interval_status_p = stats::setNames(rep_p$status, inm),
    identifiability = ident_tab, identifiable_labels = !info$globally_exchangeable,
    gap_r = gap_r,
    meiosis_count_m = stats::setNames(final_m_tot, inm),
    meiosis_count_p = stats::setNames(final_p_tot, inm),
    # sex-specific CONSENSUS maps: r_m pools maternal meioses across ALL crosses (OP +
    # full-sib), r_p pools paternal meioses across genotyped sires. NOT parent-specific.
    map_scope = "sex_specific_consensus",
    maternal_meioses_by_mother = mat_by_mom,
    paternal_meioses_by_sire = pat_by_sire,
    logLik = final_ll, q_penalty = final_pen, objective = final_obj,
    converged = converged, conv_reason = conv_reason, iters = it,
    objective_decreased = obj_dec,
    objective_trace = obj_trace, loglik_trace = ll_trace,
    max_dr_m_trace = dm_trace, max_dr_p_trace = dp_trace, max_dq_trace = dq_trace,
    epsilon = epsilon, q = q_list, posterior = post)
  out <- list(order = order, fit = fit, map_scope = "sex_specific_consensus",
              contributing_crosses = c(fs_ids, op_ids),
              contributing_mothers = moms, contributing_sires = sires,
              maternal_crosses = c(fs_ids, op_ids), paternal_crosses = fs_ids,
              family_type = stats::setNames(ftype[c(fs_ids, op_ids)], c(fs_ids, op_ids)),
              parent_phase = list(maternal = phased_m[moms], paternal = phased_p[sires]),
              dispatched = FALSE)
  class(out) <- "HSMap.mixed"
  out
}
