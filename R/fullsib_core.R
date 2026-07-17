# Known-sire full-sib mapping driver (oracle phase). See dev/known_sire_design.md.

# Build a cross's maternal/paternal 2xT allele matrices from parent genotypes and
# oracle phase vectors, aligned to `order`. Errors on a missing required parent
# genotype at a fitted marker (on_missing = "error").
.fs_build_cross <- function(cr, parent_geno, phased_m, phased_p, order, on_missing) {
  mom <- cr$mother_id; fat <- cr$father_id
  gm <- parent_geno[[mom]]; gf <- parent_geno[[fat]]
  if (is.null(gm) || all(is.na(gm)))
    stop("hmm_map_fullsib(): mother '", mom, "' has no genotype (cross '", cr$cross_id, "').", call. = FALSE)
  if (is.null(gf) || all(is.na(gf)))
    stop("hmm_map_fullsib(): sire '", fat, "' has no genotype (cross '", cr$cross_id,
         "'); full-sib fitting requires a genotyped sire.", call. = FALSE)
  names(gm) <- names(gm) %||% colnames(cr$G); names(gf) <- names(gf) %||% colnames(cr$G)
  gm <- gm[order]; gf <- gf[order]
  pvm <- phased_m[[mom]]; pvf <- phased_p[[fat]]
  if (is.null(pvm)) stop("hmm_map_fullsib(): no maternal phase supplied for mother '", mom, "'.", call. = FALSE)
  if (is.null(pvf)) stop("hmm_map_fullsib(): no paternal phase supplied for sire '", fat, "'.", call. = FALSE)
  if (length(pvm) != length(order) - 1L || length(pvf) != length(order) - 1L)
    stop("hmm_map_fullsib(): phase vectors for cross '", cr$cross_id,
         "' must have length length(order) - 1.", call. = FALSE)
  Am <- phase_to_haplotypes(gm, pvm)
  Ap <- phase_to_haplotypes(gf, pvf)
  miss_m <- which(apply(Am, 2, anyNA)); miss_p <- which(apply(Ap, 2, anyNA))
  if (length(miss_m) || length(miss_p)) {
    if (identical(on_missing, "error"))
      stop("hmm_map_fullsib(): missing required parent genotype in cross '", cr$cross_id,
           "' at ", length(unique(c(miss_m, miss_p))), " marker(s) (e.g. ",
           paste(order[utils::head(unique(c(miss_m, miss_p)), 5)], collapse = ", "),
           "). The core does not invent a parental allele; supply complete parent ",
           "genotypes or drop these markers first.", call. = FALSE)
  }
  G <- cr$G[, order, drop = FALSE]; storage.mode(G) <- "integer"
  list(G = G, Am = matrix(as.integer(Am), 2), Ap = matrix(as.integer(Ap), 2))
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
#' @param phased_m Named list mapping each mother id to an oracle maternal phase vector
#'   (length \code{z-1}, 1 = coupling, 0 = repulsion), aligned to \code{order}.
#' @param phased_p Named list mapping each sire id to an oracle paternal phase vector.
#' @param crosses Optional character vector of \code{cross_id}s to fit; default all
#'   \code{known_sire_genotyped} crosses.
#' @param order Optional marker order (character); default the common column order of
#'   the selected crosses' genotype matrices.
#' @param epsilon Symmetric genotyping error rate.
#' @param tol Convergence tolerance on the relative active objective and on \code{r}.
#' @param maxit Maximum EM iterations.
#' @param r_start Initial recombination fraction for both maps.
#' @param on_missing One of \code{"error"} (default): stop on a missing required parent
#'   genotype at a fitted marker.
#' @param return_posterior Logical; if \code{TRUE}, include per-cross posterior
#'   inheritance probabilities.
#'
#' @return An object of class \code{HSMap.fullsib} (see \code{$fit}).
#' @export
hmm_map_fullsib <- function(x, phased_m, phased_p, crosses = NULL, order = NULL,
                            epsilon = 0.05, tol = 1e-6, maxit = 1000L, r_start = 0.1,
                            on_missing = c("error"), return_posterior = FALSE) {
  if (!inherits(x, "HSMap.data")) stop("`x` must be an HSMap.data object.")
  if (is.null(x$crosses)) stop("`x` lacks cross-aware fields; re-read with read_HSMap_data().")
  on_missing <- match.arg(on_missing)
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
    .fs_build_cross(all_cross[[cid]], x$parent_genotypes, phased_m, phased_p, order, on_missing))
  names(built) <- fs_ids
  moms <- unique(vapply(fs_ids, function(cid) all_cross[[cid]]$mother_id, character(1)))
  sires <- unique(vapply(fs_ids, function(cid) all_cross[[cid]]$father_id, character(1)))

  r_m <- rep(r_start, Ti); r_p <- rep(r_start, Ti)
  r_lo <- 0.0; r_hi <- 0.5 - 1e-12
  clampr <- function(r) pmin(pmax(r, r_lo), r_hi)

  ll_trace <- numeric(0); prev_ll <- -Inf; obj_dec <- FALSE
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
    ll_trace <- c(ll_trace, ll)
    if (it > 1L && ll < prev_ll - OBJ_DEC_REL * (1 + abs(prev_ll))) obj_dec <- TRUE
    r_m_new <- clampr(ifelse(tot > 0, m_sw / tot, r_m))
    r_p_new <- clampr(ifelse(tot > 0, p_sw / tot, r_p))
    dr <- max(max(abs(r_m_new - r_m)), max(abs(r_p_new - r_p)))
    rel <- if (is.finite(prev_ll)) abs(ll - prev_ll) / (1 + abs(prev_ll)) else Inf
    r_m <- r_m_new; r_p <- r_p_new
    if (it > 1L && rel < tol && dr < tol) { converged <- TRUE; conv_reason <- "relative_objective_and_params_stable"; prev_ll <- ll; break }
    prev_ll <- ll
  }
  # final log-likelihood evaluated at the final parameters
  final_ll <- 0
  for (cid in fs_ids) { b <- built[[cid]]; final_ll <- final_ll + fs_loglik_cpp(b$G, b$Am, b$Ap, r_m, r_p, epsilon) }
  if (obj_dec) { converged <- FALSE; conv_reason <- "objective_decreased" }

  post <- NULL
  if (isTRUE(return_posterior)) {
    post <- lapply(fs_ids, function(cid) { b <- built[[cid]]
      es <- fs_estep_cpp(b$G, b$Am, b$Ap, r_m, r_p, epsilon, TRUE); es$gamma })
    names(post) <- fs_ids
  }

  fit <- list(
    r_m = stats::setNames(r_m, paste0(order[-z], "-", order[-1])),
    r_p = stats::setNames(r_p, paste0(order[-z], "-", order[-1])),
    d_m = haldane(r_m), d_p = haldane(r_p),
    logLik = final_ll, objective = final_ll,
    converged = converged, conv_reason = conv_reason, iters = it,
    objective_decreased = obj_dec, loglik_trace = ll_trace,
    epsilon = epsilon, posterior = post
  )
  out <- list(
    order = order, fit = fit,
    contributing_crosses = fs_ids,
    contributing_mothers = moms, contributing_sires = sires,
    family_type = stats::setNames(rep("known_sire_genotyped", length(fs_ids)), fs_ids),
    parent_phase = list(maternal = phased_m[moms], paternal = phased_p[sires])
  )
  class(out) <- "HSMap.fullsib"
  out
}
