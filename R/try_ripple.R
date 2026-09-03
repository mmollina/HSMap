# Local order diagnostics: placement support for one marker (try_marker) and
# likelihood support for a local order (ripple_map). Both are DIAGNOSTIC: they
# score candidate orders with the HSMap multipoint likelihood and report what
# they find. Neither searches globally, and neither modifies a map unless the
# user explicitly asks for it.
#
# Inspired by the TRY and RIPPLE commands of MAPMAKER (Lander et al. 1987),
# reimplemented for the HSMap maternal half-sib model.

`%|NA|%` <- function(a, b) if (length(a) && !is.na(a)) a else b

## ---------------------------------------------------------------------------
## Phase is a MARKER property, not an interval property.
##
## The stored phase_vec describes intervals of ONE order, so it cannot be
## carried over to a permuted order. We recover the underlying marker-level
## state h(marker) = "does homolog 1 carry the reference allele", via
##   h[k+1] = h[k] XOR (phase_vec[k] == 0),
## and then regenerate the interval phase of ANY order as
##   phase[i] = (h[m_i] == h[m_{i+1}]).
## Reconstructing the original order returns the original phase_vec exactly.
## ---------------------------------------------------------------------------
.hsmap_marker_phase <- function(order, phase_vec) {
  order <- as.character(order)
  Tn <- length(order)
  if (length(phase_vec) != Tn - 1L)
    stop("`phase_vec` must have length(order) - 1.", call. = FALSE)
  if (anyNA(phase_vec))
    stop("Local order diagnostics require a fully resolved phase; ",
         sum(is.na(phase_vec)), " interval(s) are NA. Restrict `order` to a ",
         "resolved block first.", call. = FALSE)
  h <- integer(Tn)
  h[1L] <- 1L
  for (k in seq_len(Tn - 1L))
    h[k + 1L] <- if (phase_vec[k] == 1L) h[k] else 1L - h[k]
  stats::setNames(h, order)
}

.hsmap_phase_for <- function(mk, h) {
  hv <- h[as.character(mk)]
  if (anyNA(hv))
    stop("No marker-phase state for: ",
         paste(utils::head(mk[is.na(hv)], 3), collapse = ", "), call. = FALSE)
  as.integer(hv[-length(hv)] == hv[-1L])
}

## Gap-aware local length. Intervals at/above `no_linkage_r` are no-linkage
## boundaries and are NEVER converted into large centimorgan values; they are
## excluded from the length and counted separately (same convention as
## diagnose_map_intervals()).
.hsmap_local_len <- function(r, no_linkage_r = 0.499) {
  gap <- r >= no_linkage_r
  list(cM = sum(-50 * log(1 - 2 * pmin(r[!gap], 0.49999))), n_gap = sum(gap))
}

## Resolve a dam index or name to the dam NAME. hmm_map() matches the phased
## object to the data by name, so an index must never reach the phased object.
.hsmap_dam_id <- function(x, dam) {
  nm <- names(x$G_list)
  id <- if (is.numeric(dam)) nm[as.integer(dam)] else as.character(dam)
  if (length(id) != 1L || is.na(id) || !id %in% nm)
    stop("Unknown dam: ", paste(dam, collapse = ", "), call. = FALSE)
  id
}

## One local fit through the PRODUCTION engine. `dam` must be a dam NAME.
.hsmap_fit_local <- function(x, mk, h, dam, epsilon, lambda, tol, maxit,
                             r_start, threads = NULL) {
  ph <- structure(list(dam = dam, order = as.character(mk),
                       clusters = rep(1L, length(mk)),
                       phase_vec = .hsmap_phase_for(mk, h)),
                  class = "HSMap.phased")
  m <- suppressWarnings(hmm_map(x, phased = ph, dam = dam, epsilon = epsilon,
                                lambda = lambda, r_start = r_start,
                                threads = threads, tol = tol, maxit = maxit))
  out <- list(logLik = as.numeric(m$fit[["logLik"]]),
              r = as.numeric(m$fit$r),
              converged = isTRUE(m$fit$converged))
  rm(m)
  out
}

## Score a list of candidate orders.
##
## Two-stage by design. m21 established that the EM surface can be multimodal,
## but that the multimodality is confined to badly wrong orders: among the
## top-ranked candidates the spread across starts was ~1e-4 log units, while
## non-contenders reached tens of log units. So every candidate is screened at
## a single start and only the `refine_top` best are re-fitted from all
## `r_starts`, each scored by its BEST converged likelihood. Set
## refine_top = Inf to multi-start every candidate.
.hsmap_score_orders <- function(x, cands, h, dam, epsilon, lambda, tol, maxit,
                                r_starts, refine_top, no_linkage_r,
                                threads = NULL) {
  n <- length(cands)
  screen_start <- r_starts[length(r_starts)]
  ll <- rep(NA_real_, n)
  fits <- vector("list", n)
  for (i in seq_len(n)) {
    f <- .hsmap_fit_local(x, cands[[i]], h, dam, epsilon, lambda, tol, maxit,
                          screen_start, threads)
    fits[[i]] <- f
    if (f$converged) ll[i] <- f$logLik
  }
  best_r <- rep(screen_start, n)
  spread <- rep(0, n)
  k <- if (is.finite(refine_top)) min(as.integer(refine_top), n) else n
  refine <- if (k >= n) seq_len(n) else utils::head(order(-ll, na.last = TRUE), k)
  for (i in refine) {
    lls <- rep(NA_real_, length(r_starts))
    best <- fits[[i]]
    if (best$converged) lls[length(r_starts)] <- best$logLik else best <- NULL
    for (j in seq_along(r_starts)) {
      if (r_starts[j] == screen_start) next
      f <- .hsmap_fit_local(x, cands[[i]], h, dam, epsilon, lambda, tol, maxit,
                            r_starts[j], threads)
      if (!f$converged) next
      lls[j] <- f$logLik
      if (is.null(best) || f$logLik > best$logLik) {
        best <- f; best_r[i] <- r_starts[j]
      }
    }
    if (!is.null(best)) {
      fits[[i]] <- best
      ll[i] <- best$logLik
      spread[i] <- if (sum(!is.na(lls)) > 1L) diff(range(lls, na.rm = TRUE)) else 0
    }
  }
  len <- lapply(fits, function(f) .hsmap_local_len(f$r, no_linkage_r))
  list(logLik = ll,
       converged = vapply(fits, function(f) f$converged, TRUE),
       r = lapply(fits, function(f) f$r),
       cM = vapply(len, function(z) z$cM, 0),
       n_no_linkage = vapply(len, function(z) z$n_gap, 0L),
       r_start = best_r, ll_spread = spread,
       refined = seq_len(n) %in% refine)
}

## Resolve `phased`/`map` input to order + phase_vec.
.hsmap_order_phase <- function(object) {
  if (inherits(object, "HSMap.map"))
    return(list(order = as.character(object$order),
                phase_vec = as.integer(object$phase_vec),
                dam = object$dam))
  if (inherits(object, "HSMap.phased"))
    return(list(order = as.character(object$order),
                phase_vec = as.integer(object$phase_vec),
                dam = object$dam))
  stop("`phased` must be an HSMap.phased or HSMap.map object.", call. = FALSE)
}

#' Placement support for a single marker against a fixed local order
#'
#' @description
#' Evaluates where a marker is best supported by the HSMap multipoint
#' likelihood, given an established local order. The target marker is removed,
#' the remaining scaffold is held fixed, and the marker is reinserted into
#' every gap; each placement is refitted from scratch and scored. The
#' \strong{complete placement profile} is returned, not just the winner: a flat
#' or multi-peaked profile is itself the scientific result, and often the
#' correct answer is that the position is not resolved.
#'
#' This is a local diagnostic. It does not search for a globally optimal map,
#' and it never modifies anything.
#'
#' @details
#' Every placement is a genuine refit: the recombination fractions are
#' re-estimated for that order, so distances are never inherited from the
#' baseline order. All placements use the same markers, the same offspring and
#' the same number of parameters, so their likelihoods are directly comparable.
#'
#' \strong{Phase.} Relative phase is a property of markers, not of intervals,
#' so it is reconstructed at marker level from the supplied phase and
#' regenerated for each candidate order (see \code{\link{ripple_map}}). If the
#' target marker is absent from the supplied order its orientation is unknown;
#' both orientations are then scored and the better one reported
#' (\code{orient = "auto"}).
#'
#' \strong{Unlinked markers.} A marker unlinked to the scaffold has no internal
#' placement to find. Its profile is strongly edge-favouring with both edges
#' tied, because inserting it internally severs the Markov chain (two
#' no-linkage boundaries) whereas an edge placement leaves the chromosome
#' intact. This case is reported through \code{placement}, not disguised as a
#' winning position.
#'
#' \strong{Edge placements carry one fewer interval} than interior ones, so
#' compare edge and interior likelihoods with that in mind (\code{at_edge}
#' flags them).
#'
#' \strong{Ambiguity.} Placements within \code{tie_tol} log-likelihood units of
#' the best are reported as near-optimal. The default of 2 units is
#' deliberately conservative; it is a reporting convention, not a test. Treat a
#' wide near-optimal set as "unresolved at this information content" rather
#' than as evidence for the single best position.
#'
#' @param x An \code{HSMap.data} object.
#' @param phased An \code{HSMap.phased} or \code{HSMap.map} giving the
#'   established order and its resolved phase.
#' @param marker Character scalar: the marker whose placement is tested. It may
#'   be absent from \code{phased} (a new marker placed against a framework).
#' @param framework Optional character vector of scaffold markers to test
#'   against, in order. If \code{NULL} (default) a window of \code{window}
#'   markers centred on the marker's current position is used.
#' @param window Scaffold size when \code{framework} is \code{NULL}. Default
#'   \code{21}.
#' @param dam Dam index or name. Default \code{1}.
#' @param epsilon,lambda,tol,maxit Passed to \code{\link{hmm_map}}. The
#'   defaults match \code{hmm_map()}; pass the same values used to fit the map
#'   under test.
#' @param r_starts Initial recombination fractions for the multi-start scoring.
#'   Default \code{c(0.001, 0.005, 0.01, 0.05)}.
#' @param refine_top Number of best screened placements re-fitted from all
#'   \code{r_starts}. Default \code{Inf} (all).
#' @param orient One of \code{"auto"} (default; fixed when the marker is in
#'   \code{phased}, both orientations otherwise), \code{"fixed"} or
#'   \code{"both"}.
#' @param tie_tol Log-likelihood tolerance defining the near-optimal placement
#'   set. Default \code{2}.
#' @param no_linkage_r Recombination fraction at/above which an interval is a
#'   no-linkage boundary, excluded from local length. Default \code{0.499}.
#' @param unlinked_r If, at its best placement, the target's recombination
#'   fraction to every neighbour it has is at least this value, the marker is
#'   reported as having no credible internal placement. Default \code{0.4}.
#' @param threads Worker threads for the HMM, passed to \code{\link{hmm_map}}.
#'   \code{NULL} (default) leaves the engine's current setting. These functions
#'   issue thousands of small fits, so on a busy machine set this to a modest
#'   number (2-4) rather than letting every fit fan out over all cores.
#'
#' @return An object of class \code{"HSMap.try"}: a list with \code{marker},
#'   \code{framework}, \code{current_pos} (\code{NA} if the marker was absent),
#'   \code{best_pos}, \code{displacement}, \code{profile} (one row per
#'   insertion: \code{position}, \code{logLik}, \code{dLL_best},
#'   \code{dLL_current}, \code{local_cM}, \code{n_no_linkage},
#'   \code{converged}, \code{r_start}, \code{ll_spread}, \code{at_edge},
#'   \code{orientation}), \code{orders} (the candidate order for each row),
#'   \code{near_optimal}, \code{unique_optimum}, \code{best_is_edge},
#'   \code{min_flank_r_at_best},
#'   \code{any_no_linkage}, \code{placement} (one of \code{"resolved"},
#'   \code{"unresolved (statistically equivalent positions)"} or
#'   \code{"no credible internal placement (unlinked)"}) and \code{settings}.
#'
#'
#' @section Limitations:
#' \itemize{
#'   \item \strong{Local, not global.} Both functions score a fixed set of
#'     local candidates. They cannot find an error that lies outside the window
#'     or framework supplied, and neither returns a globally optimal map.
#'   \item \strong{Finite mapping resolution.} Order is only identifiable to
#'     the extent the data support it. With few offspring, tight spacing or
#'     heavy missingness, many orders are statistically equivalent and the
#'     honest answer is a tie, not a winner.
#'   \item \strong{Ties and co-segregation.} Markers separated by no observed
#'     recombination have no identifiable relative order; the near-optimal set
#'     will legitimately contain several positions or permutations.
#'   \item \strong{Model misspecification.} Likelihoods are computed under the
#'     HSMap maternal model (paternal alleles independent across markers,
#'     a single global emission-error rate). Where those assumptions fail, the
#'     ranking can be confidently wrong.
#'   \item \strong{Genotype error.} A noisy marker can produce a displaced or
#'     irregular profile. An abnormal profile is evidence that something is
#'     wrong at that marker, not proof of either a wrong position or a bad
#'     genotype; the two are not reliably separable from the profile alone.
#'   \item \strong{Likelihood is the criterion, not map length.} A
#'     better-supported order is often \emph{longer}. Do not use these tools to
#'     shorten a map.
#' }
#'
#' @examples
#' \dontrun{
#' tr <- try_marker(dat, map, marker = "SNP123", window = 21)
#' tr
#' plot(tr)
#' subset(tr$profile, dLL_best > -2)
#' }
#' @seealso \code{\link{ripple_map}} for permuting a local order,
#'   \code{\link{diagnose_map_intervals}} for interval-level diagnostics.
#' @export
try_marker <- function(x, phased, marker, framework = NULL, window = 21L,
                       dam = 1, epsilon = 0.01, lambda = 20, tol = 1e-6,
                       maxit = 1000L,
                       r_starts = c(0.001, 0.005, 0.01, 0.05),
                       refine_top = Inf, orient = c("auto", "fixed", "both"),
                       tie_tol = 2, no_linkage_r = 0.499, unlinked_r = 0.4,
                       threads = NULL) {
  if (!inherits(x, "HSMap.data")) stop("`x` must be an HSMap.data.", call. = FALSE)
  if (!is.character(marker) || length(marker) != 1L)
    stop("`marker` must be a single marker ID.", call. = FALSE)
  orient <- match.arg(orient)
  if (!is.numeric(tie_tol) || length(tie_tol) != 1L || tie_tol < 0)
    stop("`tie_tol` must be a single non-negative number.", call. = FALSE)
  if (!length(r_starts) || any(r_starts <= 0 | r_starts >= 0.5))
    stop("`r_starts` must be in (0, 0.5).", call. = FALSE)

  op <- .hsmap_order_phase(phased)
  dam_id <- .hsmap_dam_id(x, dam)
  if (!marker %in% colnames(x$G_list[[dam_id]]))
    stop("`marker` is not present in the genotype data.", call. = FALSE)

  h_all <- .hsmap_marker_phase(op$order, op$phase_vec)
  cur_global <- match(marker, op$order)

  ## scaffold: supplied framework, or a window centred on the marker
  if (!is.null(framework)) {
    scaf <- setdiff(as.character(framework), marker)
    if (length(scaf) < 2L)
      stop("`framework` needs at least 2 markers besides `marker`.", call. = FALSE)
    if (!all(scaf %in% op$order))
      stop("All `framework` markers must be in the supplied order.", call. = FALSE)
    scaf <- op$order[sort(match(scaf, op$order))]
  } else {
    if (is.na(cur_global))
      stop("`marker` is not in the supplied order; supply `framework` ",
           "explicitly.", call. = FALSE)
    window <- as.integer(window)
    if (window < 3L) stop("`window` must be at least 3.", call. = FALSE)
    half <- (window - 1L) %/% 2L
    lo <- max(1L, cur_global - half)
    hi <- min(length(op$order), lo + window - 1L)
    lo <- max(1L, hi - window + 1L)
    scaf <- setdiff(op$order[lo:hi], marker)
  }

  ## the marker's own phase state: known if it is already in the order,
  ## otherwise both orientations are scored
  use_both <- switch(orient, auto = is.na(cur_global), fixed = FALSE, both = TRUE)
  if (!use_both && is.na(cur_global))
    stop("orient = \"fixed\" needs the marker to be in the supplied order.",
         call. = FALSE)

  npos <- length(scaf) + 1L
  cands <- lapply(seq_len(npos), function(p) append(scaf, marker, after = p - 1L))
  ## position of the marker's CURRENT placement within this scaffold
  cur_pos <- if (!is.na(cur_global)) {
    left <- sum(match(scaf, op$order) < cur_global); left + 1L
  } else NA_integer_

  run <- function(hval) {
    h <- h_all[c(scaf, marker)]
    names(h) <- c(scaf, marker)
    if (!is.null(hval)) h[marker] <- hval
    .hsmap_score_orders(x, cands, h, dam_id, epsilon, lambda, tol, maxit,
                        r_starts, refine_top, no_linkage_r, threads)
  }
  if (use_both) {
    s0 <- run(0L); s1 <- run(1L)
    take1 <- !is.na(s1$logLik) & (is.na(s0$logLik) | s1$logLik > s0$logLik)
    sc <- s0
    for (nm in c("logLik", "converged", "cM", "n_no_linkage", "r_start",
                 "ll_spread", "refined"))
      sc[[nm]][take1] <- s1[[nm]][take1]
    sc$r[take1] <- s1$r[take1]
    orientation <- ifelse(take1, 1L, 0L)
  } else {
    sc <- run(NULL)
    orientation <- rep(unname(h_all[marker]), npos)
  }

  if (all(is.na(sc$logLik)))
    stop("No placement converged; check `epsilon`/`lambda`/`maxit`.", call. = FALSE)
  best <- which.max(sc$logLik)
  ## The target's own linkage: r to its left and right neighbour at each
  ## placement (one side only at an edge). This is the direct evidence of
  ## whether the marker is linked to the scaffold at all.
  r_left <- vapply(seq_len(npos), function(p)
    if (p > 1L) sc$r[[p]][p - 1L] else NA_real_, 0)
  r_right <- vapply(seq_len(npos), function(p)
    if (p < npos) sc$r[[p]][p] else NA_real_, 0)
  min_flank <- pmin(r_left, r_right, na.rm = TRUE)
  prof <- data.frame(
    position = seq_len(npos),
    logLik = sc$logLik,
    dLL_best = sc$logLik - sc$logLik[best],
    dLL_current = if (is.na(cur_pos)) NA_real_ else sc$logLik - sc$logLik[cur_pos],
    local_cM = sc$cM,
    r_left = r_left, r_right = r_right, min_flank_r = min_flank,
    n_no_linkage = sc$n_no_linkage,
    converged = sc$converged,
    r_start = sc$r_start,
    ll_spread = sc$ll_spread,
    at_edge = seq_len(npos) %in% c(1L, npos),
    orientation = orientation,
    stringsAsFactors = FALSE)
  near <- which(prof$dLL_best >= -tie_tol)

  ## A plain verdict, because "best position" is misleading when the marker
  ## cannot be placed at all. The decisive evidence is the target's OWN
  ## linkage: if at its best placement it is unlinked to every neighbour it
  ## has, there is nothing to place. (Judging this by no-linkage interval
  ## counts alone is unreliable -- the EM commonly settles near r = 0.49,
  ## just under `no_linkage_r`, which is unlinked in substance but not by that
  ## threshold.) Such a marker is also strongly edge-favouring: an internal
  ## insertion severs the Markov chain, whereas an edge placement leaves the
  ## chromosome intact, so both edges tie at the top.
  unlinked <- !is.na(min_flank[best]) && min_flank[best] >= unlinked_r
  best_edge <- prof$at_edge[best]
  placement <- if (unlinked)
      "no credible internal placement (unlinked)"
    else if (length(near) == 1L) "resolved"
    else "unresolved (statistically equivalent positions)"

  structure(list(
    marker = marker, framework = scaf,
    current_pos = cur_pos, best_pos = best,
    displacement = if (is.na(cur_pos)) NA_integer_ else best - cur_pos,
    profile = prof, orders = cands,
    near_optimal = near, unique_optimum = length(near) == 1L,
    placement = placement, best_is_edge = best_edge,
    min_flank_r_at_best = min_flank[best],
    any_no_linkage = any(prof$n_no_linkage > 0L),
    settings = list(window = length(scaf) + 1L, dam = dam, epsilon = epsilon,
                    lambda = lambda, r_starts = r_starts, tie_tol = tie_tol,
                    refine_top = refine_top, orient = orient,
                    no_linkage_r = no_linkage_r)),
    class = "HSMap.try")
}

#' Likelihood support for a local marker order
#'
#' @description
#' Slides a window along an established order and, for every window position,
#' evaluates all \code{factorial(window)} permutations of the window's markers
#' with the markers outside it held fixed. Each candidate order is refitted from
#' scratch and scored by the HSMap multipoint likelihood, so orders are compared
#' honestly rather than through recombination fractions that are only valid for
#' the baseline.
#'
#' By default it \strong{reports} candidate improvements and changes nothing.
#'
#' @details
#' \strong{Likelihood decides, length does not.} A better-supported order may
#' well be \emph{longer} than the baseline; local map length is reported as a
#' diagnostic only and is never a selection criterion.
#'
#' \strong{Flanking context.} Each window is fitted inside a local sub-map with
#' \code{anchor} fixed markers on each side rather than on the whole
#' chromosome. This matters: information propagates from the fixed flanks into
#' the permuted block, and a narrow frame leaves the block under-determined and
#' can manufacture apparent order conflicts. The default of 10 is deliberately
#' generous; \code{anchor = Inf} uses the entire order (exact, and much
#' slower).
#'
#' \strong{Phase.} Relative phase is a marker property. The interval phase of
#' the supplied order is converted to marker-level states and regenerated for
#' each permutation; restoring the baseline order restores the baseline phase
#' exactly. Permutations therefore never reinterpret the old order's interval
#' phase.
#'
#' \strong{Scope.} This validates local orders one window at a time. It is not
#' an order-search engine: it does not walk the map accepting changes and
#' recomputing. The experimental \code{apply} mode performs a single pass over
#' non-overlapping accepted windows and is \code{FALSE} by default.
#'
#' @param x An \code{HSMap.data} object.
#' @param phased An \code{HSMap.phased} or \code{HSMap.map}.
#' @param window Number of markers permuted per window. Default \code{5}
#'   (120 permutations). Values above 6 are expensive.
#' @param by Step between successive window starts. Default \code{1}.
#' @param positions Optional integer vector of window start positions to
#'   evaluate; overrides \code{by}.
#' @param anchor Fixed flanking markers on each side of the window. Default
#'   \code{10}; use \code{Inf} for the whole order.
#' @param dam,epsilon,lambda,tol,maxit,r_starts,no_linkage_r,threads See
#'   \code{\link{try_marker}}.
#' @param refine_top Number of best screened permutations per window re-fitted
#'   from all \code{r_starts}. Default \code{10}.
#' @param tie_tol Log-likelihood tolerance defining near-optimal orders.
#'   Default \code{2}.
#' @param min_delta Log-likelihood gain required to call a window a supported
#'   improvement. Default \code{2}.
#' @param apply Experimental. If \code{TRUE}, returns a revised order applying
#'   accepted non-overlapping windows in a single pass. Default \code{FALSE}.
#' @param verbose Print per-window progress. Default \code{FALSE}.
#'
#' @return An object of class \code{"HSMap.ripple"}: a list with
#'   \code{windows} (one row per window: \code{window}, \code{start},
#'   \code{end}, \code{baseline_logLik}, \code{best_logLik}, \code{dLL},
#'   \code{identity_best}, \code{n_better}, \code{n_near_optimal},
#'   \code{baseline_cM}, \code{best_cM}, \code{dcM}, \code{converged},
#'   \code{r_start}, \code{ll_spread}, \code{supported}), \code{best_orders},
#'   \code{alternatives} (top candidates per window), \code{order} (baseline),
#'   \code{applied_order} (only when \code{apply = TRUE}) and \code{settings}.
#'
#'
#' @section Limitations:
#' \itemize{
#'   \item \strong{Local, not global.} Both functions score a fixed set of
#'     local candidates. They cannot find an error that lies outside the window
#'     or framework supplied, and neither returns a globally optimal map.
#'   \item \strong{Finite mapping resolution.} Order is only identifiable to
#'     the extent the data support it. With few offspring, tight spacing or
#'     heavy missingness, many orders are statistically equivalent and the
#'     honest answer is a tie, not a winner.
#'   \item \strong{Ties and co-segregation.} Markers separated by no observed
#'     recombination have no identifiable relative order; the near-optimal set
#'     will legitimately contain several positions or permutations.
#'   \item \strong{Model misspecification.} Likelihoods are computed under the
#'     HSMap maternal model (paternal alleles independent across markers,
#'     a single global emission-error rate). Where those assumptions fail, the
#'     ranking can be confidently wrong.
#'   \item \strong{Genotype error.} A noisy marker can produce a displaced or
#'     irregular profile. An abnormal profile is evidence that something is
#'     wrong at that marker, not proof of either a wrong position or a bad
#'     genotype; the two are not reliably separable from the profile alone.
#'   \item \strong{Likelihood is the criterion, not map length.} A
#'     better-supported order is often \emph{longer}. Do not use these tools to
#'     shorten a map.
#' }
#'
#' @examples
#' \dontrun{
#' rp <- ripple_map(dat, map, window = 5, anchor = 10)
#' rp
#' subset(rp$windows, supported)
#' }
#' @seealso \code{\link{try_marker}}
#' @export
ripple_map <- function(x, phased, window = 5L, by = 1L, positions = NULL,
                       anchor = 10L, dam = 1, epsilon = 0.01, lambda = 20,
                       tol = 1e-6, maxit = 1000L,
                       r_starts = c(0.001, 0.005, 0.01, 0.05),
                       refine_top = 10L, tie_tol = 2, min_delta = 2,
                       no_linkage_r = 0.499, apply = FALSE, threads = NULL,
                       verbose = FALSE) {
  if (!inherits(x, "HSMap.data")) stop("`x` must be an HSMap.data.", call. = FALSE)
  window <- as.integer(window)
  if (window < 2L) stop("`window` must be at least 2.", call. = FALSE)
  if (window > 7L)
    stop("`window` > 7 means ", format(factorial(8)), "+ permutations per ",
         "window; refusing. Use a smaller window.", call. = FALSE)
  dam_id <- .hsmap_dam_id(x, dam)
  op <- .hsmap_order_phase(phased)
  ord <- op$order
  Tn <- length(ord)
  if (Tn < window + 1L)
    stop("`order` is shorter than the window.", call. = FALSE)
  h_all <- .hsmap_marker_phase(ord, op$phase_vec)

  starts <- if (!is.null(positions)) as.integer(positions)
            else seq(1L, Tn - window + 1L, by = as.integer(by))
  if (any(starts < 1L | starts > Tn - window + 1L))
    stop("`positions` out of range.", call. = FALSE)

  perms <- .hsmap_perms(window)                      # identity is row 1

  rows <- vector("list", length(starts))
  best_orders <- vector("list", length(starts))
  alts <- vector("list", length(starts))
  for (w in seq_along(starts)) {
    s <- starts[w]; e <- s + window - 1L
    lo <- if (is.finite(anchor)) max(1L, s - as.integer(anchor)) else 1L
    hi <- if (is.finite(anchor)) min(Tn, e + as.integer(anchor)) else Tn
    sub <- ord[lo:hi]
    idx <- (s - lo + 1L):(e - lo + 1L)
    inner <- ord[s:e]
    cands <- lapply(seq_len(nrow(perms)), function(j) {
      o <- sub; o[idx] <- inner[perms[j, ]]; o
    })
    sc <- .hsmap_score_orders(x, cands, h_all, dam_id, epsilon, lambda, tol,
                              maxit, r_starts, refine_top, no_linkage_r,
                              threads)
    if (all(is.na(sc$logLik)))
      stop("No permutation converged in window ", w, ".", call. = FALSE)
    b <- which.max(sc$logLik)
    base_ll <- sc$logLik[1L]
    near <- which(sc$logLik >= sc$logLik[b] - tie_tol)
    rows[[w]] <- data.frame(
      window = w, start = s, end = e,
      baseline_logLik = base_ll, best_logLik = sc$logLik[b],
      dLL = sc$logLik[b] - base_ll,
      identity_best = b == 1L,
      n_better = sum(sc$logLik > base_ll + 1e-9, na.rm = TRUE),
      n_near_optimal = length(near),
      baseline_cM = sc$cM[1L], best_cM = sc$cM[b],
      dcM = sc$cM[b] - sc$cM[1L],
      converged = all(sc$converged),
      r_start = sc$r_start[b], ll_spread = sc$ll_spread[b],
      supported = (sc$logLik[b] - base_ll) > min_delta && b != 1L,
      stringsAsFactors = FALSE)
    best_orders[[w]] <- inner[perms[b, ]]
    top <- utils::head(order(-sc$logLik), min(5L, length(sc$logLik)))
    alts[[w]] <- data.frame(
      rank = seq_along(top),
      order = vapply(top, function(j) paste(inner[perms[j, ]], collapse = ">"), ""),
      logLik = sc$logLik[top], dLL_baseline = sc$logLik[top] - base_ll,
      local_cM = sc$cM[top], stringsAsFactors = FALSE)
    if (isTRUE(verbose))
      cat(sprintf("  window %d [%d-%d] dLL %.3f%s\n", w, s, e,
                  rows[[w]]$dLL, if (rows[[w]]$supported) " *" else ""))
  }
  W <- do.call(rbind, rows)

  applied <- NULL
  if (isTRUE(apply)) {
    ## experimental single pass: accept supported, non-overlapping windows,
    ## strongest first. No re-evaluation after acceptance.
    applied <- ord
    acc <- W[W$supported, , drop = FALSE]
    acc <- acc[order(-acc$dLL), , drop = FALSE]
    used <- rep(FALSE, Tn)
    for (i in seq_len(nrow(acc))) {
      rng <- acc$start[i]:acc$end[i]
      if (any(used[rng])) next
      applied[rng] <- best_orders[[acc$window[i]]]
      used[rng] <- TRUE
    }
  }

  structure(list(windows = W, best_orders = best_orders, alternatives = alts,
                 order = ord, applied_order = applied,
                 settings = list(window = window, by = by, anchor = anchor,
                                 dam = dam, epsilon = epsilon, lambda = lambda,
                                 r_starts = r_starts, refine_top = refine_top,
                                 tie_tol = tie_tol, min_delta = min_delta,
                                 no_linkage_r = no_linkage_r, apply = apply)),
            class = "HSMap.ripple")
}

## all permutations of seq_len(n), identity first
.hsmap_perms <- function(n) {
  rec <- function(k) {
    if (k == 1L) return(matrix(1L, 1L, 1L))
    sub <- rec(k - 1L)
    do.call(rbind, lapply(seq_len(k), function(i)
      cbind(i, matrix(c(seq_len(k)[-i])[sub], nrow(sub)))))
  }
  p <- rec(as.integer(n))
  storage.mode(p) <- "integer"
  ident <- which(apply(p, 1L, function(z) all(z == seq_len(n))))
  rbind(p[ident, , drop = FALSE], p[-ident, , drop = FALSE])
}

#' @rdname hsmap-print
#' @exportS3Method print HSMap.try
print.HSMap.try <- function(x, ...) {
  p <- x$profile
  cat("HSMap.try (placement support)\n")
  cat("  Marker            : ", x$marker, "\n", sep = "")
  cat("  Scaffold          : ", length(x$framework), " markers, ",
      nrow(p), " insertion positions\n", sep = "")
  cat("  Current position  : ", x$current_pos %|NA|% "not in supplied order",
      "\n", sep = "")
  cat("  Best position     : ", x$best_pos,
      if (p$at_edge[x$best_pos]) " (edge)" else "", "\n", sep = "")
  if (!is.na(x$displacement))
    cat("  Displacement      : ", x$displacement,
        sprintf("  (best - current = %.3f logLik)", p$dLL_current[x$best_pos]),
        "\n", sep = "")
  cat("  Near-optimal set  : ", length(x$near_optimal), " position(s) within ",
      x$settings$tie_tol, " logLik",
      if (length(x$near_optimal) > 1L)
        sprintf(" [%d-%d]", min(x$near_optimal), max(x$near_optimal)) else "",
      "\n", sep = "")
  cat("  Placement         : ", x$placement, "\n", sep = "")
  if (any(p$n_no_linkage > 0L))
    cat("  No-linkage        : present in ", sum(p$n_no_linkage > 0L),
        " placement(s)\n", sep = "")
  if (!all(p$converged))
    cat("  WARNING           : ", sum(!p$converged), " placement(s) did not converge\n", sep = "")
  invisible(x)
}

#' @rdname hsmap-print
#' @exportS3Method print HSMap.ripple
print.HSMap.ripple <- function(x, ...) {
  W <- x$windows
  cat("HSMap.ripple (local order support)\n")
  cat("  Order            : ", length(x$order), " markers\n", sep = "")
  cat("  Windows          : ", nrow(W), " of size ", x$settings$window,
      " (", factorial(x$settings$window), " permutations each, anchor ",
      x$settings$anchor, ")\n", sep = "")
  cat("  Baseline best    : ", sum(W$identity_best), " / ", nrow(W),
      " windows\n", sep = "")
  cat("  Supported change : ", sum(W$supported), " window(s) (dLL > ",
      x$settings$min_delta, ")\n", sep = "")
  if (any(W$supported)) {
    s <- W[W$supported, , drop = FALSE]
    s <- s[order(-s$dLL), , drop = FALSE]
    cat("  Strongest        : window ", s$window[1], " [", s$start[1], "-",
        s$end[1], "] dLL ", sprintf("%.2f", s$dLL[1]),
        ", local cM ", sprintf("%.2f -> %.2f", s$baseline_cM[1], s$best_cM[1]),
        "\n", sep = "")
  }
  cat("  Ambiguous        : ", sum(W$n_near_optimal > 1L),
      " window(s) with tied near-optimal orders\n", sep = "")
  if (!all(W$converged))
    cat("  WARNING          : non-convergence in ", sum(!W$converged),
        " window(s)\n", sep = "")
  if (!is.null(x$applied_order))
    cat("  applied_order    : present (experimental single pass)\n")
  invisible(x)
}

#' Plot a marker placement profile
#'
#' @param x An \code{HSMap.try} object.
#' @param ... Passed to \code{\link[graphics]{plot}}.
#' @return \code{x}, invisibly.
#' @exportS3Method plot HSMap.try
plot.HSMap.try <- function(x, ...) {
  p <- x$profile
  graphics::plot(p$position, p$dLL_best, type = "b", pch = 16,
                 xlab = "insertion position", ylab = "logLik - best",
                 main = paste0("try_marker: ", x$marker), ...)
  graphics::abline(h = -x$settings$tie_tol, lty = 3, col = "grey50")
  if (!is.na(x$current_pos))
    graphics::abline(v = x$current_pos, lty = 2, col = "steelblue")
  graphics::points(x$best_pos, p$dLL_best[x$best_pos], pch = 1, cex = 2.2,
                   col = "firebrick")
  invisible(x)
}
