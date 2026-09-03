#' Classify adjacent map intervals by multipoint/two-point concordance
#'
#' @description
#' The primary gap diagnostic for framework-map construction. For every
#' adjacent interval of a fitted map it reports the multipoint recombination
#' fraction alongside the two-point evidence for the same marker pair, and
#' assigns one of five descriptive statuses. It never infers \emph{why} an
#' interval is discordant and never modifies anything.
#'
#' @details
#' Statuses, in decision order:
#' \describe{
#'   \item{\code{unresolved}}{the map's relative phase at this interval is
#'     \code{NA} (no recombination estimate is meaningful).}
#'   \item{\code{no_linkage}}{multipoint \eqn{r \ge} \code{no_linkage_r}: a
#'     no-linkage boundary; its centimorgan distance is reported as \code{NA},
#'     never as a large finite number.}
#'   \item{\code{large_discordant}}{multipoint \eqn{r \ge} \code{large_r}
#'     while the two-point analysis of the \emph{same pair} indicates tight
#'     linkage with adequate support (\code{lod_r >= lod_min} and
#'     \eqn{r_{mp} - r_{2pt} \ge} \code{min_discordance}). The two analyses
#'     disagree; the reason is not inferred.}
#'   \item{\code{large_concordant}}{multipoint \eqn{r \ge} \code{large_r} and
#'     the two-point evidence does not contradict it. Such intervals may be
#'     real genomic distances and should not be "repaired" merely for being
#'     large.}
#'   \item{\code{normal}}{everything else.}
#' }
#'
#' @param map An \code{HSMap.map} (from \code{\link{hmm_map}}).
#' @param tpt An \code{HSMap.tpt} covering the map's markers.
#' @param order Optional marker order; defaults to \code{map$order}.
#' @param large_r Threshold defining an operationally large multipoint
#'   interval. Default \code{0.05}.
#' @param lod_min Minimum two-point linkage LOD for the pair to count as
#'   informative when calling discordance. Default \code{3}.
#' @param min_discordance Minimum \eqn{r_{mp} - r_{2pt}} to call discordance.
#'   Defaults to \code{large_r}.
#' @param no_linkage_r Multipoint \eqn{r} at/above which the interval is a
#'   no-linkage boundary. Default \code{0.499}.
#'
#' @return A data frame with one row per adjacent interval: \code{interval},
#'   \code{marker_left}, \code{marker_right}, \code{r_mp}, \code{cM}
#'   (Haldane; \code{NA} for \code{no_linkage}/\code{unresolved}),
#'   \code{r_2pt}, \code{lod_r}, \code{lod_ph}, \code{discordance}
#'   (\eqn{r_{mp} - r_{2pt}}) and \code{status}.
#'
#' @examples
#' \dontrun{
#' iv <- diagnose_map_intervals(map1, tpt, large_r = 0.05)
#' subset(iv, status == "large_discordant")
#' }
#' @seealso \code{\link{collapse_framework_blocks}} for the safe,
#'   user-directed repair of discordant regions.
#' @export
diagnose_map_intervals <- function(map, tpt, order = NULL,
                                   large_r = 0.05,
                                   lod_min = 3,
                                   min_discordance = large_r,
                                   no_linkage_r = 0.499) {
  if (!inherits(map, "HSMap.map")) stop("`map` must be an HSMap.map.", call. = FALSE)
  if (!inherits(tpt, "HSMap.tpt")) stop("`tpt` must be an HSMap.tpt.", call. = FALSE)
  order <- as.character(if (is.null(order)) map$order else order)
  if (!identical(order, as.character(map$order)))
    stop("`order` must match the fitted map's marker order.", call. = FALSE)
  miss <- setdiff(order, rownames(tpt$fit$r))
  if (length(miss))
    stop(length(miss), " marker(s) not present in `tpt`, e.g. ",
         paste(utils::head(miss, 3), collapse = ", "), call. = FALSE)
  if (!is.numeric(large_r) || large_r <= 0 || large_r >= 0.5)
    stop("`large_r` must be in (0, 0.5).", call. = FALSE)

  Tn <- length(order)
  r_mp <- as.numeric(map$fit$r)
  if (length(r_mp) != Tn - 1L)
    stop("length(map$fit$r) must equal length(order) - 1.", call. = FALSE)
  pv <- map$phase_vec
  ij <- cbind(seq_len(Tn - 1L), 2:Tn)
  R2 <- tpt$fit$r[order, order, drop = FALSE][ij]
  LR <- tpt$fit$lod_r[order, order, drop = FALSE][ij]
  LP <- if (!is.null(tpt$fit$lod_ph))
    tpt$fit$lod_ph[order, order, drop = FALSE][ij] else rep(NA_real_, Tn - 1L)

  unresolved <- if (length(pv) == Tn - 1L) is.na(pv) else rep(FALSE, Tn - 1L)
  nolink     <- !unresolved & r_mp >= no_linkage_r
  discord    <- !unresolved & !nolink & r_mp >= large_r &
    !is.na(R2) & !is.na(LR) & LR >= lod_min & (r_mp - R2) >= min_discordance
  large_conc <- !unresolved & !nolink & r_mp >= large_r & !discord

  status <- rep("normal", Tn - 1L)
  status[large_conc] <- "large_concordant"
  status[discord]    <- "large_discordant"
  status[nolink]     <- "no_linkage"
  status[unresolved] <- "unresolved"

  cM <- -50 * log(1 - 2 * pmin(r_mp, 0.49999))
  cM[nolink | unresolved] <- NA_real_

  data.frame(interval = seq_len(Tn - 1L),
             marker_left = order[-Tn], marker_right = order[-1L],
             r_mp = r_mp, cM = cM, r_2pt = R2, lod_r = LR, lod_ph = LP,
             discordance = r_mp - R2, status = status,
             stringsAsFactors = FALSE)
}

#' Collapse user-specified marker blocks into framework representatives
#'
#' @description
#' Safe, user-directed construction of a framework map. Each supplied block --
#' identified by the analyst from \code{\link{diagnose_map_intervals}} and the
#' pairwise-RF structure -- is represented in the multipoint HMM by a single
#' deterministic representative; the remaining block markers are recorded as
#' \code{attached} (never deleted) and inherit the representative's map
#' position. The repaired map is rebuilt from scratch
#' (\code{\link{phase_from_pairwise}} then \code{\link{hmm_map}}) and accepted
#' \strong{only} if an automatic verification passes; otherwise the original
#' map is returned unchanged.
#'
#' Block identification is deliberately \emph{not} automatic: validation on
#' real data showed automatic block detection is not reproducible enough for
#' unattended use, while collapse of correctly identified blocks is safe and
#' effective.
#'
#' @details
#' \strong{Representative rule.} Within each block, the marker with the
#' largest summed pairwise linkage LOD to the other block members; ties break
#' by original order position.
#'
#' \strong{Automatic verification.} The repaired map is accepted only if all
#' of the following hold:
#' \enumerate{
#'   \item every target interval that was large before collapse is now below
#'     \code{large_r}, or reduced to at most half its previous value;
#'   \item no interval near a block (within \code{margin} retained markers)
#'     that was below \code{large_r} before is at/above it after;
#'   \item every previously \code{large_concordant} interval still adjacent in
#'     the repaired map changes by at most \code{preserve_tol};
#'   \item the repaired EM fit converged;
#'   \item phase resolution does not degrade (no more unresolved intervals
#'     than before);
#'   \item at least \code{stability_frac} of intervals outside all blocks that
#'     are adjacent in both maps satisfy \eqn{|r_{after} - r_{before}| \le}
#'     \code{stability_tol};
#'   \item total map length does not increase by more than a factor of
#'     \code{max_len_increase}.
#' }
#' A failure of any check rejects the whole repair and returns the original
#' map, with the check table reporting what failed.
#'
#' @param x An \code{HSMap.data} object.
#' @param tpt An \code{HSMap.tpt} covering \code{order}.
#' @param order Character vector, the linkage-group marker order.
#' @param blocks A list of character vectors; each vector holds the marker IDs
#'   of one block. Every block must contain at least 2 markers, be contiguous
#'   in \code{order}, not touch the first or last marker, and blocks must not
#'   overlap.
#' @param dam Dam index/name for phasing and mapping. Default \code{1}.
#' @param epsilon Genotyping error rate for \code{\link{hmm_map}}.
#'   Default \code{0.05}.
#' @param tol,maxit EM convergence tolerance and iteration budget, passed to
#'   \code{\link{hmm_map}} for both the original and the repaired fit.
#'   Defaults \code{1e-6} and \code{1000}. Raising \code{maxit} does not weaken
#'   verification: check 4 still requires the repaired fit to have converged,
#'   it only gives the EM a larger budget to do so.
#' @param large_r,lod_min,min_discordance,no_linkage_r Passed to
#'   \code{\link{diagnose_map_intervals}} for the before/after classification.
#' @param stability_tol,stability_frac Outside-block stability requirement
#'   (check 6). Defaults \code{0.01} and \code{0.99}.
#' @param preserve_tol Maximum allowed change of a previously
#'   \code{large_concordant} interval (check 3). Default \code{0.02}.
#' @param margin Retained-marker neighbourhood searched for newly created
#'   large intervals (check 2). Default \code{10}.
#' @param max_len_increase Maximum tolerated ratio of repaired to original map
#'   length (check 7). Default \code{1.1}.
#' @param lg Optional label stamped into the marker table (e.g. \code{"LG10"}).
#' @param verbose Print progress. Default \code{TRUE}.
#'
#' @return An object of class \code{"HSMap.framework"}: a list with
#'   \code{map} (the accepted map: repaired if verification passed, otherwise
#'   the original), \code{accepted} (logical), \code{checks} (the seven
#'   verification results with values), \code{blocks} (per block: markers,
#'   representative, target interval r before/after), \code{marker_table}
#'   (\code{marker}, \code{framework_status} in
#'   \code{framework}/\code{representative}/\code{attached},
#'   \code{representative}, \code{lg}, \code{position_cM}; attached markers
#'   inherit their representative's position), \code{intervals_before},
#'   \code{intervals_after}, \code{map_original}, \code{map_repaired}, and
#'   \code{settings}.
#'
#' @examples
#' \dontrun{
#' iv <- diagnose_map_intervals(map10, tpt, large_r = 0.05)
#' subset(iv, status == "large_discordant")
#' ## blocks chosen by the analyst from the diagnostic + RF heatmap:
#' fw <- collapse_framework_blocks(dat, tpt, order = l$LG10,
#'                                 blocks = list(lg10_block1, lg10_block2))
#' fw$accepted
#' fw$marker_table
#' }
#' @export
collapse_framework_blocks <- function(x, tpt, order, blocks,
                                      dam = 1,
                                      epsilon = 0.05,
                                      tol = 1e-6,
                                      maxit = 1000L,
                                      large_r = 0.05,
                                      lod_min = 3,
                                      min_discordance = large_r,
                                      no_linkage_r = 0.499,
                                      stability_tol = 0.01,
                                      stability_frac = 0.99,
                                      preserve_tol = 0.02,
                                      margin = 10L,
                                      max_len_increase = 1.1,
                                      lg = NA_character_,
                                      verbose = TRUE) {
  if (!inherits(x, "HSMap.data")) stop("`x` must be an HSMap.data.", call. = FALSE)
  if (!inherits(tpt, "HSMap.tpt")) stop("`tpt` must be an HSMap.tpt.", call. = FALSE)
  order <- as.character(order)
  Tn <- length(order)
  if (Tn < 5L) stop("`order` needs at least 5 markers.", call. = FALSE)
  if (!is.list(blocks) || !length(blocks))
    stop("`blocks` must be a non-empty list of marker-ID vectors.", call. = FALSE)
  say <- function(...) if (isTRUE(verbose)) cat(sprintf(...), sep = "")

  ## ---- validate blocks -----------------------------------------------------
  bl_idx <- lapply(blocks, function(b) {
    b <- as.character(b)
    p <- match(b, order)
    if (anyNA(p))
      stop("Block marker(s) not in `order`: ",
           paste(utils::head(b[is.na(p)], 3), collapse = ", "), call. = FALSE)
    p <- sort(p)
    if (length(p) < 2L)
      stop("Each block needs at least 2 markers.", call. = FALSE)
    if (!all(diff(p) == 1L))
      stop("Block is not contiguous in `order` (positions ",
           paste(utils::head(p, 6), collapse = ","),
           "); split it into contiguous blocks.", call. = FALSE)
    if (p[1] <= 1L || p[length(p)] >= Tn)
      stop("A block may not include the first or last marker of the order.",
           call. = FALSE)
    p
  })
  all_idx <- unlist(bl_idx)
  if (anyDuplicated(all_idx)) stop("Blocks overlap.", call. = FALSE)

  ## ---- representatives (deterministic) -------------------------------------
  LRm <- tpt$fit$lod_r[order, order, drop = FALSE]
  reps <- vapply(bl_idx, function(p) {
    s <- rowSums(LRm[p, p, drop = FALSE], na.rm = TRUE)
    p[order(-s, p)[1L]]
  }, 0L)
  attached_idx <- sort(setdiff(all_idx, reps))

  fit_one <- function(mk) {
    ph <- phase_from_pairwise(tpt, order = mk, dam = dam)
    if (inherits(ph, "HSMap.phased.multi")) ph <- ph[[1L]]
    list(ph = ph,
         map = suppressWarnings(hmm_map(x, phased = ph, dam = dam,
                                        epsilon = epsilon, tol = tol,
                                        maxit = maxit)))
  }
  len_of <- function(m) {
    r <- as.numeric(m$fit$r)
    sum(-50 * log(1 - 2 * pmin(r[r < no_linkage_r], 0.49999)))
  }

  ## ---- original + repaired fits (full rebuilds) -----------------------------
  say("fitting original map (%d markers)...\n", Tn)
  f0 <- fit_one(order)
  iv0 <- diagnose_map_intervals(f0$map, tpt, large_r = large_r,
                                lod_min = lod_min,
                                min_discordance = min_discordance,
                                no_linkage_r = no_linkage_r)
  keep <- setdiff(seq_len(Tn), attached_idx)
  say("fitting repaired map (%d retained, %d attached)...\n",
      length(keep), length(attached_idx))
  f1 <- fit_one(order[keep])
  iv1 <- diagnose_map_intervals(f1$map, tpt, large_r = large_r,
                                lod_min = lod_min,
                                min_discordance = min_discordance,
                                no_linkage_r = no_linkage_r)
  r0 <- as.numeric(f0$map$fit$r); r1 <- as.numeric(f1$map$fit$r)

  ## span (max interval r) between nearest retained flanks of a block, in map m
  span_r <- function(m, lo, hi) {
    kp <- match(m$order, order)
    l <- kp[kp < lo]; h <- kp[kp > hi]
    if (!length(l) || !length(h)) return(NA_real_)
    p <- match(order[c(max(l), min(h))], m$order)
    max(as.numeric(m$fit$r)[min(p):(max(p) - 1L)])
  }
  before_r <- vapply(bl_idx, function(p) span_r(f0$map, min(p), max(p)), 0)
  after_r  <- vapply(bl_idx, function(p) span_r(f1$map, min(p), max(p)), 0)

  ## intervals adjacent in BOTH maps, outside every block
  in_block <- vapply(seq_len(Tn - 1L), function(i)
    any(vapply(bl_idx, function(p) i >= min(p) - 1L && i <= max(p), TRUE)), TRUE)
  both_adj <- which(!in_block)
  both_adj <- both_adj[vapply(both_adj, function(i) {
    a <- match(order[i], f1$map$order); b <- match(order[i + 1L], f1$map$order)
    !is.na(a) && !is.na(b) && b == a + 1L
  }, TRUE)]
  d_out <- vapply(both_adj, function(i)
    r1[match(order[i], f1$map$order)] - r0[i], 0)

  ## ---- the seven checks -----------------------------------------------------
  targeted <- before_r >= large_r
  chk1 <- !any(targeted) ||
    all(after_r[targeted] < large_r | after_r[targeted] <= before_r[targeted] / 2,
        na.rm = FALSE)
  near <- unique(unlist(lapply(bl_idx, function(p) {
    lo <- max(1L, min(p) - 1L - margin); hi <- min(Tn - 1L, max(p) + margin)
    intersect(both_adj, lo:hi)
  })))
  chk2 <- !length(near) || !any(
    r0[near] < large_r &
      vapply(near, function(i) r1[match(order[i], f1$map$order)], 0) >= large_r)
  lc <- intersect(which(iv0$status == "large_concordant"), both_adj)
  chk3 <- !length(lc) || all(abs(vapply(lc, function(i)
    r1[match(order[i], f1$map$order)], 0) - r0[lc]) <= preserve_tol)
  chk4 <- isTRUE(f1$map$fit$converged)
  chk5 <- sum(is.na(f1$ph$phase_vec)) <= sum(is.na(f0$ph$phase_vec))
  frac_stable <- if (length(d_out)) mean(abs(d_out) <= stability_tol) else 1
  chk6 <- frac_stable >= stability_frac
  len0 <- len_of(f0$map); len1 <- len_of(f1$map)
  chk7 <- len1 <= len0 * max_len_increase

  checks <- data.frame(
    check = c("target intervals resolved", "no new large interval nearby",
              "large_concordant preserved", "EM converged",
              "phase resolution not degraded", "outside intervals stable",
              "map length not inflated"),
    pass = c(chk1, chk2, chk3, chk4, chk5, chk6, chk7),
    value = c(paste(sprintf("%.3f->%.3f", before_r, after_r), collapse = "; "),
              sprintf("%d nearby intervals checked", length(near)),
              sprintf("%d intervals, tol %.3g", length(lc), preserve_tol),
              as.character(f1$map$fit$converged),
              sprintf("%d -> %d unresolved", sum(is.na(f0$ph$phase_vec)),
                      sum(is.na(f1$ph$phase_vec))),
              sprintf("%.4f of %d within %.3g", frac_stable, length(d_out),
                      stability_tol),
              sprintf("%.1f -> %.1f cM (limit x%.2f)", len0, len1,
                      max_len_increase)),
    stringsAsFactors = FALSE)
  accepted <- all(checks$pass)
  say("verification: %s\n", if (accepted) "ACCEPTED (7/7 checks passed)"
      else paste0("REJECTED (failed: ",
                  paste(checks$check[!checks$pass], collapse = "; "),
                  ") -- returning the original map"))

  ## ---- marker table ---------------------------------------------------------
  final <- if (accepted) f1$map else f0$map
  pos <- get_map(final, "haldane")
  status <- rep("framework", Tn)
  rep_of <- rep(NA_character_, Tn)
  if (accepted) {
    for (k in seq_along(bl_idx)) {
      status[bl_idx[[k]]] <- "attached"
      status[reps[k]] <- "representative"
      rep_of[bl_idx[[k]]] <- order[reps[k]]
    }
  }
  pcm <- unname(pos[order])
  if (accepted && length(attached_idx))
    pcm[attached_idx] <- unname(pos[rep_of[attached_idx]])
  marker_table <- data.frame(marker = order, framework_status = status,
                             representative = rep_of, lg = lg,
                             position_cM = pcm, stringsAsFactors = FALSE)

  blocks_tab <- data.frame(
    block = seq_along(bl_idx),
    start = vapply(bl_idx, min, 0L), end = vapply(bl_idx, max, 0L),
    n_markers = lengths(bl_idx),
    representative = order[reps],
    markers = vapply(blocks, paste, "", collapse = ";"),
    target_r_before = before_r, target_r_after = after_r,
    applied = accepted, stringsAsFactors = FALSE)

  out <- list(map = final, accepted = accepted, checks = checks,
              blocks = blocks_tab, marker_table = marker_table,
              intervals_before = iv0, intervals_after = iv1,
              map_original = f0$map, map_repaired = f1$map,
              settings = list(large_r = large_r, lod_min = lod_min,
                              min_discordance = min_discordance,
                              stability_tol = stability_tol,
                              stability_frac = stability_frac,
                              preserve_tol = preserve_tol, margin = margin,
                              max_len_increase = max_len_increase,
                              epsilon = epsilon, tol = tol, maxit = maxit))
  class(out) <- "HSMap.framework"
  out
}

#' @export
print.HSMap.framework <- function(x, ...) {
  cat("HSMap.framework (user-directed block collapse)\n")
  cat(sprintf("  verdict     : %s\n",
              if (x$accepted) "ACCEPTED (7/7 checks)" else "REJECTED -- original map returned"))
  st <- table(factor(x$marker_table$framework_status,
                     c("framework", "representative", "attached")))
  cat(sprintf("  markers     : %d framework, %d representative, %d attached\n",
              st[1], st[2], st[3]))
  cat(sprintf("  blocks      : %d supplied\n", nrow(x$blocks)))
  lb <- sum(x$intervals_before$status %in% c("large_concordant", "large_discordant"))
  la <- sum(x$intervals_after$status %in% c("large_concordant", "large_discordant"))
  cat(sprintf("  large intervals: %d -> %d\n", lb, la))
  if (!x$accepted)
    cat(sprintf("  failed      : %s\n",
                paste(x$checks$check[!x$checks$pass], collapse = "; ")))
  invisible(x)
}
