#' Per-child haplotype posteriors (gamma) from an HSMap map
#'
#' This wrapper builds inputs for \code{gamma_cpp()} from an \code{HSMap.data}
#' object and a map returned by \code{hmm_map()}. It supports both single-
#' dam maps and multi-dam maps (class \code{"HSMap.map.multi"} with a
#' \code{$per_dam} list).
#'
#' @param dat An \code{HSMap.data} object.
#' @param map One of: a single-dam map (class \code{"HSMap.map"}); a joint shared
#'   map (class \code{"HSMap.map.joint"} from \code{\link{hmm_map_joint}} or
#'   \code{hmm_map(method = "joint")}), decoded per dam using the shared \code{r}
#'   with each dam's phase and paternal emissions; or a multi-dam consensus map
#'   (class \code{"HSMap.map.multi"}) with element \code{$per_dam}.
#' @param dam Which dam(s) to extract. Options:
#'   \itemize{
#'     \item \code{NULL}: if \code{map} is single-dam, that dam is used;
#'           if \code{map} is multi-dam, defaults to \code{"all"}.
#'     \item character scalar: dam name (must exist in \code{names(map$per_dam)}).
#'     \item integer scalar: 1-based index into \code{map$per_dam}.
#'     \item \code{"all"}: compute for every dam in \code{map$per_dam}.
#'   }
#' @param epsilon Optional override for emission miscall rate. If \code{NULL},
#'   uses \code{fit$epsilon} from the dam’s map.
#'
#' @return
#' \itemize{
#'   \item For a single dam: an \code{HSMap.gamma} object with elements
#'         \code{gamma} \code{[2, T, n]}, \code{hap0_is_A}, \code{markers},
#'         \code{samples}, and \code{dam}.
#'   \item For multiple dams: a named list of \code{HSMap.gamma} objects with
#'         class \code{"HSMap.gamma.multi"}.
#' }
#' @export
calc_haploprob <- function(dat, map, dam = NULL, epsilon = NULL) {
  `%||%` <- function(a, b) if (is.null(a)) b else a

  is_single_map <- function(x) {
    is.list(x) && !is.null(x$order) && !is.null(x$phase_vec) && !is.null(x$fit)
  }
  is_multi_map <- function(x) {
    inherits(x, "HSMap.map.multi") || (!is.null(x$per_dam) && is.list(x$per_dam))
  }
  is_joint_map <- function(x) {
    inherits(x, "HSMap.map.joint") ||
      (is.list(x) && !is.null(x$phase_list) && !is.null(x$fit) && !is.null(x$order))
  }

  # dispatch on map type
  if (is_joint_map(map)) {
    return(.calc_haploprob_joint(dat = dat, map = map, dam = dam, epsilon = epsilon))
  }

  if (is_single_map(map)) {
    # single-dam path (ignore dam argument if given)
    return(.calc_haploprob_one(dat = dat, map_entry = map, dam_label = map$dam %||% "Dam1", epsilon = epsilon))
  }

  if (!is_multi_map(map)) {
    stop("`map` must be a single-dam HSMap.map or a multi-dam HSMap.map.multi with $per_dam.")
  }

  per <- map$per_dam
  dam_names <- names(per)

  # default: if multi-dam and dam not provided, run for all
  if (is.null(dam)) dam <- "all"

  if (identical(dam, "all")) {
    out <- lapply(seq_along(per), function(j) {
      dj <- dam_names[j] %||% paste0("Dam", j)
      .calc_haploprob_one(dat = dat, map_entry = per[[j]], dam_label = dj, epsilon = epsilon)
    })
    names(out) <- dam_names %||% paste0("Dam", seq_along(per))
    class(out) <- "HSMap.gamma.multi"
    return(out)
  }

  # single dam by name or index
  if (is.character(dam) && length(dam) == 1L) {
    j <- match(dam, dam_names)
    if (is.na(j)) stop("Dam '", dam, "' not found in map$per_dam. Valid: ", paste(dam_names, collapse = ", "))
    return(.calc_haploprob_one(dat = dat, map_entry = per[[j]], dam_label = dam, epsilon = epsilon))
  }
  if (is.numeric(dam) && length(dam) == 1L) {
    j <- as.integer(dam)
    if (j < 1L || j > length(per)) stop("`dam` index out of range [1..", length(per), "].")
    dj <- dam_names[j] %||% paste0("Dam", j)
    return(.calc_haploprob_one(dat = dat, map_entry = per[[j]], dam_label = dj, epsilon = epsilon))
  }

  stop("`dam` must be NULL, 'all', a single name, or a single 1-based index.")
}

# Internal worker: compute gamma for a single dam entry (HSMap.map)
.calc_haploprob_one <- function(dat, map_entry, dam_label, epsilon = NULL) {
  if (!inherits(dat, "HSMap.data")) stop("`dat` must be an HSMap.data object.")
  if (!is.list(map_entry) || is.null(map_entry$order) || is.null(map_entry$phase_vec) || is.null(map_entry$fit)) {
    stop("`map_entry` must be a single-dam HSMap.map with fields `order`, `phase_vec`, and `fit`.")
  }

  ord <- map_entry$order
  ph  <- map_entry$phase_vec
  fit <- map_entry$fit

  T <- length(ord)
  if (length(ph) != T - 1L) stop("`phase_vec` length must equal length(order) - 1.")

  # choose emissions used in the last E-step
  pi_emis <- fit$pi_emission
  if (is.null(pi_emis)) pi_emis <- fit$pi
  if (is.null(pi_emis)) stop("`fit` must contain `pi` or `pi_emission` (3 x T).")

  # r and epsilon
  r   <- as.numeric(fit$r)
  eps <- epsilon %||% fit$epsilon %||% 1e-3
  if (length(r) != T - 1L) stop("`fit$r` must have length length(order) - 1.")

  # find dam matrices in HSMap.data
  if (is.null(dat$G_list) || is.null(dat$M_list)) stop("`dat` must have G_list and M_list.")
  dam_idx <- if (!is.null(names(dat$G_list)) && dam_label %in% names(dat$G_list)) {
    match(dam_label, names(dat$G_list))
  } else if (length(dat$G_list) == 1L) {
    1L
  } else {
    stop("Dam '", dam_label, "' not found in `dat$G_list`. Available: ",
         paste(names(dat$G_list), collapse = ", "))
  }

  G_fam <- dat$G_list[[dam_idx]]   # offspring x markers
  M_fam <- dat$M_list[[dam_idx]]   # named vector
  if (is.null(colnames(G_fam))) stop("`dat$G_list[[dam]]` must have marker column names.")
  if (is.null(names(M_fam)))     stop("`dat$M_list[[dam]]` must be a named vector by marker.")

  # alignment checks
  miss_G <- setdiff(ord, colnames(G_fam))
  miss_M <- setdiff(ord, names(M_fam))
  if (length(miss_G)) stop("Markers missing in G for dam '", dam_label, "': ",
                           paste(utils::head(miss_G, 10), collapse = ", "),
                           if (length(miss_G) > 10) " ...")
  if (length(miss_M)) stop("Markers missing in M for dam '", dam_label, "': ",
                           paste(utils::head(miss_M, 10), collapse = ", "),
                           if (length(miss_M) > 10) " ...")

  # subset and align
  G <- as.matrix(G_fam[, ord, drop = FALSE])   # n x T
  storage.mode(G) <- "integer"
  M <- as.integer(M_fam[ord])

  # align pi_emis to ord if it has colnames
  if (is.matrix(pi_emis) && !is.null(colnames(pi_emis))) {
    miss_pi <- setdiff(ord, colnames(pi_emis))
    if (length(miss_pi)) stop("Markers missing in fit$pi/_emission: ",
                              paste(utils::head(miss_pi, 10), collapse = ", "),
                              if (length(miss_pi) > 10) " ...")
    pi_emis <- pi_emis[, ord, drop = FALSE]
  } else {
    if (ncol(pi_emis) != length(ord))
      stop("`pi`/`pi_emission` must have T columns (T = length(order)).")
  }

  # call C++ (returns [2, T, n] with attr 'hap0_is_A')
  gam <- gamma_cpp(G = G, M = M, phase_vec = as.integer(ph), r = r,
                   pi_emis = pi_emis, epsilon = eps)

  out <- list(
    gamma       = gam,
    hap0_is_A   = attr(gam, "hap0_is_A"),
    markers     = ord,
    samples     = rownames(G),
    dam         = dam_label
  )
  class(out) <- "HSMap.gamma"
  out
}

# Internal: per-child haplotype posteriors for each dam of a joint shared map
# (class "HSMap.map.joint"). The map carries one shared `r`, one shared `order`,
# and per-dam `phase_list` / paternal emissions; we decode each dam using the
# shared `r` with that dam's phase and paternal emission table.
.calc_haploprob_joint <- function(dat, map, dam = NULL, epsilon = NULL) {
  `%||%` <- function(a, b) if (is.null(a)) b else a
  if (!inherits(dat, "HSMap.data")) stop("`dat` must be an HSMap.data object.")
  if (is.null(map$order) || is.null(map$fit) || is.null(map$phase_list))
    stop("`map` must be a joint shared map with fields `order`, `phase_list`, `fit`.")

  fit <- map$fit
  dam_names <- map$dams %||% names(map$phase_list)
  if (is.null(dam_names)) stop("Joint map has no dam names.")

  # resolve requested dam(s)
  if (is.null(dam)) dam <- "all"
  if (identical(dam, "all")) {
    dsel <- dam_names
  } else if (is.character(dam)) {
    if (!all(dam %in% dam_names))
      stop("Dam(s) not in joint map: ", paste(setdiff(dam, dam_names), collapse = ", "))
    dsel <- dam
  } else if (is.numeric(dam)) {
    idx <- as.integer(dam)
    if (any(idx < 1L | idx > length(dam_names))) stop("`dam` index out of range.")
    dsel <- dam_names[idx]
  } else stop("`dam` must be NULL, 'all', a name/vector of names, or index/indices.")

  # per-dam emission table (3 x T): two_locus uses pi_emission_list, else pi_list
  emis_for_dam <- function(d_name) {
    em <- if (!is.null(fit$pi_emission_list)) fit$pi_emission_list[[d_name]]
          else if (!is.null(fit$pi_list))     fit$pi_list[[d_name]]
          else stop("Joint fit must contain `pi_list` or `pi_emission_list`.")
    em <- as.matrix(em)
    if (is.null(colnames(em))) colnames(em) <- map$order
    em
  }

  one <- function(d_name) {
    entry <- list(
      order     = map$order,
      phase_vec = as.integer(map$phase_list[[d_name]]),
      fit = list(r          = as.numeric(fit$r),
                 pi_emission = emis_for_dam(d_name),
                 epsilon     = fit$epsilon),
      dam = d_name
    )
    .calc_haploprob_one(dat = dat, map_entry = entry, dam_label = d_name, epsilon = epsilon)
  }

  if (length(dsel) == 1L) return(one(dsel))
  res <- lapply(dsel, one)
  names(res) <- dsel
  class(res) <- "HSMap.gamma.multi"
  res
}

