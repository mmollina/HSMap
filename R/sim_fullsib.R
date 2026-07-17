#' Simulate known-sire full-sib and mixed (open-pollinated + full-sib) data
#'
#' @description
#' Generate crosses with genotyped mothers and (optionally) genotyped sires, with
#' independent maternal and paternal recombination maps. For a known-sire cross the
#' paternal alleles are transmitted through a linked paternal meiosis (NOT drawn
#' independently across markers); for an open-pollinated cross the paternal allele is
#' drawn independently per marker from a pollen-pool frequency (the current OP model).
#' Parents referenced by the same ID across crosses are simulated once, so a shared
#' mother (or sire) carries one genotype, one phase, and one map.
#'
#' @param n_markers Integer number of markers (z).
#' @param crosses A data.frame describing crosses with columns \code{mother},
#'   \code{father} (\code{NA} = open-pollinated / unknown sire), and \code{n}
#'   (offspring count). If \code{NULL}, one full-sib cross \code{M1 x S1} with 100
#'   offspring is simulated.
#' @param r_m,r_p Optional numeric length \code{z-1} maternal / paternal recombination
#'   vectors. If \code{NULL}, filled with \code{r_const_m} / \code{r_const_p}.
#' @param r_const_m,r_const_p Constant maternal / paternal recombination fraction used
#'   when \code{r_m}/\code{r_p} are \code{NULL}.
#' @param epsilon Symmetric genotyping error rate applied to offspring calls.
#' @param missing_rate Probability an offspring call is set missing (NA).
#' @param maternal_mode,paternal_mode Parent genotype model, \code{"all_het"} (default;
#'   maximally informative) or \code{"hwe"} at allele-A frequency \code{parent_pA}.
#' @param parent_pA Allele-A frequency for HWE parent genotypes.
#' @param repulsion_rate Fraction of intervals simulated in repulsion (per parent phase).
#' @param op_paternal_pA Pollen-pool allele-A frequency for open-pollinated crosses
#'   (independent per marker), scalar or length-\code{z}.
#' @param seed Optional RNG seed.
#'
#' @return A list with:
#'   \itemize{
#'     \item \code{data}: an \code{HSMap.data} object (cross-aware fields populated:
#'       \code{G_list}, \code{M_list}, \code{F_list}, \code{crosses},
#'       \code{parent_genotypes}, \code{cross_table}, \code{alleles}, \code{pedigree}).
#'     \item \code{truth}: \code{r_m}, \code{r_p}; \code{parents} (per id: genotype,
#'       phase_vec, hap = 2 x z); \code{maternal_paths} and \code{paternal_paths}
#'       (per cross: offspring x z homolog indices in \{1,2\}; paternal is \code{NA}
#'       for OP crosses); \code{cross_table}.
#'     \item \code{markers}, \code{params}.
#'   }
#' @section Lifecycle - experimental (not part of the current paper):
#' The known-sire / full-sib functions are \strong{experimental} and are \strong{not}
#' part of the published open-pollinated method. They currently support only
#' \strong{oracle parental haplotypes} and are \strong{not ready for automatic real-data
#' mapping}; the API may change without a deprecation cycle. See
#' \code{dev/known_sire_design.md}.
#' @importFrom stats runif rbinom
#' @export
sim_fullsib <- function(n_markers = 20L,
                        crosses = NULL,
                        r_m = NULL, r_p = NULL,
                        r_const_m = 0.1, r_const_p = 0.1,
                        epsilon = 0.0, missing_rate = 0.0,
                        maternal_mode = c("all_het", "hwe"),
                        paternal_mode = c("all_het", "hwe"),
                        parent_pA = 0.5, repulsion_rate = 0.3,
                        op_paternal_pA = 0.5, seed = NULL) {
  maternal_mode <- match.arg(maternal_mode)
  paternal_mode <- match.arg(paternal_mode)
  if (!is.null(seed)) set.seed(seed)
  z <- as.integer(n_markers)
  if (z < 2L) stop("`n_markers` must be >= 2.")
  markers <- sprintf("M%04d", seq_len(z))

  if (is.null(crosses))
    crosses <- data.frame(mother = "M1", father = "S1", n = 100L,
                          stringsAsFactors = FALSE)
  crosses$mother <- as.character(crosses$mother)
  crosses$father <- ifelse(is.na(crosses$father), NA_character_, as.character(crosses$father))
  if (is.null(crosses$n)) crosses$n <- 100L

  r_m <- if (is.null(r_m)) rep(r_const_m, z - 1L) else r_m
  r_p <- if (is.null(r_p)) rep(r_const_p, z - 1L) else r_p
  if (length(r_m) != z - 1L || length(r_p) != z - 1L)
    stop("`r_m`/`r_p` must have length n_markers - 1.")
  op_pA <- if (length(op_paternal_pA) == 1L) rep(op_paternal_pA, z) else op_paternal_pA

  # ---- simulate each referenced parent ONCE ---------------------------------
  gen_geno <- function(mode) {
    if (mode == "all_het") rep(1L, z)
    else as.integer(rbinom(z, 2, parent_pA))         # HWE
  }
  gen_phase <- function() as.integer(rbinom(z - 1L, 1, 1 - repulsion_rate))  # 1=coupling

  mothers <- unique(crosses$mother)
  sires   <- unique(stats::na.omit(crosses$father))
  parents <- list()
  make_parent <- function(id, mode) {
    if (!is.null(parents[[id]])) return(parents[[id]])
    g  <- gen_geno(mode)
    pv <- gen_phase()
    list(genotype = stats::setNames(g, markers), phase_vec = pv,
         hap = phase_to_haplotypes(g, pv))
  }
  for (m in mothers) parents[[m]] <- make_parent(m, maternal_mode)
  for (s in sires)   parents[[s]] <- make_parent(s, paternal_mode)

  # ---- meiosis: transmitted homolog path along the chromosome ---------------
  meiosis_path <- function(nn, rr) {
    P <- matrix(0L, nn, z)
    P[, 1L] <- sample.int(2L, nn, replace = TRUE)     # 1 or 2 at first marker
    for (k in 2:z) {
      sw <- rbinom(nn, 1, rr[k - 1L]) == 1L           # recombination -> switch homolog
      P[, k] <- ifelse(sw, 3L - P[, k - 1L], P[, k - 1L])
    }
    P
  }
  err_missing <- function(G) {
    if (epsilon > 0) {
      hit <- matrix(runif(length(G)) < epsilon, nrow(G), ncol(G))
      if (any(hit)) {
        alt <- matrix(sample(0:2, length(G), replace = TRUE), nrow(G), ncol(G))
        # ensure a DIFFERENT genotype
        same <- hit & (alt == G)
        while (any(same, na.rm = TRUE)) {
          alt[same] <- sample(0:2, sum(same), replace = TRUE)
          same <- hit & (alt == G)
        }
        G[hit] <- alt[hit]
      }
    }
    if (missing_rate > 0) {
      miss <- matrix(runif(length(G)) < missing_rate, nrow(G), ncol(G))
      G[miss] <- NA_integer_
    }
    G
  }

  # ---- build crosses --------------------------------------------------------
  crs <- list(); F_list <- list(); G_list <- list(); M_list <- list()
  mat_paths <- list(); pat_paths <- list(); ct_rows <- list()
  ped_rows <- list(); pr <- 1L
  # founders in pedigree
  parent_ids <- names(parents)
  for (pid in parent_ids) {
    ped_rows[[pr]] <- data.frame(id = pid, mother = NA_character_, father = NA_character_,
                                 generation = 1L, family_id = pid, stringsAsFactors = FALSE); pr <- pr + 1L
  }

  for (i in seq_len(nrow(crosses))) {
    mom <- crosses$mother[i]; fat <- crosses$father[i]; nn <- as.integer(crosses$n[i])
    known <- !is.na(fat)
    cid <- if (known) paste(mom, fat, sep = "__x__") else mom
    Hm <- parents[[mom]]$hap

    Pm <- meiosis_path(nn, r_m)                       # maternal homolog path
    mat_allele <- matrix(0L, nn, z)
    for (k in seq_len(z)) mat_allele[, k] <- Hm[cbind(Pm[, k], k)]

    if (known) {
      Hp <- parents[[fat]]$hap
      Pp <- meiosis_path(nn, r_p)                     # LINKED paternal meiosis
      pat_allele <- matrix(0L, nn, z)
      for (k in seq_len(z)) pat_allele[, k] <- Hp[cbind(Pp[, k], k)]
      ftype <- "known_sire_genotyped"
    } else {
      Pp <- NULL
      pat_allele <- matrix(rbinom(nn * z, 1, rep(op_pA, each = nn)), nn, z)  # independent per marker
      ftype <- "open_pollinated"
    }

    G <- mat_allele + pat_allele                      # dosage 0/1/2
    storage.mode(G) <- "integer"
    G <- err_missing(G)
    kids <- sprintf("%s_off%03d", cid, seq_len(nn))
    rownames(G) <- kids; colnames(G) <- markers

    Mvec <- parents[[mom]]$genotype
    Fvec <- if (known) parents[[fat]]$genotype else stats::setNames(rep(NA_integer_, z), markers)

    crs[[cid]] <- list(cross_id = cid, mother_id = mom,
                       father_id = if (known) fat else HSMAP_UNKNOWN_SIRE,
                       family_type = ftype, offspring = kids,
                       M = Mvec, F = Fvec, G = G)
    F_list[[cid]] <- Fvec; G_list[[cid]] <- G; M_list[[cid]] <- Mvec
    mat_paths[[cid]] <- Pm; pat_paths[[cid]] <- Pp
    ct_rows[[i]] <- data.frame(cross_id = cid, mother_id = mom,
                               father_id = if (known) fat else HSMAP_UNKNOWN_SIRE,
                               family_type = ftype, n_offspring = nn,
                               mother_genotyped = TRUE, father_genotyped = known,
                               stringsAsFactors = FALSE)
    for (kk in kids) {
      ped_rows[[pr]] <- data.frame(id = kk, mother = mom,
                                   father = if (known) fat else NA_character_,
                                   generation = 2L, family_id = cid,
                                   stringsAsFactors = FALSE); pr <- pr + 1L
    }
  }
  cross_table <- do.call(rbind, ct_rows)
  pedigree <- do.call(rbind, ped_rows)
  parent_genotypes <- lapply(parents, function(p) as.integer(p$genotype))
  alleles <- data.frame(marker_id = markers, REF = "A", ALT = "B",
                        chrom = 1L, position = seq_len(z), stringsAsFactors = FALSE)

  dat <- list(G_list = G_list, M_list = M_list, alleles = alleles, pedigree = pedigree,
              stats = NULL, cross_table = cross_table, crosses = crs,
              parent_genotypes = parent_genotypes, F_list = F_list)
  class(dat) <- "HSMap.data"

  truth <- list(
    r_m = r_m, r_p = r_p,
    parents = parents,
    maternal_paths = mat_paths, paternal_paths = pat_paths,
    cross_table = cross_table
  )
  list(data = dat, truth = truth, markers = markers,
       params = list(n_markers = z, epsilon = epsilon, missing_rate = missing_rate,
                     maternal_mode = maternal_mode, paternal_mode = paternal_mode,
                     repulsion_rate = repulsion_rate))
}
