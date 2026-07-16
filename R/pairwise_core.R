#' Pairwise recombination fraction (RF) estimation across markers (parallel)
#'
#' @description
#' Estimates pairwise recombination fractions in a maternal-only setting across multiple families,
#' wrapping the C++ backend \code{pairwise_rf_estimation_multi_parallel_cpp()}.
#'
#' @param G_list List of integer genotype matrices (offspring x markers), values in \code{\{0,1,2,NA\}}.
#' @param M_list List of integer maternal genotype vectors (length T), values in \code{\{0,1,2,NA\}}.
#'   \code{M_list[[1]]} must be named to define the canonical marker order.
#' @param lambda Numeric scalar (default \code{20}); total pseudocount for the
#'   dam-specific paternal gametic-frequency prior (\code{alpha = lambda * q0},
#'   \code{beta = lambda * (1 - q0)}).
#' @param q0 Numeric prior mean for the paternal gametic frequency (default
#'   \code{0.5}); the pseudocount target that \code{q_k^(d)} is shrunk toward.
#' @param tol Numeric tolerance for the r optimizer's bounded local refinement
#'   (default \code{1e-6}).
#' @param maxit Integer maximum iterations for the local refinement (default \code{200}).
#' @param tiny Numeric floor to avoid \code{log(0)} (default \code{1e-12}).
#' @param share_q_across_dams Logical; if \code{FALSE} (default) \code{q_k^(d)} is
#'   estimated per dam; if \code{TRUE} the AA/aa counts are pooled across dams to a
#'   single per-marker \code{q}.
#' @param return_diagnostics Logical; if \code{TRUE}, additionally return a
#'   \code{diagnostics} list of per-pair optimizer and count matrices.
#' @param verbose Logical; print a one-line header (default \code{FALSE}).
#' @param n_threads Integer number of threads (positive).
#'
#' @return A list with matrices \code{r}, \code{lod_r} (raw LOD vs the exact
#' \code{r = 0.5} null), \code{lod_ph}, \code{logLik}, and \code{no_linkage} (all
#' T x T); \code{mom_phase_list}; \code{lod_ph_list} (per-dam phase-LOD matrices,
#' whose elementwise sum is \code{lod_ph}); \code{q_list} (per-dam per-marker
#' \code{q_k^(d)}, \code{NA} where the dam is not heterozygous); \code{optimizer},
#' \code{n_grid}; and, if \code{return_diagnostics}, a \code{diagnostics} list.
#'
#' @keywords internal
#' @noRd
cpp_pairwise_rf <- function(G_list,
                            M_list,
                            lambda = 20,
                            q0 = 0.5,
                            tol = 1e-6,
                            maxit = 200,
                            tiny = 1e-12,
                            share_q_across_dams = FALSE,
                            return_diagnostics = FALSE,
                            verbose = FALSE,
                            n_threads = getOption("HSMap.n_threads", 1L)) {

  # basic check
  if (!is.list(G_list) || !is.list(M_list))
    stop("`G_list` and `M_list` must be lists of the same length.")
  if (length(G_list) != length(M_list))
    stop("`G_list` and `M_list` must have the same length (one entry per dam/population).")

  # --- validate scalar hyper-parameters --------------------------------------
  if (!is.numeric(lambda) || length(lambda) != 1L || !is.finite(lambda) || lambda < 0)
    stop("`lambda` must be a single finite, non-negative number.")
  if (!is.numeric(q0) || length(q0) != 1L || !is.finite(q0) || q0 < 0 || q0 > 1)
    stop("`q0` must be a single finite number in [0, 1].")
  if (!is.numeric(tol) || length(tol) != 1L || !is.finite(tol) || tol <= 0)
    stop("`tol` must be a single finite, positive number.")
  if (!is.numeric(maxit) || length(maxit) != 1L || !is.finite(maxit) ||
      maxit < 1 || maxit != round(maxit))
    stop("`maxit` must be a single positive integer.")
  if (!is.numeric(tiny) || length(tiny) != 1L || !is.finite(tiny) || tiny <= 0)
    stop("`tiny` must be a single finite, positive number.")

  if (!is.numeric(n_threads) || length(n_threads) != 1L || !is.finite(n_threads) || n_threads < 1)
    stop("`n_threads` must be a positive finite integer.")
  n_threads <- as.integer(n_threads)
  RcppParallel::setThreadOptions(numThreads = as.integer(n_threads))

  # ensure M_list[[1]] has names (marker IDs) — required by the C++ code
  if (length(M_list) == 0L) stop("`M_list` cannot be empty.")
  if (is.null(names(M_list[[1]])))
    stop("`M_list[[1]]` must be a named integer vector with marker IDs as names.")

  # (light) coercions
  G_list <- lapply(G_list, function(G) {
    Gm <- as.matrix(G)
    storage.mode(Gm) <- "integer"
    Gm
  })
  M_list <- lapply(M_list, function(M) {
    Mi <- as.integer(M)
    names(Mi) <- names(M)
    Mi
  })
  # optional header
  if (isTRUE(verbose)) {
    mk <- names(M_list[[1]])
    cat(sprintf(
      "[cpp_pairwise_rf] dams: %d | markers: %d | threads: %d\n",
      length(G_list), length(mk), n_threads
    ))
  }

  # ---- call C++ core --------------------------------------------------------
  res <- pairwise_rf_estimation_multi_parallel_cpp(
    G_list = G_list,
    M_list = M_list,
    lambda = lambda,
    q0 = q0,
    tol = tol,
    maxit = maxit,
    tiny = tiny,
    share_q_across_dams = share_q_across_dams,
    return_diagnostics = return_diagnostics,
    verbose = verbose
  )

  # ---- post-process / class -------------------------------------------------
  markers <- colnames(res$r)
  out <- list(
    r              = res$r,
    lod_r          = res$lod_r,
    lod_ph         = res$lod_ph,
    logLik         = res$logLik,
    mom_phase_list = res$mom_phase_list,
    lod_ph_list    = res$lod_ph_list,
    q_list         = res$q_list,
    no_linkage     = res$no_linkage,
    optimizer      = res$optimizer,
    n_grid         = res$n_grid,
    markers        = markers,
    n_dams         = length(G_list)
  )
  if (isTRUE(return_diagnostics)) out$diagnostics <- res$diagnostics
  class(out) <- "cpp_pairwise_rf"
  out
}





#' Two-Point Recombination Fraction Estimation for HSMap.data Objects
#'
#' @description
#' This function estimates **pairwise recombination fractions (RFs)** for all
#' marker pairs by aggregating genotype data across multiple nuclear families
#' (dams). It leverages HSMap’s efficient C++ backend for parallel computation.
#'
#' The function first selects a common set of markers based on a user-defined
#' **presence threshold**, which controls whether to use the union, intersection,
#' or a subset of markers present across families. It then prepares the data
#' and invokes the backend to estimate RFs, calculate LOD scores, and infer
#' maternal linkage phase.
#'
#' @section Data Model:
#' The estimation process assumes the following data structure:
#' \itemize{
#'   \item Data consists of \eqn{G} independent nuclear families, each with a known
#'     \strong{maternal genotype vector} \eqn{M_g} and an \strong{offspring
#'     genotype matrix} \eqn{G_g}.
#'   \item Genotypes are coded as 0 (aa), 1 (Aa), 2 (AA), or \code{NA} (missing).
#'   \item In each \eqn{G_g} matrix, markers are in columns and offspring are in rows.
#'   \item Paternal contributions are modeled as a mixture of genotypes and are
#'     statistically integrated out using a robust likelihood function.
#'   \item A family is informative for a marker pair \eqn{(i,j)} only if the dam
#'     is heterozygous for both markers (\emph{double-heterozygous}).
#' }
#'
#' @section Estimation Process:
#' For each marker pair \eqn{(i,j)}, the function:
#' \enumerate{
#'   \item Aggregates the \eqn{3 \times 3} contingency tables of offspring
#'     genotype combinations \eqn{(Y_i, Y_j)} across all families.
#'   \item Estimates paternal allele frequencies from the marginal offspring
#'     genotype distributions.
#'   \item Computes the log-likelihood of the data under both coupling and
#'     repulsion phases for a given recombination fraction \eqn{r}.
#'   \item Maximizes the combined objective function over \eqn{r \in (10^{-6}, 0.49)}
#'     using a golden-section search to find the optimal \eqn{\hat{r}}:
#'     \deqn{ \mathrm{Obj}(r) = \sum_{g=1}^{G} \max\{\ell_g(\mathrm{Coupling}; r), \ell_g(\mathrm{Repulsion}; r)\} }
#' }
#'
#' @section Marker Selection and Presence Threshold:
#' A marker \eqn{t} is considered "present" in family \eqn{g} if its maternal
#' genotype is non-missing and at least one offspring has a non-missing genotype
#' at that marker. The set of markers used for analysis is determined by the
#' `presence_threshold`, \eqn{\tau}.
#'
#' A marker \eqn{t} is retained if its presence frequency, \eqn{f(t)}, meets the
#' threshold:
#' \deqn{ f(t) = \frac{1}{G}\sum_{g=1}^{G} \mathbf{1}\{t \text{ is present in family } g\} \ge \tau }
#' Key threshold values:
#' \itemize{
#'   \item \eqn{\tau = 1}: Strict **intersection** of markers present in all families.
#'   \item \eqn{\tau = 0}: Full **union** of markers present in at least one family.
#'   \item Intermediate \eqn{\tau} (e.g., 0.8) balances marker coverage with data quality.
#' }
#'
#' @section Output Statistics:
#' The function's output includes several \eqn{T \times T} matrices, where \eqn{T}
#' is the number of markers analyzed:
#' \itemize{
#'   \item \code{r}: The matrix of estimated pairwise recombination fractions, \eqn{\hat{r}}.
#'   \item \code{lod_r}: The LOD score comparing the likelihood of linkage at \eqn{\hat{r}}
#'     versus no linkage (\eqn{r = 0.5}).
#'     \deqn{ \mathrm{LOD}_r = \frac{\ell(\hat r) - \ell(0.5)}{\log(10)} }
#'   \item \code{lod_ph}: The LOD score for the inferred maternal phase (coupling vs. repulsion).
#'   \item \code{logLik}: The maximized log-likelihood value, \eqn{\ell(\hat r)}, for each pair.
#'   \item \code{mom_phase_list}: A list of matrices (one per family) indicating the
#'     inferred phase: \code{1} for coupling, \code{0} for repulsion, and
#'     \code{NA} if uninformative.
#' }
#'
#' @section Regularization (`lambda`):
#' The likelihood calculation includes a regularization parameter `lambda` that
#' adds pseudo-counts to the paternal genotype mixture probabilities. This acts
#' as a Dirichlet-type shrinkage, stabilizing estimates for pairs with few
#' informative offspring. The default value of `20` provides a balance between
#' bias and variance.
#'
#' @section Performance and Parallelism:
#' \itemize{
#'   \item **Complexity**: The algorithm scales as \eqn{O(T^2 N)}, where \eqn{T} is
#'     the number of markers and \eqn{N} is the total number of offspring.
#'   \item **Memory**: The primary memory usage comes from the output matrices,
#'     requiring approximately \eqn{32 \cdot T^2} bytes. For \eqn{T=4000}, this
#'     is about 640 MB.
#'   \item **Parallelism**: Computation is parallelized over marker pairs using the
#'     \pkg{RcppParallel} library. The number of threads can be set via the
#'     `threads` argument.
#' }
#'
#' @param x An object of class `HSMap.data`. See \code{\link{read_HSMap_data}}
#'   for details on creating this object.
#' @param snps An optional character vector of marker IDs to include. The final
#'   marker set will be the intersection of this vector and the set determined
#'   by `presence_threshold`.
#' @param presence_threshold A numeric value between 0 and 1 specifying the minimum
#'   proportion of families a marker must be present in to be included. For
#'   convenience, values between 1 and 100 are automatically scaled by 1/100.
#'   Defaults to `1` (strict intersection).
#' @param parent_label An optional character vector for labeling families. If `NULL`
#'   (the default), names are taken from `names(x$G_list)`. If those are also
#'   null, generic labels (`"Pop1"`, `"Pop2"`, etc.) are created.
#' @param threads An integer specifying the number of parallel threads to use.
#'   If `NULL`, \pkg{RcppParallel} uses its default setting.
#' @param lambda A numeric scalar for the pseudo-count weight used in `q`
#'   regularization. Default `20`, **retained for backward compatibility only**;
#'   this value has **not** been selected through a formal sensitivity study and may
#'   change once one is available. See the \emph{Regularization} section.
#' @param q0 Numeric prior target (in `[0, 1]`) for the dam-specific paternal
#'   gametic frequencies \eqn{q_k^{(d)}}: the value each `q` is shrunk toward, via
#'   pseudocounts `alpha = lambda * q0` and `beta = lambda * (1 - q0)`. Default `0.5`.
#' @param tol Numeric tolerance for the r optimizer's bounded local refinement.
#'   Default is `1e-6`.
#' @param maxit The maximum number of iterations for the local refinement.
#'   Default is `200`.
#' @param tiny A small numeric floor to prevent `log(0)` errors in the likelihood
#'   calculation. Default is `1e-12`.
#' @param return_diagnostics A logical value. If `TRUE`, the returned `fit`
#'   includes a `diagnostics` list of per-pair optimizer and missing-data count
#'   matrices. Default `FALSE` (compact output).
#' @param return_input A logical value. If `TRUE`, the aligned input lists
#'   (`G_list` and `M_list`) sent to the C++ backend are included in the output,
#'   which is useful for debugging.
#' @param ... Reserved for deprecated/removed arguments. Passing `r_start` or
#'   `share_pi_across_dams` triggers a deprecation warning (they no longer have an
#'   effect); any other unknown argument name raises an error.
#'
#' @return An object of class `"HSMap.tpt"`, which is a list containing:
#' \itemize{
#'   \item \code{fit}: A list containing the core result matrices. See the
#'     \emph{Output Statistics} section for details.
#'   \item \code{markers}: A character vector of the marker IDs included in the
#'     analysis, in the order they appear in the output matrices.
#'   \item \code{presence_threshold}: The numeric presence threshold used.
#'   \item \code{parent_label}: The character vector of family labels used.
#'   \item \code{threads}: The number of threads used for computation.
#'   \item \code{time_sec}: The total wall-clock time for the analysis in seconds.
#'   \item \code{inputs} (optional): If `return_input = TRUE`, this list contains
#'     the aligned `G_list` and `M_list`.
#' }
#'
#' @section Troubleshooting:
#' \itemize{
#'   \item \strong{All \code{NA} results for a pair}: This typically means no families
#'     were informative for that marker pair (e.g., no double-heterozygous dams).
#'     Try lowering `presence_threshold` to include more data.
#'   \item \strong{High memory usage or slow performance}: The \eqn{T \times T} matrices
#'     can become very large. Reduce the number of markers by filtering, analyzing
#'     chromosomes separately, or using the `snps` argument.
#'   \item \strong{Unexpected marker order}: The order of markers in the output
#'     matrices is given by the `$markers` component of the returned object.
#' }
#'
#' @examples
#' \dontrun{
#' # Load example data
#' dat <- read_HSMap_data(pedigree = "ped.csv", genotypes = "geno.csv")
#'
#' # 1. Analyze with the strict intersection of markers across all families
#' tpt_intersect <- pairwise_rf(dat, presence_threshold = 1, threads = 8)
#' image(tpt_intersect$fit$r, main = "Pairwise RF (Marker Intersection)")
#'
#' # 2. Analyze with markers present in at least 70% of families
#' tpt_union <- pairwise_rf(dat, presence_threshold = 0.7, threads = 8)
#'
#' # 3. Analyze a specific subset of markers and save inputs for debugging
#' my_snps <- tpt_union$markers[1:500]
#' tpt_subset <- pairwise_rf(
#'   dat,
#'   snps = my_snps,
#'   parent_label = paste0("Dam", seq_along(dat$G_list)),
#'   lambda = 30,
#'   threads = 4,
#'   return_input = TRUE
#' )
#'
#' # Inspect the structure of the results
#' str(tpt_subset$fit)
#' }
#'
#' @export
pairwise_rf <- function(
    x,
    snps = NULL,
    presence_threshold = 1,
    parent_label = NULL,
    threads = NULL,
    lambda = 20,
    q0 = 0.5,
    tol = 1e-6,
    maxit = 200,
    tiny = 1e-12,
    return_diagnostics = FALSE,
    return_input = FALSE,
    ...
){
  if (!inherits(x, "HSMap.data")) stop("`x` must be an HSMap.data object.")
  # Deprecated / removed arguments: warn rather than silently ignore.
  .dots <- list(...)
  if (length(.dots)) {
    dep <- intersect(names(.dots), c("r_start", "share_pi_across_dams"))
    for (nm in dep) {
      if (identical(nm, "share_pi_across_dams"))
        warning("`share_pi_across_dams` is deprecated and ignored here; the ",
                "public default is dam-specific q. Use the internal ",
                "`share_q_across_dams` if pooling is required.", call. = FALSE)
      else
        warning(sprintf("`%s` is deprecated and has no effect; it was removed.", nm),
                call. = FALSE)
    }
    unknown <- setdiff(names(.dots), dep)
    if (length(unknown))
      stop("unused argument(s): ", paste(unknown, collapse = ", "), call. = FALSE)
  }
  # Resolve threads here (NULL is allowed at the public API): use the package option
  # when it is a valid positive integer, else a safe default of 1 (deterministic).
  if (is.null(threads)) {
    opt <- getOption("HSMap.n_threads", 1L)
    threads <- if (is.numeric(opt) && length(opt) == 1L && is.finite(opt) && opt >= 1)
      as.integer(opt) else 1L
  }
  if (!is.numeric(threads) || length(threads) != 1L || !is.finite(threads) || threads < 1)
    stop("`threads` must be NULL or a single positive integer.")
  threads <- as.integer(threads)
  G_list_in <- x$G_list
  M_list_in <- x$M_list
  fam_ids   <- names(G_list_in) %||% paste0("Pop", seq_along(G_list_in))

  # --- parent_label handling ---
  if (is.null(parent_label)) {
    parent_label <- fam_ids
  } else if (length(parent_label) == 1L) {
    parent_label <- rep(parent_label, length(fam_ids))
  } else if (length(parent_label) != length(fam_ids)) {
    stop("`parent_label` must be NULL, length 1, or length equal to number of families.")
  }
  names(G_list_in) <- names(M_list_in) <- parent_label

  # --- normalize threshold ---
  thr <- ifelse(presence_threshold > 1, presence_threshold/100, presence_threshold)
  if (thr < 0 || thr > 1) stop("`presence_threshold` must be in [0,1] or [0,100].")

  # --- collect marker presence per family ---
  present_list <- lapply(seq_along(G_list_in), function(g) {
    Gg <- as.matrix(G_list_in[[g]])
    Mg <- M_list_in[[g]]
    m_present <- names(Mg)[!is.na(Mg)]
    o_present <- colnames(Gg)[colSums(!is.na(Gg)) > 0]
    intersect(m_present, o_present)
  })

  union_all <- sort(unique(unlist(present_list)))
  if (!is.null(snps)) union_all <- intersect(union_all, snps)
  if (!length(union_all)) stop("No markers after initial union / snp filter.")

  freq <- table(factor(unlist(present_list), levels = union_all))
  frac_present <- as.numeric(freq) / length(present_list)
  keep <- union_all[frac_present >= thr]
  if (!length(keep)) stop("No markers pass `presence_threshold`.")

  markers <- keep
  Tm <- length(markers)

  # --- align each family ---
  align_one <- function(Gg, Mg) {
    Mv <- rep(NA_integer_, Tm); names(Mv) <- markers
    commonM <- intersect(names(Mg), markers)
    if (length(commonM)) Mv[commonM] <- as.integer(Mg[commonM])

    Gout <- matrix(NA_integer_, nrow = nrow(Gg), ncol = Tm,
                   dimnames = list(rownames(Gg), markers))
    commonG <- intersect(colnames(Gg), markers)
    if (length(commonG)) Gout[, commonG] <- Gg[, commonG, drop = FALSE]
    storage.mode(Gout) <- "integer"
    list(G = Gout, M = Mv)
  }

  G_list <- list(); M_list <- list()
  for (g in seq_along(G_list_in)) {
    aligned <- align_one(G_list_in[[g]], M_list_in[[g]])
    G_list[[g]] <- aligned$G
    M_list[[g]] <- aligned$M
  }
  names(G_list) <- names(M_list) <- parent_label

  # --- run C++ (thread options are set inside cpp_pairwise_rf) ---
  t0 <- proc.time()[["elapsed"]]
  fit <- cpp_pairwise_rf(
    G_list = G_list, M_list = M_list,
    lambda = lambda, q0 = q0, tol = tol, maxit = maxit,
    tiny = tiny, return_diagnostics = return_diagnostics, n_threads = threads
  )
  t1 <- proc.time()[["elapsed"]]

  out <- list(
    fit = fit,
    markers = markers,
    presence_threshold = thr,
    parent_label = parent_label,
    threads = as.integer(threads),
    time_sec = unname(t1 - t0)
  )
  if (isTRUE(return_input)) out$inputs <- list(G_list = G_list, M_list = M_list)
  class(out) <- "HSMap.tpt"
  out
}

# tiny helper
`%||%` <- function(a, b) if (is.null(a)) b else a

