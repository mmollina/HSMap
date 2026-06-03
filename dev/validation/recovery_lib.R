## ---------------------------------------------------------------------------
## HSMap simulation-recovery validation library
##
## Reusable functions to quantify how well the pipeline recovers recombination
## (r), phase, and map length from simulated single-dam data. Separates the
## ORACLE (true order + true phase) from the INFERRED pipeline (true order +
## pairwise/graph phase) so phasing error can be told apart from HMM error.
##
## Robust metrics (see Finding O): a single mis-phased interval drives r_hat->0.5,
## and inv_haldane(0.5)~1577 cM, so naive summed map length is fragile. We report
## bounded recombination count sum(r), gap-flagged Haldane length, and metrics on
## non-gap intervals, alongside the raw length.
##
## Usage:
##   source("dev/validation/recovery_lib.R")    # after devtools::load_all(".")
## ---------------------------------------------------------------------------

`%||%` <- function(a, b) if (is.null(a)) b else a

## Haldane cM with the same clamping the package uses
.haldane_cm <- function(r) {
  r <- pmin(pmax(as.numeric(r), 1e-8), 0.499999)
  -50 * log(1 - 2 * r)
}

## Build a minimal HSMap.data from a sim_multi_pop object
make_hsmap_data <- function(sim) {
  structure(list(G_list = sim$G_list, M_list = sim$M_list), class = "HSMap.data")
}

## True HMM phase_vec for dam g: 1 = coupling, 0 = repulsion (= 1 - v_true).
## phase_vec is invariant to a global homolog flip, so it can be compared directly.
phase_vec_truth <- function(sim, g = 1L) as.integer(1L - sim$truth$v_true[[g]])

## A ready-to-use HSMap.phased object carrying a known phase_vec (oracle).
oracle_phased <- function(markers, phase_vec, dam = "P1") {
  structure(list(dam = dam, order = markers,
                 clusters = integer(length(markers)),
                 phase_vec = as.integer(phase_vec)),
            class = "HSMap.phased")
}

## Recombination-recovery metrics. gap_thr flags intervals the HMM pushed toward
## "unlinked" (likely phase errors).
metrics_r <- function(r_true, r_hat, gap_thr = 0.45) {
  stopifnot(length(r_true) == length(r_hat))
  is_gap <- r_hat > gap_thr
  ng     <- sum(!is_gap)
  data.frame(
    n_int          = length(r_true),
    n_gap          = sum(is_gap),                                   # spurious "unlinked" intervals
    max_r_hat      = max(r_hat),
    cor_all        = if (stats::sd(r_true) > 0) stats::cor(r_true, r_hat) else NA_real_,
    cor_nogap      = if (ng > 2 && stats::sd(r_true[!is_gap]) > 0)
                       stats::cor(r_true[!is_gap], r_hat[!is_gap]) else NA_real_,
    rmse_all       = sqrt(mean((r_true - r_hat)^2)),
    rmse_nogap     = if (ng > 0) sqrt(mean((r_true[!is_gap] - r_hat[!is_gap])^2)) else NA_real_,
    bias_all       = mean(r_hat - r_true),
    bias_nogap     = if (ng > 0) mean(r_hat[!is_gap] - r_true[!is_gap]) else NA_real_,
    recomb_true    = sum(r_true),                                   # expected #recombinations (bounded, robust)
    recomb_hat     = sum(r_hat),
    hald_true      = sum(.haldane_cm(r_true)),
    hald_hat_raw   = sum(.haldane_cm(r_hat)),                       # can blow up at r~0.5
    hald_hat_capped= sum(.haldane_cm(pmin(r_hat, gap_thr))),        # robust length (gaps capped)
    stringsAsFactors = FALSE
  )
}

## Run ONE scenario replicate: simulate single-dam data, then fit ORACLE and
## INFERRED. Returns a 2-row data.frame (one per method) of params + metrics.
run_scenario <- function(scen, rep_id) {
  seed <- scen$base_seed + rep_id
  set.seed(seed)

  r_true <- stats::runif(scen$Tm - 1L, scen$r_lo, scen$r_hi)
  sim <- sim_multi_pop(
    T_markers = scen$Tm, n_pops = 1L, n_ind_per_pop = scen$n_ind,
    marker_intersection = 1, r_vec = r_true,
    phase_mode = "random", repulsion_rate = scen$rep_rate,
    maternal_geno_mode = "all_het",
    paternal_pA_base = scen$pA, error_rate = scen$err, seed = seed
  )
  dat     <- make_hsmap_data(sim)
  markers <- sim$truth$markers_union
  ptrue   <- phase_vec_truth(sim, 1L)

  fit_eps <- scen$eps
  pat     <- scen$pat

  ## ORACLE: true order + true phase
  t0  <- proc.time()[["elapsed"]]
  mor <- hmm_map(dat, phased = oracle_phased(markers, ptrue), dam = 1,
                 epsilon = fit_eps, paternal_mode = pat, tol = 1e-6, maxit = 500)
  t_or <- proc.time()[["elapsed"]] - t0

  ## INFERRED: true order + pairwise/graph phase
  t0  <- proc.time()[["elapsed"]]
  tpt <- pairwise_rf(dat, presence_threshold = 1, threads = 1)
  ph  <- phase_from_pairwise(tpt, order = markers, dam = "all")
  min <- hmm_map(dat, phased = ph, dam = 1,
                 epsilon = fit_eps, paternal_mode = pat, tol = 1e-6, maxit = 500)
  t_in <- proc.time()[["elapsed"]] - t0
  pacc <- mean(ph$phase_vec == ptrue)              # flip-invariant; direct compare

  base <- data.frame(
    scenario = scen$name, Tm = scen$Tm, n_ind = scen$n_ind,
    r_lo = scen$r_lo, r_hi = scen$r_hi, rep_rate = scen$rep_rate,
    pA = scen$pA, err = scen$err, eps = fit_eps, pat = pat,
    rep = rep_id, seed = seed, stringsAsFactors = FALSE
  )
  row_or <- cbind(base, method = "oracle",   phase_acc = NA_real_,
                  metrics_r(r_true, as.numeric(mor$fit$r)),
                  iters = mor$fit$iters, secs = round(t_or, 2))
  row_in <- cbind(base, method = "inferred", phase_acc = pacc,
                  metrics_r(r_true, as.numeric(min$fit$r)),
                  iters = min$fit$iters, secs = round(t_in, 2))
  rbind(row_or, row_in)
}

## ---------------------------------------------------------------------------
## Multi-dam (joint EM) mode
## ---------------------------------------------------------------------------

## Oracle per-dam phase (HSMap.phased.multi) from sim truth.
oracle_multi <- function(sim, markers) {
  res <- lapply(seq_along(sim$G_list), function(g)
    oracle_phased(markers, 1L - sim$truth$v_true[[g]], dam = names(sim$G_list)[g]))
  names(res) <- names(sim$G_list)
  class(res) <- "HSMap.phased.multi"
  res
}

## Run ONE multi-dam scenario replicate. Simulates D dams sharing one map with
## dam-specific phase and unknown (heterogeneous) sires, then fits:
##   - joint EM with ORACLE phase    (estimation ceiling)
##   - joint EM with INFERRED phase  (full pipeline; the headline estimator)
##   - consensus with ORACLE phase   (diagnostic; shrinkage-like, r_start-dependent)
## Returns a tidy 3-row data.frame (params + metrics). The joint EM is the
## likelihood-based estimator; the consensus row is reported only for comparison.
run_scenario_joint <- function(scen, rep_id) {
  seed <- scen$base_seed + rep_id
  set.seed(seed)
  D     <- scen$n_dams
  n_ind <- if (length(scen$n_ind) == 1L) rep(scen$n_ind, D) else scen$n_ind
  rs    <- scen$r_start %||% 0.05
  matm  <- scen$mat_mode %||% "HWE"

  r_true <- stats::runif(scen$Tm - 1L, scen$r_lo, scen$r_hi)
  sim <- sim_multi_pop(
    T_markers = scen$Tm, n_pops = D, n_ind_per_pop = n_ind,
    marker_intersection = 1, r_vec = r_true,
    phase_mode = "random", repulsion_rate = scen$rep_rate,
    maternal_geno_mode = matm, maternal_pA = 0.5,
    paternal_pA_base = scen$pA, paternal_pA_sd = scen$pA_sd %||% 0.05,
    error_rate = scen$err, seed = seed
  )
  dat <- make_hsmap_data(sim); mk <- sim$truth$markers_union

  ## interval informativeness (>=1 dam het at both flanking markers)
  Mhet <- sapply(sim$M_list, function(m) as.integer(m == 1))   # Tm x D
  ndh  <- sapply(seq_len(scen$Tm - 1L), function(k) sum(Mhet[k, ] & Mhet[k + 1, ], na.rm = TRUE))
  frac_informed <- mean(ndh >= 1)

  oph <- oracle_multi(sim, mk)
  tpt <- pairwise_rf(dat, presence_threshold = 1, threads = 1)
  iph <- phase_from_pairwise(tpt, order = mk, dam = "all")

  ## phase accuracy: per dam, on that dam's informed intervals; averaged
  pacc <- mean(vapply(seq_len(D), function(g) {
    ptrue_g <- as.integer(1L - sim$truth$v_true[[g]])
    inf_g   <- (Mhet[-scen$Tm, g] & Mhet[-1, g]) == 1
    if (sum(inf_g, na.rm = TRUE) < 1) return(NA_real_)
    mean(iph[[g]]$phase_vec[inf_g] == ptrue_g[inf_g], na.rm = TRUE)
  }, numeric(1)), na.rm = TRUE)

  base <- function(method, phase, pa) data.frame(
    scenario = scen$name, Tm = scen$Tm, n_dams = D,
    n_ind = paste(n_ind, collapse = "/"), mat_mode = matm, r_start = rs,
    rep = rep_id, seed = seed, frac_informed = round(frac_informed, 3),
    method = method, phase = phase, phase_acc = pa, stringsAsFactors = FALSE)

  fit_joint <- function(phased)
    hmm_map(dat, phased = phased, dam = "all", method = "joint", epsilon = 0.01,
            paternal_mode = "per_marker", r_start = rs, tol = 1e-6, maxit = 300)$fit$r
  cons_r <- hmm_map(dat, phased = oph, dam = "all", method = "consensus", epsilon = 0.01,
                    paternal_mode = "per_marker", r_start = rs, tol = 1e-6, maxit = 300)$consensus$r

  rbind(
    cbind(base("joint",     "oracle",   NA_real_), metrics_r(r_true, fit_joint(oph))),
    cbind(base("joint",     "inferred", pacc),     metrics_r(r_true, fit_joint(iph))),
    cbind(base("consensus", "oracle",   NA_real_), metrics_r(r_true, cons_r))
  )
}

## Run a list of multi-dam scenarios x reps.
run_joint_grid <- function(scenarios, reps = 3L, verbose = TRUE) {
  out <- list(); k <- 0L
  for (scen in scenarios) for (rp in seq_len(reps)) {
    k <- k + 1L
    d <- run_scenario_joint(scen, rp)
    out[[k]] <- d
    if (verbose) {
      ji <- d[d$method == "joint" & d$phase == "inferred", ]
      co <- d[d$method == "consensus", ]
      cat(sprintf("[%s rep%d] informed=%.2f phase=%.3f | joint-inf RMSE=%.4f bias=%+.4f recomb=%.1f/%.1f | cons bias=%+.4f\n",
                  scen$name, rp, ji$frac_informed, ji$phase_acc, ji$rmse_all, ji$bias_all,
                  ji$recomb_hat, ji$recomb_true, co$bias_all))
    }
  }
  do.call(rbind, out)
}

## Run a list of scenarios x reps; rbind all rows.
run_grid <- function(scenarios, reps = 5L, verbose = TRUE) {
  out <- list(); k <- 0L
  for (scen in scenarios) {
    for (rp in seq_len(reps)) {
      k <- k + 1L
      out[[k]] <- run_scenario(scen, rp)
      if (verbose) {
        r <- out[[k]]
        cat(sprintf("[%s rep%d] inferred: phase=%.3f n_gap=%d cor_nogap=%.3f recomb hat/true=%.1f/%.1f\n",
                    scen$name, rp,
                    r$phase_acc[r$method == "inferred"],
                    r$n_gap[r$method == "inferred"],
                    r$cor_nogap[r$method == "inferred"],
                    r$recomb_hat[r$method == "inferred"],
                    r$recomb_true[r$method == "inferred"]))
      }
    }
  }
  do.call(rbind, out)
}
