## Headline: does joint pooling beat per-dam-average under REALISTIC (HWE) maternal
## genotypes, where dams are heterozygous at different markers?
## Run: Rscript dev/validation/test_joint_pooling.R [reps]
suppressMessages(devtools::load_all(".", quiet = TRUE, export_all = TRUE))
RcppParallel::setThreadOptions(numThreads = 1)
source("dev/validation/recovery_lib.R")

oracle_multi <- function(sim, markers) {
  res <- lapply(seq_along(sim$G_list), function(g)
    oracle_phased(markers, 1L - sim$truth$v_true[[g]], dam = names(sim$G_list)[g]))
  names(res) <- names(sim$G_list); class(res) <- "HSMap.phased.multi"; res
}
rmse <- function(a, b) sqrt(mean((a - b)^2))

args <- commandArgs(trailingOnly = TRUE)
reps <- if (length(args) >= 1) as.integer(args[[1]]) else 6L
Tm <- 80L

rows <- list()
for (rp in seq_len(reps)) {
  seed <- 300L + rp; set.seed(seed)
  r_true <- runif(Tm - 1, 0.005, 0.12)
  sim <- sim_multi_pop(T_markers = Tm, n_pops = 3, n_ind_per_pop = c(60,60,60),
                       marker_intersection = 1, r_vec = r_true,
                       phase_mode = "random", repulsion_rate = 0.3,
                       maternal_geno_mode = "HWE", maternal_pA = 0.5,
                       paternal_pA_base = 0.4, paternal_pA_sd = 0.05,
                       error_rate = 0.01, seed = seed)
  dat <- make_hsmap_data(sim); mk <- sim$truth$markers_union
  oph <- oracle_multi(sim, mk)

  ## informative dams per interval (het at both flanking markers)
  Mhet <- sapply(sim$M_list, function(m) as.integer(m == 1))
  ndh  <- sapply(1:(Tm-1), function(k) sum(Mhet[k,] & Mhet[k+1,], na.rm = TRUE))

  mj <- hmm_map_joint(dat, phased = oph, dam = "all", epsilon = 0.01,
                      paternal_mode = "gametic", tol = 1e-6, maxit = 500)
  mc <- hmm_map(dat, phased = oph, dam = "all", epsilon = 0.01,
                paternal_mode = "gametic", tol = 1e-6, maxit = 500)
  r_joint <- mj$fit$r
  r_cons  <- mc$consensus$r

  rows[[rp]] <- data.frame(
    rep = rp,
    frac_int_informed = mean(ndh >= 1),
    joint_rmse = rmse(r_true, r_joint),
    cons_rmse  = rmse(r_true, r_cons),
    joint_bias = mean(r_joint - r_true),
    cons_bias  = mean(r_cons  - r_true),
    joint_recomb = sum(r_joint), cons_recomb = sum(r_cons), true_recomb = sum(r_true),
    joint_hald = sum(.haldane_cm(r_joint)), cons_hald = sum(.haldane_cm(r_cons)),
    true_hald  = sum(.haldane_cm(r_true))
  )
  cat(sprintf("rep%d informed=%.2f | RMSE joint=%.4f cons=%.4f | bias joint=%+.4f cons=%+.4f | cM joint=%.0f cons=%.0f true=%.0f\n",
              rp, rows[[rp]]$frac_int_informed, rows[[rp]]$joint_rmse, rows[[rp]]$cons_rmse,
              rows[[rp]]$joint_bias, rows[[rp]]$cons_bias,
              rows[[rp]]$joint_hald, rows[[rp]]$cons_hald, rows[[rp]]$true_hald))
}
res <- do.call(rbind, rows)
cat("\n== means over reps ==\n")
cat(sprintf("RMSE   joint=%.4f  cons=%.4f   (lower is better)\n", mean(res$joint_rmse), mean(res$cons_rmse)))
cat(sprintf("bias   joint=%+.4f cons=%+.4f  (0 is best)\n",      mean(res$joint_bias), mean(res$cons_bias)))
cat(sprintf("length joint=%.0f cons=%.0f true=%.0f cM\n",       mean(res$joint_hald), mean(res$cons_hald), mean(res$true_hald)))
outdir <- "dev/validation/results"; dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
utils::write.csv(res, file.path(outdir, "joint_vs_consensus.csv"), row.names = FALSE)
cat("wrote", file.path(outdir, "joint_vs_consensus.csv"), "\nDONE\n")
