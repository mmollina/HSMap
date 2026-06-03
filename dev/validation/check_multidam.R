## Check: can we simulate the multi-dam case for the joint EM?
##   3 dams x 100 offspring, unknown (heterogeneous) sires, ONE shared map,
##   dam-specific phase. Confirms the sim output is well-formed for joint EM.
## Run: Rscript dev/validation/check_multidam.R
suppressMessages(devtools::load_all(".", quiet = TRUE, export_all = TRUE))
RcppParallel::setThreadOptions(numThreads = 1)
source("dev/validation/recovery_lib.R")
set.seed(123)

Tm <- 60
r_true <- runif(Tm - 1, 0.01, 0.15)             # ONE shared recombination map

sim <- sim_multi_pop(
  T_markers = Tm, n_pops = 3, n_ind_per_pop = c(100, 100, 100),
  marker_intersection = 1, r_vec = r_true,      # all dams share all markers + map
  phase_mode = "random", repulsion_rate = 0.3,  # each dam gets its OWN phase
  maternal_geno_mode = "HWE", maternal_pA = 0.5,# realistic: dams het at different markers
  paternal_pA_base = 0.4, paternal_pA_sd = 0.05,# heterogeneous UNKNOWN sire pools
  known_sire_prob_per_pop = 0, error_rate = 0.01, seed = 123
)

cat("== structure ==\n")
cat("dams:", length(sim$G_list), "| names:", paste(names(sim$G_list), collapse = ","), "\n")
for (g in 1:3) cat(sprintf("  %s: G = %d x %d | M len %d | maternal het rate %.2f\n",
  names(sim$G_list)[g], nrow(sim$G_list[[g]]), ncol(sim$G_list[[g]]),
  length(sim$M_list[[g]]), mean(sim$M_list[[g]] == 1, na.rm = TRUE)))
cat("shared r_true: length", length(sim$truth$r_true),
    "| head", paste(round(head(sim$truth$r_true, 5), 3), collapse = ","), "\n")

cat("\n== unknown sires ==\n")
cat("father_geno all-NA per dam:",
    paste(sapply(sim$father_geno_list, function(f) all(is.na(f))), collapse = ","), "\n")
cat("pi_true[,1] dam1 (AA,Aa,aa):", paste(round(sim$pi_true_list[[1]][, 1], 3), collapse = ","),
    "  <- HWE mixture, not one-hot\n")
cat("pi_true colSums dam1 head:",
    paste(round(head(colSums(sim$pi_true_list[[1]]), 5), 3), collapse = ","), "\n")
pat_p <- sapply(1:3, function(g) sim$pi_true_list[[g]][1, 1] + 0.5 * sim$pi_true_list[[g]][2, 1])
cat("per-dam sire A-freq at m1 (heterogeneous):", paste(round(pat_p, 3), collapse = ","), "\n")

cat("\n== dam-specific phase differs ==\n")
for (g in 1:3) cat(sprintf("  dam %d v_true head: %s | repulsion frac %.2f\n",
  g, paste(head(sim$truth$v_true[[g]], 12), collapse = ""), mean(sim$truth$v_true[[g]])))
agree12 <- mean(sim$truth$v_true[[1]] == sim$truth$v_true[[2]])
cat(sprintf("  phase agreement dam1 vs dam2: %.2f (should be ~0.5 = independent)\n", agree12))

cat("\n== interval informativeness for the SHARED map ==\n")
Mhet <- sapply(sim$M_list, function(m) as.integer(m == 1))      # Tm x 3 het indicator
adj_dh <- sapply(1:(Tm - 1), function(k) sum(Mhet[k, ] & Mhet[k + 1, ], na.rm = TRUE))
cat("adjacent double-het dam count per interval: min", min(adj_dh),
    "| mean", round(mean(adj_dh), 2), "| intervals with 0 informative dams:",
    sum(adj_dh == 0), "of", Tm - 1, "\n")
cat("  (with HWE maternal some intervals rely on only 1-2 dams -> exactly where pooling helps)\n")

cat("\n== upstream multi-dam two-point works ==\n")
dat <- make_hsmap_data(sim)
tpt <- pairwise_rf(dat, presence_threshold = 1, threads = 1)
cat("tpt$fit$r:", paste(dim(tpt$fit$r), collapse = "x"),
    "| mom_phase_list length:", length(tpt$fit$mom_phase_list), "\n")
cat("OK\n")
