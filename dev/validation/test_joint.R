## Validate the joint multi-dam EM:
##  A) single-dam equivalence (joint with D=1 ~ existing single-dam hmm_map)
##  B) clean multi-dam recovery (all_het, large families)
##  C) pooling benefit: joint vs per-dam-average (consensus) on SMALL families
## Run: Rscript dev/validation/test_joint.R
suppressMessages(devtools::load_all(".", quiet = TRUE, export_all = TRUE))
RcppParallel::setThreadOptions(numThreads = 1)
source("dev/validation/recovery_lib.R")

oracle_multi <- function(sim, markers) {
  D <- length(sim$G_list)
  res <- lapply(seq_len(D), function(g)
    oracle_phased(markers, 1L - sim$truth$v_true[[g]], dam = names(sim$G_list)[g]))
  names(res) <- names(sim$G_list)
  class(res) <- "HSMap.phased.multi"
  res
}
rmse <- function(a, b) sqrt(mean((a - b)^2))

## ---------- A) single-dam equivalence ----------
cat("== A) single-dam equivalence ==\n")
set.seed(11); Tm <- 60; r_true <- runif(Tm - 1, 0.01, 0.15)
sim1 <- sim_multi_pop(T_markers = Tm, n_pops = 1, n_ind_per_pop = 250,
                      marker_intersection = 1, r_vec = r_true,
                      phase_mode = "random", repulsion_rate = 0.3,
                      maternal_geno_mode = "all_het",
                      paternal_pA_base = 0.4, error_rate = 0.01, seed = 11)
dat1 <- make_hsmap_data(sim1); markers <- sim1$truth$markers_union
oph1 <- oracle_phased(markers, 1L - sim1$truth$v_true[[1]], dam = "P1")
m_single <- hmm_map(dat1, phased = oph1, dam = 1, epsilon = 0.01,
                    paternal_mode = "per_marker", tol = 1e-6, maxit = 1000)
m_joint1 <- hmm_map_joint(dat1, phased = oph1, dam = "all", epsilon = 0.01,
                          paternal_mode = "per_marker", tol = 1e-6, maxit = 1000)
cat(sprintf("max|r_joint - r_single| = %.3e   (joint iters=%d, single iters=%d)\n",
            max(abs(m_joint1$fit$r - m_single$fit$r)), m_joint1$fit$iters, m_single$fit$iters))
cat(sprintf("recomb: true=%.2f single=%.2f joint=%.2f\n",
            sum(r_true), sum(m_single$fit$r), sum(m_joint1$fit$r)))

## ---------- B) clean multi-dam recovery (all_het, large families) ----------
cat("\n== B) multi-dam recovery (3 dams x 150, all_het) ==\n")
set.seed(22); r_true <- runif(Tm - 1, 0.01, 0.15)
simB <- sim_multi_pop(T_markers = Tm, n_pops = 3, n_ind_per_pop = c(150,150,150),
                      marker_intersection = 1, r_vec = r_true,
                      phase_mode = "random", repulsion_rate = 0.3,
                      maternal_geno_mode = "all_het",
                      paternal_pA_base = 0.4, paternal_pA_sd = 0.05,
                      error_rate = 0.01, seed = 22)
datB <- make_hsmap_data(simB); mk <- simB$truth$markers_union
ophB <- oracle_multi(simB, mk)
mjB <- hmm_map_joint(datB, phased = ophB, dam = "all", epsilon = 0.01,
                     paternal_mode = "per_marker", tol = 1e-6, maxit = 500)
cat(sprintf("joint: RMSE=%.4f recomb true/hat=%.2f/%.2f iters=%d conv=%s\n",
            rmse(r_true, mjB$fit$r), sum(r_true), sum(mjB$fit$r),
            mjB$fit$iters, mjB$fit$converged))

## ---------- C) pooling benefit: joint vs per-dam-average ----------
cat("\n== C) pooling benefit (3 dams x 40, all_het, small r) ==\n")
set.seed(33); r_true <- runif(Tm - 1, 0.005, 0.12)
simC <- sim_multi_pop(T_markers = Tm, n_pops = 3, n_ind_per_pop = c(40,40,40),
                      marker_intersection = 1, r_vec = r_true,
                      phase_mode = "random", repulsion_rate = 0.3,
                      maternal_geno_mode = "all_het",
                      paternal_pA_base = 0.4, paternal_pA_sd = 0.05,
                      error_rate = 0.01, seed = 33)
datC <- make_hsmap_data(simC); mk <- simC$truth$markers_union
ophC <- oracle_multi(simC, mk)

mjC  <- hmm_map_joint(datC, phased = ophC, dam = "all", epsilon = 0.01,
                      paternal_mode = "per_marker", tol = 1e-6, maxit = 500)
mcC  <- hmm_map(datC, phased = ophC, dam = "all", epsilon = 0.01,
                paternal_mode = "per_marker", tol = 1e-6, maxit = 500)
r_joint <- mjC$fit$r
r_cons  <- mcC$consensus$r                       # offspring-weighted avg of per-dam MLEs
perdam  <- sapply(mcC$per_dam, function(m) as.numeric(m$fit$r))   # K x 3
bnd <- function(v) sum(v <= 1e-5 | v >= 0.5 - 1e-5)

cat(sprintf("true recomb = %.2f\n", sum(r_true)))
cat(sprintf("JOINT     : RMSE=%.4f  recomb=%.2f  boundary-hits=%d\n",
            rmse(r_true, r_joint), sum(r_joint), bnd(r_joint)))
cat(sprintf("CONSENSUS : RMSE=%.4f  recomb=%.2f  boundary-hits(avg)=%.1f\n",
            rmse(r_true, r_cons), sum(r_cons), mean(apply(perdam, 2, bnd))))
for (g in 1:3)
  cat(sprintf("  per-dam %d: RMSE=%.4f recomb=%.2f boundary-hits=%d\n",
              g, rmse(r_true, perdam[, g]), sum(perdam[, g]), bnd(perdam[, g])))
cat("DONE\n")
