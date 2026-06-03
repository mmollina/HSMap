## ---------------------------------------------------------------------------
## HSMap joint multi-dam validation grid
##
## Run from the package root:
##   Rscript dev/validation/run_joint_grid.R [reps]
##
## Reports bias, RMSE, map length, and phase accuracy for the JOINT (likelihood-
## based) shared-map estimator across multi-dam scenarios. The consensus row is
## included only as a diagnostic (shrinkage-like, r_start-dependent) -- the joint
## EM is NOT presented as the RMSE-minimizing estimator, but as the estimator
## implied by the model (proper likelihood; LRT / heterogeneity test / SEs).
## ---------------------------------------------------------------------------

suppressMessages(devtools::load_all(".", quiet = TRUE, export_all = TRUE))
RcppParallel::setThreadOptions(numThreads = 1)
source("dev/validation/recovery_lib.R")

args <- commandArgs(trailingOnly = TRUE)
reps <- if (length(args) >= 1) as.integer(args[[1]]) else 3L

outdir <- "dev/validation/results"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

baseline <- list(name = "baseline", Tm = 80L, n_dams = 3L, n_ind = 60L,
                 r_lo = 0.005, r_hi = 0.12, rep_rate = 0.3, pA = 0.40, pA_sd = 0.05,
                 err = 0.01, mat_mode = "HWE", r_start = 0.05, base_seed = 1000L)
mk <- function(name, base_seed, ...) modifyList(baseline, c(list(name = name, base_seed = base_seed), list(...)))

scenarios <- list(
  baseline,                                                # 3 dams x 60, HWE
  mk("dams_5",       2000L, n_dams = 5L),                  # more dams
  mk("small_fam",    3000L, n_ind = 30L),                  # small families (pooling matters)
  mk("large_fam",    4000L, n_ind = 150L),                 # large families
  mk("imbalanced",   5000L, n_ind = c(150L, 40L, 40L)),    # one big, two small
  mk("clean_allhet", 6000L, mat_mode = "all_het")          # every interval informed
)

cat(sprintf("== joint validation: %d scenarios x %d reps ==\n", length(scenarios), reps))
t0  <- proc.time()[["elapsed"]]
res <- run_joint_grid(scenarios, reps = reps, verbose = TRUE)
cat(sprintf("== done in %.1fs ==\n\n", proc.time()[["elapsed"]] - t0))

csv <- file.path(outdir, "joint_validation.csv")
utils::write.csv(res, csv, row.names = FALSE)

## ---- summary: means over reps, by scenario x method x phase -----------------
key <- with(res, paste(scenario, method, phase, sep = " | "))
agg <- do.call(rbind, lapply(split(seq_len(nrow(res)), key), function(ix) {
  d <- res[ix, ]
  data.frame(
    scenario   = d$scenario[1], method = d$method[1], phase = d$phase[1],
    informed    = round(mean(d$frac_informed), 2),
    phase_acc   = round(mean(d$phase_acc), 3),
    n_gap       = round(mean(d$n_gap), 2),            # spurious r~0.5 (phase-error) intervals
    RMSE        = round(mean(d$rmse_all), 4),
    RMSE_nogap  = round(mean(d$rmse_nogap), 4),       # excluding the evident phase-error gaps
    bias        = round(mean(d$bias_all), 4),
    recomb_hat  = round(mean(d$recomb_hat), 2),
    recomb_true = round(mean(d$recomb_true), 2),
    cM_capped   = round(mean(d$hald_hat_capped), 0),
    cM_true     = round(mean(d$hald_true), 0),
    stringsAsFactors = FALSE
  )
}))
agg <- agg[order(agg$scenario, agg$method, agg$phase), ]
cat("== summary (means over reps) ==\n")
cat("  joint = likelihood-based estimator; consensus = diagnostic (shrinkage-like)\n\n")
print(agg, row.names = FALSE)
utils::write.csv(agg, file.path(outdir, "joint_validation_summary.csv"), row.names = FALSE)

## ---- plots ------------------------------------------------------------------
png(file.path(outdir, "joint_validation.png"), width = 1200, height = 520, res = 120)
op <- par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))

# (1) calibration of joint (inferred) recombination count
ji <- res[res$method == "joint" & res$phase == "inferred", ]
lim <- range(c(ji$recomb_true, ji$recomb_hat))
plot(ji$recomb_true, ji$recomb_hat, pch = 16, col = "#1f78b4", xlim = lim, ylim = lim,
     xlab = "true expected #recombinations", ylab = "estimated (joint, inferred phase)",
     main = "Joint EM: map-size calibration")
abline(0, 1, lty = 2)

# (2) bias by method (joint vs consensus) -- shows consensus r_start-shrinkage
jb <- res$bias_all[res$method == "joint" & res$phase == "inferred"]
cb <- res$bias_all[res$method == "consensus"]
boxplot(list(`joint\n(inferred)` = jb, `consensus\n(oracle)` = cb),
        ylab = "bias in r", main = "Bias: joint vs consensus", col = c("#1f78b4", "#fdbf6f"))
abline(h = 0, lty = 2)
par(op); dev.off()

cat("\nwrote", csv, "\nwrote", file.path(outdir, "joint_validation_summary.csv"),
    "\nwrote", file.path(outdir, "joint_validation.png"), "\nOK\n")
