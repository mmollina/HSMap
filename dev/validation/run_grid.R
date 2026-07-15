## ---------------------------------------------------------------------------
## HSMap simulation-recovery driver (pilot grid)
##
## Run from the package root:
##   Rscript dev/validation/run_grid.R [reps]
##
## Builds the package from source (load_all), runs a small grid of scenarios x
## replicates, writes a tidy CSV, prints an aggregated summary, and saves a
## diagnostic PNG. Scale up by editing `scenarios` / `reps`.
## ---------------------------------------------------------------------------

suppressMessages(devtools::load_all(".", quiet = TRUE, export_all = TRUE))
RcppParallel::setThreadOptions(numThreads = 1)
source("dev/validation/recovery_lib.R")

args <- commandArgs(trailingOnly = TRUE)
reps <- if (length(args) >= 1) as.integer(args[[1]]) else 3L

outdir <- "dev/validation/results"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

## Baseline + one-axis variations. base_seed differs per scenario so reps are
## independent across scenarios but reproducible.
baseline <- list(name = "baseline", Tm = 80L, n_ind = 300L, r_lo = 0.01, r_hi = 0.15,
                 rep_rate = 0.25, pA = 0.40, err = 0.01, eps = 0.01,
                 pat = "gametic", base_seed = 1000L)

mk <- function(name, base_seed, ...) modifyList(baseline, c(list(name = name, base_seed = base_seed), list(...)))

scenarios <- list(
  baseline,
  mk("n_ind_100",   2000L, n_ind = 100L),
  mk("n_ind_600",   3000L, n_ind = 600L),
  mk("err_0",       4000L, err = 0.00),
  mk("err_0.05",    5000L, err = 0.05),
  mk("repuls_0.50", 6000L, rep_rate = 0.50),
  mk("dense_lowr",  7000L, r_lo = 0.005, r_hi = 0.05),
  mk("pat_pA_0.1",  8000L, pA = 0.10)
)

cat(sprintf("== running %d scenarios x %d reps ==\n", length(scenarios), reps))
t0  <- proc.time()[["elapsed"]]
res <- run_grid(scenarios, reps = reps, verbose = TRUE)
cat(sprintf("== done in %.1fs ==\n", proc.time()[["elapsed"]] - t0))

csv <- file.path(outdir, "recovery_pilot.csv")
utils::write.csv(res, csv, row.names = FALSE)
cat("wrote", csv, "\n\n")

## ---- aggregated summary (inferred vs oracle, by scenario) ------------------
agg <- function(df) {
  s <- split(df, df$scenario)
  do.call(rbind, lapply(names(s), function(nm) {
    d  <- s[[nm]]
    f  <- function(method, col) mean(d[[col]][d$method == method], na.rm = TRUE)
    data.frame(
      scenario      = nm,
      phase_acc     = round(f("inferred", "phase_acc"), 3),
      n_gap_inf     = round(f("inferred", "n_gap"), 2),
      cor_nogap_inf = round(f("inferred", "cor_nogap"), 3),
      cor_nogap_or  = round(f("oracle",   "cor_nogap"), 3),
      recomb_true   = round(f("oracle",   "recomb_true"), 1),
      recomb_inf    = round(f("inferred", "recomb_hat"), 1),
      recomb_or     = round(f("oracle",   "recomb_hat"), 1),
      hald_raw_inf  = round(f("inferred", "hald_hat_raw"), 0),
      hald_cap_inf  = round(f("inferred", "hald_hat_capped"), 0),
      hald_true     = round(f("oracle",   "hald_true"), 0),
      stringsAsFactors = FALSE
    )
  }))
}
summary_tbl <- agg(res)
cat("== summary (means over reps) ==\n")
print(summary_tbl, row.names = FALSE)
utils::write.csv(summary_tbl, file.path(outdir, "recovery_summary.csv"), row.names = FALSE)

## ---- diagnostic plot: recomb-count recovery (robust) inferred vs oracle ----
png(file.path(outdir, "recovery_recomb.png"), width = 1100, height = 520, res = 120)
op <- par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
inf <- res[res$method == "inferred", ]
orc <- res[res$method == "oracle", ]
lim <- range(c(inf$recomb_true, inf$recomb_hat, orc$recomb_hat), na.rm = TRUE)
plot(orc$recomb_true, orc$recomb_hat, pch = 16, col = "#1f78b4", xlim = lim, ylim = lim,
     xlab = "true expected #recombinations", ylab = "estimated", main = "Oracle phase")
abline(0, 1, lty = 2)
plot(inf$recomb_true, inf$recomb_hat, pch = 16, col = "#e31a1c", xlim = lim, ylim = lim,
     xlab = "true expected #recombinations", ylab = "estimated", main = "Inferred phase")
abline(0, 1, lty = 2)
par(op); dev.off()
cat("wrote", file.path(outdir, "recovery_recomb.png"), "\n")
cat("OK\n")
