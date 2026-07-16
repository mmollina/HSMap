#!/usr/bin/env Rscript
# =============================================================================
# HSMap real-data pilot  (m5-real-data-pilot)
# -----------------------------------------------------------------------------
# A CONTROLLED, REPRODUCIBLE pilot of the HSMap blockwise mapping workflow on a
# real half-sib family. This is an IMPLEMENTATION / SOFTWARE-EXERCISING pilot,
# NOT the final biological analysis. Provisional QC thresholds and nuisance
# parameters are labelled as such throughout and must not be read as
# recommended defaults.
#
# All input paths and knobs come from the environment (or --key=value CLI args);
# there are NO hard-coded user paths. Nothing this script writes belongs in Git:
# point HSMAP_PILOT_OUTPUT at a location outside the repository (or under
# analysis/output/, which is git-ignored).
#
# Environment variables (CLI --key=value overrides them):
#   HSMAP_PEDIGREE      path to pedigree CSV                        (required)
#   HSMAP_GENOTYPES     path to genotype CSV                        (required)
#   HSMAP_PILOT_OUTPUT  output directory (created if missing)       (required)
#   HSMAP_THREADS       integer thread count      (default: safe auto-detect)
#   HSMAP_PILOT_N       backbone size: 300 = smoke, 1500 = pilot    (default 1500)
#   HSMAP_SEED          RNG seed                                    (default 2026)
#   HSMAP_PILOT_STAGE   qc | backbone | smoke | pilot | sensitivity | all
#                                                                   (default all)
#   HSMAP_ANCHORS_G1    optional path to a G1 anchor-marker vector (RDS/txt/csv)
#   HSMAP_ANCHORS_G2    optional path to a G2 anchor-marker vector (RDS/txt/csv)
#
# Usage:
#   Rscript analysis/real_data_pilot.R                       # env-driven
#   Rscript analysis/real_data_pilot.R --HSMAP_PILOT_N=300 --HSMAP_PILOT_STAGE=smoke
# =============================================================================

suppressWarnings(suppressMessages({
  ok <- requireNamespace("HSMap", quietly = TRUE)
}))

# ---- 0. Config -------------------------------------------------------------
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0L || is.na(a[1])) b else a

parse_cli <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  kv <- list()
  for (a in args) {
    m <- regmatches(a, regexec("^--([A-Za-z0-9_]+)=(.*)$", a))[[1]]
    if (length(m) == 3L) kv[[m[2]]] <- m[3]
  }
  kv
}
CLI <- parse_cli()
cfg <- function(key, default = NULL) {
  v <- CLI[[key]]
  if (is.null(v)) v <- Sys.getenv(key, unset = NA)
  if (is.na(v) || !nzchar(v)) default else v
}

PED    <- cfg("HSMAP_PEDIGREE")
GENO   <- cfg("HSMAP_GENOTYPES")
OUTDIR <- cfg("HSMAP_PILOT_OUTPUT")
SEED   <- as.integer(cfg("HSMAP_SEED", "2026"))
N_BACK <- as.integer(cfg("HSMAP_PILOT_N", "1500"))
STAGE  <- tolower(cfg("HSMAP_PILOT_STAGE", "all"))
ANCH1  <- cfg("HSMAP_ANCHORS_G1")
ANCH2  <- cfg("HSMAP_ANCHORS_G2")

safe_threads <- function() {
  t <- cfg("HSMAP_THREADS")
  if (!is.null(t) && grepl("^[0-9]+$", t) && as.integer(t) >= 1L) return(as.integer(t))
  nc <- tryCatch(parallel::detectCores(logical = TRUE), error = function(e) NA_integer_)
  if (is.na(nc) || nc < 2L) return(1L)
  as.integer(max(1L, min(4L, nc - 1L)))            # conservative, deterministic cap
}
THREADS <- safe_threads()

if (is.null(PED) || is.null(GENO) || is.null(OUTDIR))
  stop("HSMAP_PEDIGREE, HSMAP_GENOTYPES, and HSMAP_PILOT_OUTPUT must all be set ",
       "(env vars or --KEY=value). See analysis/README.md.")
for (p in c(PED, GENO)) if (!file.exists(p)) stop("input not found: ", p)
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# provisional pilot parameters (NOT recommended defaults)
PILOT <- list(
  lambda = 20,        # compatibility-only pseudocount total
  epsilon = 0.05,     # provisional error rate
  gap_r = 0.499,
  tol = 1e-6,
  maxit = 1000L,
  group_k = 12,
  phase_lods_summary = c(0, 1, 3, 5, 10),   # phase sensitivity grid (summaries)
  block_lods = c(3, 5),                       # thresholds used for blockwise fits
  qc = list(mother_het = 1L, min_call_rate = 0.95,
            min_hom_calls = 20L, max_abs_z = 8)   # provisional backbone QC rule
)

set.seed(SEED)
options(warn = 1)             # surface warnings immediately (pilot diagnostics)
LOGCON <- file(file.path(OUTDIR, "pilot_log.txt"), open = "at")
log <- function(...) {
  msg <- sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(..., collapse = ""))
  cat(msg, "\n"); cat(msg, "\n", file = LOGCON); flush(LOGCON)
}
on.exit(close(LOGCON), add = TRUE)

# ---- 1. Provenance ---------------------------------------------------------
git_rev <- tryCatch(system("git rev-parse HEAD", intern = TRUE, ignore.stderr = TRUE),
                    error = function(e) NA_character_)
git_branch <- tryCatch(system("git rev-parse --abbrev-ref HEAD", intern = TRUE, ignore.stderr = TRUE),
                       error = function(e) NA_character_)
pkg_ver <- tryCatch(as.character(utils::packageVersion("HSMap")), error = function(e) NA_character_)
prov <- c(
  "HSMap real-data pilot — provenance",
  paste0("timestamp:        ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste0("git_commit:       ", git_rev[1] %||% NA),
  paste0("git_branch:       ", git_branch[1] %||% NA),
  paste0("HSMap_version:    ", pkg_ver),
  paste0("R_version:        ", R.version.string),
  paste0("seed:             ", SEED),
  paste0("threads:          ", THREADS),
  paste0("backbone_N:       ", N_BACK),
  paste0("stage:            ", STAGE),
  paste0("pedigree:         ", PED),
  paste0("genotypes:        ", GENO),
  paste0("output:           ", OUTDIR),
  paste0("cli_args:         ", paste(commandArgs(trailingOnly = TRUE), collapse = " ")),
  "",
  "PILOT parameters (provisional; NOT recommended defaults):",
  paste0("  ", names(unlist(PILOT[c('lambda','epsilon','gap_r','tol','maxit','group_k')])), " = ",
         unlist(PILOT[c('lambda','epsilon','gap_r','tol','maxit','group_k')])),
  paste0("  phase_lods_summary = ", paste(PILOT$phase_lods_summary, collapse = ",")),
  paste0("  block_lods = ", paste(PILOT$block_lods, collapse = ",")),
  paste0("  backbone QC rule = mother_het==1, call_rate>=", PILOT$qc$min_call_rate,
         ", hom_calls>=", PILOT$qc$min_hom_calls, ", |z|<=", PILOT$qc$max_abs_z)
)
writeLines(prov, file.path(OUTDIR, "provenance.txt"))
writeLines(capture.output(utils::sessionInfo()), file.path(OUTDIR, "sessionInfo.txt"))
log("provenance written; git=", substr(git_rev[1] %||% "NA", 1, 8),
    " pkg=", pkg_ver, " threads=", THREADS, " N=", N_BACK, " stage=", STAGE)

library(HSMap)
options(HSMap.n_threads = THREADS)
if (requireNamespace("RcppParallel", quietly = TRUE))
  RcppParallel::setThreadOptions(numThreads = THREADS)

saveRDS_z <- function(obj, name) {
  f <- file.path(OUTDIR, paste0(name, ".rds"))
  saveRDS(obj, f, compress = "xz"); f
}

# ---- 2. Load + independent raw inspection ----------------------------------
log("reading data via read_HSMap_data() ...")
dat <- read_HSMap_data(PED, GENO)
dam_names <- names(dat$G_list)
log("families: ", length(dam_names), " | markers: ", nrow(dat$alleles),
    " | offspring(fam1): ", nrow(dat$G_list[[1]]))

# independent raw read (Phase 2 cross-check)
ped_raw  <- utils::read.csv(PED, stringsAsFactors = FALSE)
geno_raw <- utils::read.csv(GENO, stringsAsFactors = FALSE, check.names = FALSE)
meta_cols <- intersect(c("marker_id","REF","ALT","chrom","position"), names(geno_raw))
sample_cols <- setdiff(names(geno_raw), meta_cols)
Xms <- as.matrix(geno_raw[, sample_cols, drop = FALSE])   # markers x samples
storage.mode(Xms) <- "integer"
rownames(Xms) <- geno_raw$marker_id

# =============================================================================
# PHASE 2 — QC
# =============================================================================
run_qc <- function() {
  log("PHASE 2: QC")
  # ---- sample-level -------------------------------------------------------
  ped_ids   <- ped_raw$id
  gsamples  <- sample_cols
  dup_ped   <- ped_ids[duplicated(ped_ids)]
  dup_geno  <- gsamples[duplicated(gsamples)]
  ped_not_in_geno <- setdiff(unique(ped_ids), gsamples)
  geno_not_in_ped <- setdiff(gsamples, unique(ped_ids))

  role <- with(ped_raw, ifelse(generation == 1, "mother",
                        ifelse(generation == 2, "offspring", "other")))
  role_map <- setNames(role, ped_ids); fam_map <- setNames(ped_raw$family_id, ped_ids)
  samp_miss <- colMeans(is.na(Xms))
  sample_qc <- data.frame(
    sample = gsamples,
    in_pedigree = gsamples %in% ped_ids,
    role = unname(role_map[gsamples]),
    family_id = unname(fam_map[gsamples]),
    n_markers = nrow(Xms),
    n_called = colSums(!is.na(Xms))[gsamples],
    missingness = round(samp_miss[gsamples], 5),
    stringsAsFactors = FALSE
  )
  utils::write.csv(sample_qc, file.path(OUTDIR, "sample_qc.csv"), row.names = FALSE)

  # ---- marker-level (family 1; single-dam dataset) -----------------------
  fam <- dam_names[1]
  G <- dat$G_list[[fam]]                     # offspring x markers
  M <- dat$M_list[[fam]]                     # named length-markers
  mk <- colnames(G)
  n_off <- nrow(G)
  n0 <- colSums(G == 0L, na.rm = TRUE); n1 <- colSums(G == 1L, na.rm = TRUE)
  n2 <- colSums(G == 2L, na.rm = TRUE); nm <- colSums(is.na(G))
  ncall <- n0 + n1 + n2
  call_rate <- ncall / n_off
  Mv <- M[mk]
  mother_het <- Mv == 1L
  # Mendelian-incompatible calls (maternal-homozygous markers only)
  incompat <- rep(NA_integer_, length(mk))
  incompat[Mv == 0L] <- n2[Mv == 0L]         # mother aa -> offspring cannot be AA
  incompat[Mv == 2L] <- n0[Mv == 2L]         # mother AA -> offspring cannot be aa
  incompat[Mv == 1L] <- 0L                    # het mother -> all offspring genotypes allowed
  # het-mother diagnostics
  het_frac <- ifelse(mother_het & ncall > 0, n1 / ncall, NA_real_)
  # P(offspring het | het mother) = 0.5 exactly (marginal over paternal q); z of Binom(ncall,0.5)
  z_het <- ifelse(mother_het & ncall > 0,
                  (n1 - 0.5 * ncall) / (0.5 * sqrt(ncall)), NA_real_)
  # raw paternal-q from offspring homozygotes (identifiable for het mother, n2+n0>0)
  hom <- n2 + n0
  q_raw <- ifelse(mother_het & hom > 0, n2 / hom, NA_real_)
  # exact duplicate offspring-genotype profile (incl. NA pattern)
  prof <- vapply(seq_len(ncol(G)), function(j) {
    v <- G[, j]; paste0(ifelse(is.na(v), ".", as.character(v)), collapse = "")
  }, character(1))
  is_dup <- duplicated(prof) | duplicated(prof, fromLast = TRUE)

  marker_qc <- data.frame(
    marker_id = mk,
    chrom = dat$alleles$chrom[match(mk, dat$alleles$marker_id)],
    position = dat$alleles$position[match(mk, dat$alleles$marker_id)],
    maternal_geno = Mv,
    mother_het = mother_het,
    n0 = n0, n1 = n1, n2 = n2, n_miss = nm,
    call_rate = round(call_rate, 5),
    mendelian_incompat = incompat,
    het_frac = round(het_frac, 5),
    z_het_dev = round(z_het, 4),
    q_raw = round(q_raw, 5),
    n_hom = hom,
    dup_profile = is_dup,
    profile_hash = prof,
    stringsAsFactors = FALSE
  )
  utils::write.csv(marker_qc[, setdiff(names(marker_qc), "profile_hash")],
                   file.path(OUTDIR, "marker_qc.csv"), row.names = FALSE)
  saveRDS_z(marker_qc, "marker_qc")

  # ---- summary ------------------------------------------------------------
  s <- c(
    "HSMap pilot QC summary",
    paste0("families: ", length(dam_names), "  markers: ", length(mk),
           "  offspring: ", n_off),
    "",
    "== Samples ==",
    paste0("genotype samples: ", length(gsamples)),
    paste0("pedigree ids: ", length(unique(ped_ids))),
    paste0("pedigree ids missing from genotypes: ", length(ped_not_in_geno),
           if (length(ped_not_in_geno)) paste0(" [", paste(utils::head(ped_not_in_geno,10), collapse=","), "]") else ""),
    paste0("genotype samples absent from pedigree: ", length(geno_not_in_ped),
           if (length(geno_not_in_ped)) paste0(" [", paste(utils::head(geno_not_in_ped,10), collapse=","), "]") else ""),
    paste0("duplicate pedigree ids: ", length(dup_ped)),
    paste0("duplicate genotype sample ids: ", length(dup_geno)),
    paste0("sample missingness: min=", round(min(samp_miss),4),
           " median=", round(stats::median(samp_miss),4),
           " max=", round(max(samp_miss),4)),
    "",
    "== Markers ==",
    paste0("maternal genotype: AA(2)=", sum(Mv==2, na.rm=TRUE),
           " Aa(1)=", sum(Mv==1, na.rm=TRUE),
           " aa(0)=", sum(Mv==0, na.rm=TRUE),
           " NA=", sum(is.na(Mv))),
    paste0("maternal heterozygous markers: ", sum(mother_het, na.rm=TRUE)),
    paste0("offspring call rate: min=", round(min(call_rate),3),
           " median=", round(stats::median(call_rate),3),
           " max=", round(max(call_rate),3)),
    paste0("markers with call_rate>=0.95: ", sum(call_rate >= 0.95)),
    paste0("Mendelian-incompatible calls (hom-mother markers): total=",
           sum(incompat, na.rm=TRUE),
           " markers_with>=1=", sum(incompat > 0, na.rm=TRUE)),
    paste0("het-mother het_frac: median=", round(stats::median(het_frac, na.rm=TRUE),3)),
    paste0("|z_het_dev| among het-mother markers: median=",
           round(stats::median(abs(z_het), na.rm=TRUE),2),
           " q95=", round(stats::quantile(abs(z_het), 0.95, na.rm=TRUE),2),
           " max=", round(max(abs(z_het), na.rm=TRUE),2)),
    paste0("q_raw identifiable markers: ", sum(!is.na(q_raw))),
    paste0("markers in an exact-duplicate profile group: ", sum(is_dup),
           " (unique profiles: ", length(unique(prof)), ")"),
    "",
    "NOTE: no markers are deleted here; QC flags are recorded for downstream",
    "backbone selection (Phase 3). Thresholds used later are PILOT-only."
  )
  writeLines(s, file.path(OUTDIR, "qc_summary.txt"))
  log(paste(utils::tail(s, 1)))

  # ---- diagnostic plots ---------------------------------------------------
  grDevices::pdf(file.path(OUTDIR, "qc_plots.pdf"), width = 9, height = 6.5)
  op <- graphics::par(mfrow = c(2, 3), mar = c(4,4,2,1))
  graphics::hist(samp_miss, breaks = 40, col = "grey80",
                 main = "Sample missingness", xlab = "fraction NA")
  graphics::hist(call_rate, breaks = 60, col = "grey80",
                 main = "Marker call rate", xlab = "call rate")
  graphics::barplot(c(AA=sum(Mv==2,na.rm=TRUE), Aa=sum(Mv==1,na.rm=TRUE),
                      aa=sum(Mv==0,na.rm=TRUE), NA.=sum(is.na(Mv))),
                    main = "Maternal genotype", col = "steelblue")
  graphics::hist(het_frac, breaks = 50, col = "grey80",
                 main = "Offspring het frac (het mothers)", xlab = "het fraction")
  graphics::abline(v = 0.5, col = "red", lwd = 2)
  zz <- z_het[is.finite(z_het)]
  graphics::hist(pmax(pmin(zz, 30), -30), breaks = 80, col = "grey80",
                 main = "z of het count vs 0.5", xlab = "standardized deviation (clipped +/-30)")
  graphics::abline(v = c(-8, 8), col = "red", lty = 2)
  graphics::hist(q_raw, breaks = 50, col = "grey80",
                 main = "Raw paternal q (het mothers)", xlab = "q_raw")
  graphics::par(op); grDevices::dev.off()
  log("QC plots -> qc_plots.pdf")

  invisible(list(marker_qc = marker_qc, sample_qc = sample_qc,
                 ped_not_in_geno = ped_not_in_geno, geno_not_in_ped = geno_not_in_ped))
}

# =============================================================================
# PHASE 3 — deterministic pilot backbone
# =============================================================================
load_anchor_vec <- function(path) {
  if (is.null(path) || !file.exists(path)) return(NULL)
  ext <- tolower(tools::file_ext(path))
  out <- tryCatch({
    if (ext == "rds") { v <- readRDS(path); as.character(unlist(v, use.names = FALSE)) }
    else if (ext == "csv") { d <- utils::read.csv(path, stringsAsFactors = FALSE); as.character(d[[1]]) }
    else as.character(readLines(path))
  }, error = function(e) NULL)
  out
}

build_backbone <- function(n_target) {
  log("PHASE 3: backbone (N=", n_target, ")")
  mqc_f <- file.path(OUTDIR, "marker_qc.rds")
  mqc <- if (file.exists(mqc_f)) readRDS(mqc_f) else run_qc()$marker_qc
  q <- PILOT$qc

  # provisional QC rule (pilot-only) with recorded reasons
  reason <- rep("", nrow(mqc))
  fail_het  <- !(mqc$maternal_geno == q$mother_het) | is.na(mqc$maternal_geno)
  fail_call <- mqc$call_rate < q$min_call_rate
  fail_hom  <- mqc$n_hom < q$min_hom_calls | is.na(mqc$n_hom)
  fail_z    <- !is.finite(mqc$z_het_dev) | abs(mqc$z_het_dev) > q$max_abs_z
  add <- function(reason, mask, tag) ifelse(mask, ifelse(nzchar(reason), paste(reason, tag, sep=";"), tag), reason)
  reason <- add(reason, fail_het,  "mother_not_het")
  reason <- add(reason, fail_call, "call_rate<0.95")
  reason <- add(reason, fail_hom,  "hom_calls<20")
  reason <- add(reason, fail_z,    "abs_z>8")
  eligible <- !nzchar(reason)
  log("eligible after provisional QC rule: ", sum(eligible), " / ", nrow(mqc))

  # exact-duplicate binning: deterministic representative = first by (chrom, position, marker_id)
  elig_idx <- which(eligible)
  ord_key <- order(mqc$chrom[elig_idx], mqc$position[elig_idx], mqc$marker_id[elig_idx])
  elig_idx <- elig_idx[ord_key]
  hashes <- mqc$profile_hash[elig_idx]
  rep_mask <- !duplicated(hashes)
  bin_id <- match(hashes, hashes[rep_mask])       # bin index into representatives
  bins <- data.frame(
    representative = mqc$marker_id[elig_idx][bin_id][ifelse(TRUE, seq_along(elig_idx), NA)],
    marker_id = mqc$marker_id[elig_idx],
    bin = bin_id,
    is_representative = rep_mask,
    stringsAsFactors = FALSE
  )
  bins$representative <- mqc$marker_id[elig_idx][match(bins$bin, bin_id)]
  # only multi-member bins are informative to report, but save all for traceability
  utils::write.csv(bins, file.path(OUTDIR, "marker_bins.csv"), row.names = FALSE)
  reps <- mqc$marker_id[elig_idx][rep_mask]
  n_multi <- sum(table(bin_id) > 1)
  reason[match(mqc$marker_id[elig_idx][!rep_mask], mqc$marker_id)] <-
    add(reason[match(mqc$marker_id[elig_idx][!rep_mask], mqc$marker_id)], TRUE, "duplicate_of_representative")
  log("post-dup representatives: ", length(reps), " (", n_multi, " multi-marker bins)")

  # anchors (optional, non-fatal)
  anchors <- unique(c(load_anchor_vec(ANCH1), load_anchor_vec(ANCH2)))
  anchor_in_reps <- intersect(anchors, reps)
  if (length(anchors)) {
    log("anchors supplied: ", length(anchors),
        " | overlap with eligible representatives: ", length(anchor_in_reps),
        " | overlap with all markers: ", length(intersect(anchors, mqc$marker_id)))
  } else log("no anchor vectors supplied (skipping anchor forcing)")

  # stratified deterministic sampling across deciles of q_raw
  rep_qc <- mqc[match(reps, mqc$marker_id), ]
  qv <- rep_qc$q_raw
  brks <- stats::quantile(qv, probs = seq(0, 1, 0.1), na.rm = TRUE, type = 7)
  brks[1] <- -Inf; brks[length(brks)] <- Inf
  decile <- cut(qv, breaks = unique(brks), include.lowest = TRUE, labels = FALSE)
  decile[is.na(decile)] <- 0L
  n_need <- min(n_target, length(reps))
  # force anchors first
  chosen <- character(0)
  forced <- intersect(anchor_in_reps, reps)
  chosen <- c(chosen, forced)
  remaining <- setdiff(reps, chosen)
  # proportional allocation across deciles, deterministic within decile by q then id
  rem_dec <- decile[match(remaining, reps)]
  per_dec <- max(1L, floor((n_need - length(chosen)) / length(unique(rem_dec[rem_dec>0]))))
  set.seed(SEED)
  for (d in sort(unique(rem_dec))) {
    pool <- remaining[rem_dec == d]
    pool <- pool[order(rep_qc$q_raw[match(pool, reps)], pool)]
    take <- min(per_dec, length(pool))
    if (take > 0) chosen <- c(chosen, sample(pool, take))
  }
  # top up (or trim) to exactly n_need, deterministically
  if (length(chosen) < n_need) {
    pool <- setdiff(reps, chosen)
    pool <- pool[order(rep_qc$position[match(pool, reps)], pool)]
    chosen <- c(chosen, utils::head(pool, n_need - length(chosen)))
  }
  chosen <- unique(chosen)[seq_len(min(n_need, length(unique(chosen))))]
  # order backbone by genomic position for reporting (mapping ignores this)
  chosen <- chosen[order(rep_qc$chrom[match(chosen, reps)],
                         rep_qc$position[match(chosen, reps)], chosen)]

  sel <- data.frame(
    marker_id = chosen,
    chrom = mqc$chrom[match(chosen, mqc$marker_id)],
    position = mqc$position[match(chosen, mqc$marker_id)],
    q_raw = mqc$q_raw[match(chosen, mqc$marker_id)],
    decile = decile[match(chosen, reps)],
    is_anchor = chosen %in% forced,
    stringsAsFactors = FALSE
  )
  utils::write.csv(sel, file.path(OUTDIR, "backbone_selected.csv"), row.names = FALSE)
  excl <- data.frame(marker_id = mqc$marker_id, reason = reason,
                     eligible = eligible, stringsAsFactors = FALSE)
  utils::write.csv(excl[nzchar(excl$reason), ],
                   file.path(OUTDIR, "backbone_exclusions.csv"), row.names = FALSE)
  log("backbone selected: ", nrow(sel), " markers (", sum(sel$is_anchor), " forced anchors)")
  saveRDS_z(list(selected = chosen, table = sel, bins = bins,
                 anchors = anchors, decile_counts = table(sel$decile)),
            paste0("backbone_N", n_target))
  invisible(chosen)
}

# =============================================================================
# PHASE 4 / 5 — workflow (gap-safe blockwise)
# =============================================================================
run_workflow <- function(mode = c("smoke", "pilot"), n_target) {
  mode <- match.arg(mode)
  log("PHASE ", if (mode=="smoke") "4 (smoke)" else "5 (pilot)", ": workflow on N=", n_target)
  bb_f <- file.path(OUTDIR, paste0("backbone_N", n_target, ".rds"))
  backbone <- if (file.exists(bb_f)) readRDS(bb_f)$selected else build_backbone(n_target)
  backbone <- intersect(backbone, colnames(dat$G_list[[1]]))
  log("backbone markers present in data: ", length(backbone))

  gc(reset = TRUE); t0 <- Sys.time()
  log("pairwise_rf() on ", length(backbone), " markers, threads=", THREADS, " ...")
  tpt <- pairwise_rf(dat, snps = backbone, threads = THREADS, lambda = PILOT$lambda)
  t_pw <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  mem_pw <- sum(gc()[, "max used"] * c(56, 8)) / 1e6   # rough Mb (Ncells*56 + Vcells*8)
  log("pairwise done in ", round(t_pw,1), "s (~", round(mem_pw), " Mb peak)")

  grp <- group_markers(tpt, k = PILOT$group_k, inter = FALSE)
  gvec <- grp$groups.snp %||% grp$groups
  lg_sizes <- sort(table(gvec), decreasing = TRUE)
  log("group_markers(k=", PILOT$group_k, "): ", length(lg_sizes), " groups; sizes(top): ",
      paste(utils::head(lg_sizes, 8), collapse=","))

  ord <- mds_order(grp, tpt, plot_each = FALSE)
  # ord is a list per LG (or object); normalize to named list of ordered marker vectors
  lg_orders <- .extract_orders(ord, gvec)
  log("mds_order: ordered ", length(lg_orders), " linkage groups")

  # ---- phase-threshold summaries (Phase 5) --------------------------------
  phase_grid <- if (mode == "pilot") PILOT$phase_lods_summary else c(0, 3)
  phase_summ <- list()
  for (lg in names(lg_orders)) {
    o <- lg_orders[[lg]]
    if (length(o) < 3) next
    for (thr in phase_grid) {
      ph <- tryCatch(phase_from_pairwise(tpt, order = o, dam = dam_names[1],
                                         min_phase_lod = thr),
                     error = function(e) { log("  phase LG", lg, " thr", thr, " ERROR: ", conditionMessage(e)); NULL })
      if (is.null(ph)) next
      pv <- ph$phase_vec
      res_via <- ph$resolved_via
      phase_summ[[length(phase_summ)+1]] <- data.frame(
        lg = lg, min_phase_lod = thr, n_markers = length(o),
        n_components = ph$n_components %||% NA_integer_,
        comp_sizes = paste(ph$component_sizes %||% NA, collapse="/"),
        n_direct = sum(res_via == "direct", na.rm=TRUE),
        n_path = sum(res_via == "path", na.rm=TRUE),
        n_unresolved = sum(res_via == "unresolved", na.rm=TRUE),
        prop_resolved = round(mean(!is.na(pv)), 4),
        objective = round(ph$objective %||% NA_real_, 3),
        converged = isTRUE(ph$converged),
        stringsAsFactors = FALSE
      )
    }
  }
  phase_tab <- if (length(phase_summ)) do.call(rbind, phase_summ) else data.frame()
  if (nrow(phase_tab)) utils::write.csv(phase_tab, file.path(OUTDIR, paste0(mode, "_phase_threshold.csv")), row.names = FALSE)

  # ---- blockwise fits at block_lods (Phase 5: 3 & 5) ----------------------
  block_lods <- if (mode == "smoke") 3 else PILOT$block_lods
  block_results <- list(); block_tabs <- list()
  for (thr in block_lods) {
    log("blockwise fits at min_phase_lod=", thr, " ...")
    per_lg <- list()
    for (lg in names(lg_orders)) {
      o <- lg_orders[[lg]]
      if (length(o) < 3) next
      ph <- tryCatch(phase_from_pairwise(tpt, order = o, dam = dam_names[1], min_phase_lod = thr),
                     error = function(e) NULL)
      if (is.null(ph)) next
      gc(reset = TRUE); tb0 <- Sys.time()
      blk <- tryCatch(hmm_map_blocks(dat, ph, epsilon = PILOT$epsilon, lambda = PILOT$lambda,
                                     tol = PILOT$tol, maxit = PILOT$maxit, gap_r = PILOT$gap_r),
                      error = function(e) { log("  block LG", lg, " ERROR: ", conditionMessage(e)); NULL })
      if (is.null(blk)) next
      bm <- tryCatch(get_block_map(blk, "haldane", gap_r = PILOT$gap_r),
                     error = function(e) { log("  get_block_map LG", lg, " ERROR: ", conditionMessage(e)); NULL })
      dt <- as.numeric(difftime(Sys.time(), tb0, units = "secs"))
      per_lg[[lg]] <- list(lg = lg, phase = ph, blocks = blk, block_map = bm, secs = dt)
      block_tabs[[length(block_tabs)+1]] <- .block_report(lg, thr, blk, bm)
    }
    block_results[[paste0("lod", thr)]] <- per_lg
  }
  block_tab <- if (length(block_tabs)) do.call(rbind, block_tabs) else data.frame()
  if (nrow(block_tab)) utils::write.csv(block_tab, file.path(OUTDIR, paste0(mode, "_block_report.csv")), row.names = FALSE)

  # ---- one plot_block_map for the largest fitted LG (gap-safe) ------------
  .plot_largest(block_results, mode)

  t_all <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  timing <- data.frame(mode = mode, n_markers = length(backbone),
                       n_groups = length(lg_sizes), threads = THREADS,
                       pairwise_secs = round(t_pw,1), total_secs = round(t_all,1),
                       pairwise_peak_mb = round(mem_pw))
  utils::write.csv(timing, file.path(OUTDIR, paste0(mode, "_timing.csv")), row.names = FALSE)
  saveRDS_z(list(tpt = tpt, groups = grp, orders = lg_orders,
                 phase_tab = phase_tab, block_results = block_results,
                 block_tab = block_tab, timing = timing), paste0(mode, "_results"))
  log(mode, " workflow complete in ", round(t_all,1), "s; ",
      nrow(block_tab), " fitted blocks reported")
  invisible(list(phase_tab = phase_tab, block_tab = block_tab, timing = timing,
                 orders = lg_orders, tpt = tpt, block_results = block_results))
}

.extract_orders <- function(ord, gvec) {
  # Normalize mds_order output into a named list of ordered marker-id vectors.
  if (is.list(ord) && !is.null(ord$orders)) ord <- ord$orders
  out <- list()
  if (is.list(ord)) {
    for (nm in names(ord)) {
      el <- ord[[nm]]
      v <- if (is.character(el)) el
           else if (is.list(el) && !is.null(el$order)) as.character(el$order)
           else if (is.list(el) && !is.null(el$snp)) as.character(el$snp)
           else if (!is.null(names(el))) names(el) else NULL
      if (!is.null(v) && length(v)) out[[nm]] <- v
    }
  }
  if (!length(out)) {           # fallback: group by gvec, arbitrary order
    sp <- split(names(gvec), gvec)
    for (nm in names(sp)) out[[nm]] <- sp[[nm]]
  }
  out
}

.block_report <- function(lg, thr, blk, bm) {
  blocks <- blk$blocks
  rows <- lapply(seq_along(blocks), function(i) {
    b <- blocks[[i]]
    hm <- b$fit                         # hmm_map object (HSMap.map)
    fit <- if (!is.null(hm)) hm$fit else NULL   # nested numeric fit
    rr <- if (!is.null(fit) && !is.null(fit$r)) as.numeric(fit$r) else numeric(0)
    qq <- if (!is.null(fit) && !is.null(fit$q)) as.numeric(fit$q) else numeric(0)
    data.frame(
      lg = lg, min_phase_lod = thr, block = b$block,
      n_markers = length(b$markers),
      converged = if (!is.null(fit)) isTRUE(fit$converged) else NA,
      conv_reason = if (!is.null(fit)) (fit$conv_reason %||% NA) else NA,
      objective_decreased = if (!is.null(fit)) isTRUE(fit$objective_decreased) else NA,
      iters = if (!is.null(fit)) (fit$iters %||% NA) else NA,
      logLik = if (!is.null(fit)) round(fit$logLik %||% NA_real_, 3) else NA,
      penalized_obj = if (!is.null(fit)) round(fit$penalized_obj %||% NA_real_, 3) else NA,
      r_min = if (length(rr)) round(min(rr, na.rm=TRUE),4) else NA,
      r_med = if (length(rr)) round(stats::median(rr, na.rm=TRUE),4) else NA,
      r_max = if (length(rr)) round(max(rr, na.rm=TRUE),4) else NA,
      q_min = if (length(qq)) round(min(qq, na.rm=TRUE),4) else NA,
      q_med = if (length(qq)) round(stats::median(qq, na.rm=TRUE),4) else NA,
      q_max = if (length(qq)) round(max(qq, na.rm=TRUE),4) else NA,
      n_gaps = if (!is.null(bm)) sum(bm$interval_table$block == b$block &
                                     bm$interval_table$status %in%
                                     c("no_linkage_boundary"), na.rm=TRUE) else NA,
      stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  if (!is.null(bm)) {
    out$total_linked_cM <- round(bm$total_linked_length %||% NA_real_, 2)
    out$n_segments <- bm$n_segments %||% NA
  }
  out
}

.plot_largest <- function(block_results, mode) {
  best <- NULL; best_n <- -1
  for (lod in names(block_results)) for (lg in names(block_results[[lod]])) {
    b <- block_results[[lod]][[lg]]
    nmk <- length(b$phase$order)
    if (!is.null(b$blocks) && nmk > best_n) { best <- b; best_n <- nmk }
  }
  if (is.null(best)) return(invisible())
  p <- tryCatch(plot_block_map(best$blocks, map.function = "haldane"),
                error = function(e) NULL)
  if (!is.null(p) && requireNamespace("ggplot2", quietly = TRUE)) {
    ggplot2::ggsave(file.path(OUTDIR, paste0(mode, "_block_map_largestLG.pdf")),
                    p, width = 9, height = 5)
    log("plot_block_map -> ", mode, "_block_map_largestLG.pdf (LG ", best$lg, ", ", best_n, " markers)")
  }
}

# =============================================================================
# PHASE 6 — nuisance-parameter sensitivity (largest LG only)
# =============================================================================
run_sensitivity <- function(n_target) {
  log("PHASE 6: nuisance sensitivity (largest LG)")
  res_f <- file.path(OUTDIR, "pilot_results.rds")
  if (!file.exists(res_f)) { log("no pilot_results.rds; running pilot first"); run_workflow("pilot", n_target) }
  R <- readRDS(res_f)
  tpt <- R$tpt; orders <- R$orders
  # pick largest LG by markers
  lg <- names(orders)[which.max(vapply(orders, length, integer(1)))]
  o <- orders[[lg]]
  log("largest LG = ", lg, " (", length(o), " markers)")

  lambdas <- c(0, 2, 20); epsilons <- c(0.01, 0.05, 0.10); lods <- c(1, 3, 5, 10)
  rows <- list()
  for (thr in lods) {
    ph <- tryCatch(phase_from_pairwise(tpt, order = o, dam = dam_names[1], min_phase_lod = thr),
                   error = function(e) NULL)
    if (is.null(ph)) next
    ncomp <- ph$n_components %||% NA
    presolved <- round(mean(!is.na(ph$phase_vec)), 4)
    for (lam in lambdas) for (eps in epsilons) {
      blk <- tryCatch(hmm_map_blocks(dat, ph, epsilon = eps, lambda = lam,
                                     tol = PILOT$tol, maxit = PILOT$maxit, gap_r = PILOT$gap_r),
                      error = function(e) NULL)
      bm <- if (!is.null(blk)) tryCatch(get_block_map(blk, "haldane", gap_r = PILOT$gap_r),
                                        error = function(e) NULL) else NULL
      rr <- unlist(lapply(blk$blocks, function(b) if (!is.null(b$fit$fit)) as.numeric(b$fit$fit$r) else NULL))
      qq <- unlist(lapply(blk$blocks, function(b) if (!is.null(b$fit$fit)) as.numeric(b$fit$fit$q) else NULL))
      conv <- all(vapply(blk$blocks, function(b) if (!is.null(b$fit$fit)) isTRUE(b$fit$fit$converged) else TRUE, logical(1)))
      ngap <- if (!is.null(bm)) sum(bm$interval_table$status == "no_linkage_boundary", na.rm=TRUE) else NA
      rows[[length(rows)+1]] <- data.frame(
        lg = lg, min_phase_lod = thr, n_components = ncomp, prop_resolved = presolved,
        lambda = lam, epsilon = eps,
        r_med = if (length(rr)) round(stats::median(rr, na.rm=TRUE),4) else NA,
        map_cM = if (!is.null(bm)) round(bm$total_linked_length %||% NA,2) else NA,
        q_med = if (length(qq)) round(stats::median(qq, na.rm=TRUE),4) else NA,
        converged = conv, n_gaps = ngap, stringsAsFactors = FALSE)
    }
  }
  tab <- if (length(rows)) do.call(rbind, rows) else data.frame()
  if (nrow(tab)) utils::write.csv(tab, file.path(OUTDIR, "sensitivity.csv"), row.names = FALSE)
  saveRDS_z(tab, "sensitivity")
  log("sensitivity complete: ", nrow(tab), " combinations")
  invisible(tab)
}

# ---- dispatch --------------------------------------------------------------
if (STAGE %in% c("qc", "all")) run_qc()
if (STAGE %in% c("backbone", "all")) build_backbone(N_BACK)
if (STAGE == "smoke") run_workflow("smoke", N_BACK)
if (STAGE %in% c("pilot", "all")) run_workflow("pilot", N_BACK)
if (STAGE %in% c("sensitivity", "all")) run_sensitivity(N_BACK)
log("DONE (stage=", STAGE, ")")
