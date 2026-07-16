# HSMap real-data pilot (`analysis/`)

This directory holds **scripts and documentation only**. No data, no per-individual
outputs, and no derived haplotypes are tracked here. `analysis/output/` and
`analysis/private/` are git-ignored; keep all pilot outputs there or elsewhere
outside the repository.

The pilot is a **controlled software-exercising run**, not the final biological
analysis. Every QC threshold and nuisance parameter it uses is **provisional** and
labelled as such. Do not read any number it produces as a recommended default or a
biological conclusion.

## What the pilot does

`real_data_pilot.R` runs, in stages:

1. **QC** (Phase 2) — pedigree/sample checks and per-marker diagnostics.
2. **Backbone** (Phase 3) — a deterministic, provisional QC rule; exact-duplicate
   marker binning; stratified deterministic sampling across deciles of the raw
   paternal-`q` estimate; optional anchor forcing.
3. **Smoke test** (Phase 4) — the full gap-safe blockwise workflow on 300 markers.
4. **Pilot** (Phase 5) — the same workflow on the 1,500-marker backbone, with a
   phase-threshold sensitivity table and blockwise fits at `min_phase_lod` 3 and 5.
5. **Sensitivity** (Phase 6) — `lambda`/`epsilon`/`min_phase_lod` sweep on the
   largest fitted linkage group only.

The workflow uses only the **new gap-safe / blockwise** functions:
`pairwise_rf()` → `group_markers(k = 12, inter = FALSE)` → `mds_order(plot_each = FALSE)`
→ `phase_from_pairwise()` (per linkage group) → `hmm_map_blocks()` → `get_block_map()`
→ `plot_block_map()`. It deliberately does **not** use `sum(inv_haldane(r))`, the old
`get_map()`, `drop_gap_markers()` as a substitute for unresolved-block handling, or a
direct `hmm_map()` on unresolved phase.

## Configuration (no hard-coded paths)

All inputs come from environment variables (or `--KEY=value` CLI overrides):

| Variable | Meaning | Default |
|---|---|---|
| `HSMAP_PEDIGREE`     | pedigree CSV path            | (required) |
| `HSMAP_GENOTYPES`    | genotype CSV path            | (required) |
| `HSMAP_PILOT_OUTPUT` | output directory (kept out of Git) | (required) |
| `HSMAP_THREADS`      | thread count                 | safe auto-detect (≤4) |
| `HSMAP_PILOT_N`      | backbone size (300 smoke / 1500 pilot) | 1500 |
| `HSMAP_SEED`         | RNG seed                     | 2026 |
| `HSMAP_PILOT_STAGE`  | `qc`/`backbone`/`smoke`/`pilot`/`sensitivity`/`all` | `all` |
| `HSMAP_ANCHORS_G1`   | optional anchor-marker vector (rds/csv/txt) | — |
| `HSMAP_ANCHORS_G2`   | optional anchor-marker vector (rds/csv/txt) | — |

## Running

```sh
export HSMAP_PEDIGREE=/path/to/ped_HSMap.csv
export HSMAP_GENOTYPES=/path/to/geno_HSMap.csv
export HSMAP_PILOT_OUTPUT=/path/to/output        # keep OUTSIDE Git
export HSMAP_THREADS=4

# smoke test (300 markers)
Rscript analysis/real_data_pilot.R --HSMAP_PILOT_N=300 --HSMAP_PILOT_STAGE=smoke

# staged main pilot (1,500 markers)
Rscript analysis/real_data_pilot.R --HSMAP_PILOT_STAGE=qc
Rscript analysis/real_data_pilot.R --HSMAP_PILOT_STAGE=backbone
Rscript analysis/real_data_pilot.R --HSMAP_PILOT_STAGE=pilot
Rscript analysis/real_data_pilot.R --HSMAP_PILOT_STAGE=sensitivity
```

## Provenance

Every run writes to the output directory: `provenance.txt` (git commit, package
version, CLI args, seed, threads, paths), `sessionInfo.txt`, and `pilot_log.txt`.
Large fitted objects are saved as xz-compressed `.rds` files. None of these belong
in Git.

## Outputs (in `HSMAP_PILOT_OUTPUT`)

`sample_qc.csv`, `marker_qc.csv`, `qc_summary.txt`, `qc_plots.pdf`,
`backbone_selected.csv`, `backbone_exclusions.csv`, `marker_bins.csv`,
`{smoke,pilot}_phase_threshold.csv`, `{smoke,pilot}_block_report.csv`,
`{smoke,pilot}_timing.csv`, `{smoke,pilot}_block_map_largestLG.pdf`,
`sensitivity.csv`, and the corresponding `.rds` bundles.
