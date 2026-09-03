# HSMap 0.1.0

First public release candidate.

## Scope

HSMap builds a **maternal linkage map** for a **single open-pollinated / unknown-sire
diploid half-sib family**: one known dam, many offspring, and unknown fathers. The
unknown paternal contribution is integrated out through a per-marker paternal gametic
frequency, estimated jointly with the recombination map by EM. Pooling several dams under
one shared map, and known-sire / full-sib crosses, are extensions developed in a
**companion paper (forthcoming)**; both are included in the software (see *Extensions*
below) but are not part of this release's paper.

## Stable public API

- **Reading data:** `read_HSMap_data()`.
- **Two-point analysis:** `pairwise_rf()`, `tpt_filter()`, `pairwise_heatmap()`.
- **Grouping and ordering:** `group_markers()`, `mds_order()`.
- **Phasing:** `phase_from_pairwise()`, `plot_phase()`.
- **Multipoint mapping:** `hmm_map()` (single family), `hmm_map_blocks()` (blockwise
  fitting across unresolved phase); `hmm_map_joint()` (joint multi-dam — *extension*).
- **Map reporting and plotting:** `get_block_map()`, `plot_block_map()`,
  `plot_map_list()`, and the map-distance functions `haldane()`/`kosambi()`/`morgan()`
  with their `inv_*` inverses.
- **Decoding and diagnostics:** `calc_haploprob()`; `test_map_heterogeneity()`
  (multi-dam — *extension*).
- **Simulation and I/O:** `sim_multi_pop()`, `sim_multi_chrom()`, `make_map()`,
  `write_sim_genotypes()`, `write_sim_pedigree()`.
- **Utilities:** `hs_pal()`, `aggregate_matrix()`, `drop_gap_markers()` (superseded for
  gap handling by the blockwise workflow) and `print` methods for the main result
  classes.

## Highlights

- Blockwise multipoint fitting: the EM is fitted only within **resolved phase blocks**,
  so unresolved phase never forces an imputed map.
- **Gap-safe map reporting:** intervals with no linkage (recombination fraction at 0.5)
  or unresolved phase are reported as gaps (`NA` distance), never as large centimorgan
  distances; within-block map segments reset after each gap.
- A small **simulated** example dataset ships in `inst/extdata/` so all examples,
  the README, and the vignette run without any private files.

## Extensions (companion paper, forthcoming)

These build on the single-family core and are included in the software, but are **not part
of this paper**:

- **Multiple dams under one shared map.** A joint EM pools recombination information across
  families (`hmm_map_joint()`, `hmm_map_blocks()` on several dams), with an optional
  conditional global-scale test for map homogeneity across dams
  (`test_map_heterogeneity()`).
- **Known-sire / full-sib crosses** (`hmm_map_fullsib()`, `hmm_map_mixed()`,
  `sim_fullsib()`, and helpers): a four-state model in which paternal transmission is
  itself a linked hidden path. It is **oracle-phase only** (no automatic full-sib
  two-point or parental-phase inference) and lives on a separate development branch; its
  API may change without a deprecation cycle.

## New features

- Framework-map construction, in two deliberately separated steps:
  `diagnose_map_intervals()` classifies every adjacent interval of a fitted map
  by multipoint/two-point concordance (`normal`, `large_concordant`,
  `large_discordant`, `unresolved`, `no_linkage`; thresholds configurable) and
  never modifies anything; `collapse_framework_blocks()` collapses
  **user-specified** marker blocks to one deterministic representative each
  (largest summed within-block linkage LOD, ties by order), records the other
  block markers as `attached` (never deleted; they inherit the
  representative's position), rebuilds phase and map from scratch, and accepts
  the repair **only** if an automatic 7-point verification passes -- otherwise
  the original map is returned unchanged. Block identification is deliberately
  manual: validation showed automatic block detection is not reproducible
  enough for unattended use, while collapse of correctly identified blocks is
  safe and effective.

- `get_map()` is now exported. It converts a fitted map to gap-aware marker
  positions in cM (no-linkage and unresolved-phase intervals do not become
  huge distances; positions restart after each gap, and `n_gaps`/`block`
  attributes say where they are). It was previously internal, which left
  users hand-rolling `cumsum(inv_haldane(c(0, map$fit$r)))` -- the very
  computation that turns a gap into a spurious distance.

- `plot_map_lmv()`: publication-quality linkage-map figures via
  **LinkageMapView** (`lmv.linkage.plot()`), with marker names or a
  marker-density heat map (`denmap = TRUE`). Positions use the same gap-aware
  rule as `plot_map_list()`; a map containing no-linkage or unresolved-phase
  gaps is drawn as one bar per resolved segment (`split_gaps = TRUE`, default)
  instead of a misleading continuous bar. LinkageMapView is a `Suggests`
  dependency, used conditionally.

## Performance

- Two-point analysis: for a **single family**, `pairwise_rf()` now maximizes each
  phase-specific likelihood exactly, via a safeguarded Newton solve on its
  concave log-likelihood (every two-locus cell probability is affine in *r*),
  keeping the larger phase maximum — mathematically the same objective, measured
  ~10x faster end-to-end (~35x at the optimizer level). Multi-family analyses
  keep the original grid + local-refinement search (the phase decomposition does
  not extend to a sum of per-dam maxima), and the grid implementation remains in
  the package as the reference for regression tests. At `r = 0.5` exactly, the
  phase call is now explicitly `NA` with phase LOD 0, since coupling and
  repulsion coincide at the no-linkage null.

## Fixes

- `tpt_filter()` documentation: `thresh.LOD.ph` is now documented prominently
  as **not** being a phase-LOD filter -- the mask screens the linkage LOD only,
  with the effective cutoff `max(thresh.LOD.ph, thresh.LOD.rf)`. Behaviour is
  unchanged for backward compatibility; the argument's misleading name is
  flagged for a rename with deprecation in a future release.
- `phase_from_pairwise()` is dramatically faster on realistic linkage groups.
  The internal connected-components search labelled vertices when popped and
  rebuilt its stack on every pop, so on a dense phase graph (the normal case
  when `min_phase_lod` admits most marker pairs) the stack grew quadratically
  and the cost grew like the fourth power of the number of markers. Vertices are
  now labelled when pushed and the stack is preallocated. A 400-marker group
  went from 26 s to 0.10 s (260x); component labelling and all phase output are
  unchanged.
- `sim_multi_chrom()` no longer fails when `miss_rate` is left at its default (the
  formal default was self-referential).
