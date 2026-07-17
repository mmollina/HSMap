# HSMap 0.1.0

First public release candidate.

## Scope

HSMap builds **maternal linkage maps** for **open-pollinated / unknown-sire diploid
half-sib families**: one or more known dams, many offspring, and unknown fathers. The
unknown paternal contribution is integrated out through a per-marker, dam-specific
paternal gametic frequency; several dams are combined by a joint EM that estimates one
shared recombination map while keeping phase and the paternal model dam-specific.

## Stable public API

- **Reading data:** `read_HSMap_data()`.
- **Two-point analysis:** `pairwise_rf()`, `tpt_filter()`, `pairwise_heatmap()`.
- **Grouping and ordering:** `group_markers()`, `mds_order()`.
- **Phasing:** `phase_from_pairwise()`, `plot_phase()`.
- **Multipoint mapping:** `hmm_map()` (single dam), `hmm_map_joint()` (joint multi-dam),
  `hmm_map_blocks()` (blockwise fitting across unresolved phase).
- **Map reporting and plotting:** `get_block_map()`, `plot_block_map()`,
  `plot_map_list()`, and the map-distance functions `haldane()`/`kosambi()`/`morgan()`
  with their `inv_*` inverses.
- **Decoding and diagnostics:** `calc_haploprob()`, `test_map_heterogeneity()`.
- **Simulation and I/O:** `sim_multi_pop()`, `sim_multi_chrom()`, `make_map()`,
  `write_sim_genotypes()`, `write_sim_pedigree()`.
- **Utilities:** `hs_pal()`, `aggregate_matrix()`, `drop_gap_markers()` (superseded for
  gap handling by the blockwise workflow) and `print` methods for the main result
  classes.

## Highlights

- Blockwise multipoint fitting: the joint EM is fitted only within **resolved phase
  blocks**, so unresolved phase never forces an imputed map.
- **Gap-safe map reporting:** intervals with no linkage (recombination fraction at 0.5)
  or unresolved phase are reported as gaps (`NA` distance), never as large centimorgan
  distances; within-block map segments reset after each gap.
- Conditional global-scale heterogeneity test for a shared-map assumption across dams.
- A small **simulated** example dataset ships in `inst/extdata/` so all examples,
  the README, and the vignette run without any private files.

## Experimental (not part of this release)

Known-sire / full-sib support (`hmm_map_fullsib()`, `hmm_map_mixed()`, `sim_fullsib()`,
and helpers) lives on a separate development branch. It is **oracle-phase only** (no
automatic full-sib two-point or parental-phase inference) and **not ready for automatic
real-data mapping**; its API may change without a deprecation cycle.

## Fixes

- `sim_multi_chrom()` no longer fails when `miss_rate` is left at its default (the
  formal default was self-referential).
