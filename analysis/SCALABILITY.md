# HSMap scalability assessment and proposed full-data workflow

*Pilot milestone m5. This is an implementation/design note, not a biological result.
All thresholds referenced here are the pilot's provisional values.*

## 1. Measured cost of the current all-pairs workflow

`pairwise_rf()` computes a dense `m × m` two-point result (recombination `r`,
`lod_r`, phase LOD, no-linkage flags, plus per-dam phase matrices), so both time
and memory scale as **O(m²)** in the number of markers `m`.

| backbone `m` | pairwise time | rough peak | groups |
|---:|---:|---:|---:|
| 300  |  0.9 s |  ~330 MB | 12 |
| 1500 | 24.1 s |  ~375 MB | 12 |

Empirical exponent: `log(24.1/0.9) / log(1500/300) ≈ 2.04` → essentially quadratic
(the small excess over 2 is per-pair EM/overhead). Downstream stages are cheap
relative to `pairwise_rf()`: on the 1,500-marker pilot, grouping + MDS ordering +
phase (5 thresholds) + blockwise HMM fits (LOD 3 and 5, 10–12 groups) added
~56 s on top of the 24 s pairwise step (80.5 s total, 4 threads).

## 2. Extrapolation to the full data (NOT executed)

This family has **19,609 markers**; **~10,239 pass the provisional backbone QC**
(maternal-het, call rate ≥ 0.95, ≥ 20 homozygote calls, |z| ≤ 8) after exact-duplicate
binning. Using the fitted `~O(m²)` model and a dense-matrix footprint of
~7 `double` `m × m` matrices (≈ 56 bytes per marker pair):

| target set | `m` | est. pairwise time | est. dense-matrix memory* |
|---|---:|---:|---:|
| pilot backbone | 1,500 | 24 s | ~0.13 GB |
| all QC-passing | 10,239 | **~20 min** | **~5.9 GB** |
| all markers | 19,609 | **~75 min** | **~21.5 GB** |

\* r + lod_r + lod_ph + no_linkage + mom_phase + lod_ph_list + q_list, doubles only;
transient copies during assembly push true peak higher (the 19.6k case realistically
needs 30 GB+).

**The full all-pairs run was deliberately not launched.** At ~10 k markers it is
borderline (~6–8 GB, ~20 min); at ~19.6 k markers the O(m²) memory (~21.5 GB dense,
30 GB+ peak) is not safe on a typical workstation and is not justified for a pilot.
The design below removes the O(m²) wall.

## 3. Proposed scalable full-data workflow (future milestone)

Replace the single global all-pairs step with a **backbone-anchored** pipeline whose
dominant cost is linear in the number of non-backbone markers.

1. **Exact co-segregating bins.** Collapse markers with identical offspring-genotype
   profiles (incl. missing pattern) into bins; carry one representative forward, keep
   the bin membership. (The pilot already computes these: `marker_bins.csv`.)
   Cost O(m) with hashing.
2. **High-quality backbone map.** Select a bounded backbone (`b ≈ 1–2 k` well-behaved
   representatives, stratified as in the pilot), run the existing all-pairs workflow
   on it only: `pairwise_rf → group_markers → mds_order → phase_from_pairwise →
   hmm_map_blocks → get_block_map`. Cost O(b²), which the pilot shows is ~25 s at
   b = 1500.
3. **Assign remaining markers to linkage groups** using pairwise evidence **against
   backbone markers only** — an `(M−b) × b` two-point computation, never the full
   `M × M`. Each non-backbone marker joins the LG of its strongest, phase-consistent
   backbone linkage (with an ambiguity/no-assignment class). Cost O((M−b)·b), linear
   in `M`, and streamable so the full `M × M` matrix is never materialized.
4. **Insertion / local ordering** of assigned markers within their LG: place each new
   marker into the backbone order by local two-point/multipoint likelihood over a
   small window of neighbors; optionally re-run MDS per LG on backbone + inserted
   markers. Cost O(M · w) for a small window `w`, per-LG.
5. **Final blockwise HMM refinement** per LG on the augmented order using the existing
   gap-safe `hmm_map_blocks → get_block_map` (unchanged), which already handles
   unresolved phase and no-linkage gaps.

### Proposed API (sketch — not implemented in the pilot)

```r
# Step 1 (exists in pilot as script logic; promote to a helper)
bin_markers(x)                      # -> list(representatives, bins)

# Step 3
assign_to_backbone(x, backbone_tpt, backbone_groups, snps = remaining,
                   min_lod = 3, threads = NULL)
  # -> data.frame(marker, lg, best_backbone, lod, phase, status)
  #    status in {assigned, ambiguous, unassigned}; streams (M-b)x b, no m x m matrix

# Step 4
insert_markers(x, lg_order, new_markers, window = 10, threads = NULL)
  # -> updated per-LG marker order (+ per-insertion support)
```

### Complexity summary

| stage | current | proposed |
|---|---|---|
| two-point | O(M²) time & memory | O(b²) backbone + O((M−b)·b) streamed assignment |
| memory peak | dense M×M (~21 GB at 19.6 k) | dense b×b (~0.1 GB) + O(b) per streamed marker |
| ordering | global | per-LG local windows |

For `M = 19,609`, `b = 1,500`: assignment is `18,109 × 1,500 ≈ 2.7×10⁷` pairs vs the
full `≈ 1.9×10⁸` — ~7× fewer pairs **and** no `M × M` allocation, so memory drops from
tens of GB to ~0.1 GB resident plus streaming buffers.

### Required tests before adopting this extension

- `bin_markers()`: identical profiles (incl. NA pattern) bin together; a unique profile
  is its own bin; representative choice is deterministic.
- `assign_to_backbone()`: a marker simulated inside a known LG assigns to that LG;
  an unlinked marker is `unassigned`; a marker linked to two LGs is `ambiguous`;
  result never materializes an `M × M` matrix (peak-memory / size assertion); phase
  sign is consistent with the backbone.
- `insert_markers()`: a marker with known map position inserts at the correct interval
  within tolerance; order is unchanged when a marker is already present; window edges
  behave at LG boundaries.
- End-to-end: backbone-anchored map on simulated multi-LG data recovers the same LGs,
  orders (up to reflection), and interval `r` as the all-pairs workflow, within
  tolerance, at materially lower peak memory.
- Determinism: identical results across thread counts (seted RNG for any sampling).

This extension should be added only as a small, independently tested milestone, not
folded into the pilot.
