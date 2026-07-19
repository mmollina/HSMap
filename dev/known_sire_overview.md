# Known-sire / full-sib support — development overview

> **DEVELOPMENT STAGE — NOT VALIDATED FOR PUBLICATION.** This documents the core
> oracle-phase full-sib extension on branch `m6-known-sire-core`. It has passed
> analytical (brute-force) and simulation-recovery checks only. The main manuscript's
> published claims are unchanged and must stay unchanged until this extension is
> independently reviewed and simulated at scale. Companion notes:
> `dev/known_sire_design.md` (math spec) and `dev/known_sire_perf.md` (benchmarks).

## Open-pollinated vs. full-sib likelihoods

- **Open-pollinated (OP, existing):** a two-state maternal HMM. The sire is unknown, so
  the paternal allele is marginalized **independently at each marker** at a dam-specific
  gametic frequency `q_k`. There is no paternal linkage and no paternal map.
- **Full-sib (new):** a genotyped sire transmits **linked** paternal haplotypes, so the
  paternal contribution is a hidden inheritance path, not a per-marker frequency. The
  model is a four-state HMM over `(maternal homolog, paternal homolog)`. Its likelihood
  is genuinely different from the OP likelihood on the same offspring (test 12), and it
  identifies a paternal recombination map that the OP model cannot.

## The four-state hidden process

State `Z = (h_m, h_p)`, `h_m,h_p ∈ {1,2}` (which maternal and which paternal homolog was
transmitted); 4 states, uniform `1/4` initial. Transition is the **phase-free Kronecker
product** `T = T_m(r_m) ⊗ T_p(r_p)`. Phase is applied **once**, in the emission, via the
parents' phased homolog allele labels: expected offspring dosage `= Am[h_m] + Ap[h_p]`,
with symmetric error `ε` (`1-ε` on the expected genotype, `ε/2` on each other class, `1`
on a missing offspring call). No paternal `q` enters the full-sib emission.

## Separate maternal and paternal maps

`r_m` and `r_p` are estimated **separately** and never assumed equal (recovered
independently in tests 4/16). `hmm_map_fullsib()` / `hmm_map_mixed()` return both maps
and their cM distances. The paternal map is a genuine **sire recombination map**; it is
**not** a pollen-pool map.

## Known-sire genotype requirements

Full-sib fitting requires the sire (and mother) genotype at every fitted marker so the
phased allele labels are defined. The core **does not invent** a parental allele and
**does not** substitute `q=0.5`: a missing required parent genotype at a fitted marker
**stops with an error** (`on_missing = "error"`, test 15). Per-cross blockwise exclusion
of missing-sire markers is deferred to a later branch.

## Cross vs. parent identity

Data are organized by **cross** (`cross_id = mother × father`), with an explicit
`"__unknown_sire__"` token for OP crosses. Parent genotypes are stored **once per parent
ID** (`parent_genotypes`); crosses reference them. Thus:

- **A mother used with several sires** contributes all her crosses to **one** maternal
  phase and **one** maternal map (test 8); each cross adds its own maternal meioses to
  the shared maternal recombination counts.
- **A sire used across several dams** contributes all his crosses to **one** paternal
  phase and **one** paternal map (test 9).

## Mixed-family estimation

`hmm_map_mixed()` fits datasets containing both OP and full-sib crosses: it estimates one
shared maternal map from **all** crosses (OP two-state maternal counts pool with full-sib
four-state maternal counts, since both count the same physical maternal crossover) and a
separate shared paternal map from the full-sib crosses; OP crosses keep their dam-specific
`q`. The joint likelihood is the sum of the family-specific likelihoods (test 10). With
only OP crosses present it **dispatches to the existing engine**, so OP-only results are
byte-identical to the current implementation (test 11). `known_sire_untyped` crosses are
never silently treated as full-sib: the user must opt into the OP fallback
(`untyped_sire = "open_pollinated"`) or supply the sire genotype.

## Current limitations (this core branch)

- **Oracle phase only.** Parental phase (maternal and paternal) is supplied and held
  fixed; full-sib two-point estimation and automatic sire-phase inference are **not**
  implemented here.
- **No full-sib linkage-group construction / ordering** — the marker order is supplied.
- **Missing sire genotype** at a fitted marker errors (no blockwise exclusion yet).
- **One shared maternal map and one shared paternal map** — no parent-instance-specific
  or region-specific maps, and no test of maternal-vs-paternal map heterogeneity.
- **Serial engine** (no RcppParallel yet); no polyploid support.
- **Not validated for biological use** — recovery checks only, no real-data analysis and
  no publication simulation study.

## Next branch (full-sib pairwise + phase inference)

See the exact task list in the PR description / final report: full-sib two-point
estimation and LOD; automatic maternal **and** paternal phase inference; full-sib linkage
grouping and ordering; blockwise handling of missing-sire markers; parallelization; and a
publication-scale simulation study before any manuscript claim is added.
