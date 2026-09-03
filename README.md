# HSMap

<!-- badges: start -->
[![R-CMD-check](https://github.com/mmollina/HSMap/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/mmollina/HSMap/actions/workflows/R-CMD-check.yaml)
[![Lifecycle: maturing](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://lifecycle.r-lib.org/articles/stages.html#maturing)
<!-- badges: end -->

**HSMap** performs maternal linkage mapping in **diploid half-sib /
open-pollinated families** in which the maternal parent is known and genotyped
while paternal genotypes are unavailable: one dam, many offspring, unobserved
fathers.

From offspring SNP dosages and the dam's genotype, HSMap:

- estimates **pairwise maternal recombination fractions** and phase support for
  every marker pair;
- assembles the dam's **linkage phase** along an ordered marker set;
- **groups and orders** markers into linkage groups;
- models the **paternal population structure** of the pollen cloud
  (per-marker paternal gametic frequencies and local cross-marker dependence);
- fits a **source-aware hidden Markov model** whose hidden state pairs the
  transmitted maternal homolog with a paternal *source* state, so that
  persistent paternal haplotype **sharing** (as arises when some pollen parents
  are related to the dam) is distinguished from maternal recombination rather
  than absorbed into it;
- reconstructs per-offspring **maternal inheritance paths** (posterior homolog
  probabilities); and
- produces **maternal linkage maps**, with no-linkage or unresolved-phase
  intervals reported as gaps rather than as large centimorgan distances.

The sharing states of the source-aware model are a *statistical* description of
shared-haplotype tracts. They are not father identification, not pedigree
reconstruction, and not proof of identity by descent. Setting the sharing entry
rate to zero (`alpha = 0`) recovers the simpler population-mode estimator
exactly.

**Scope.** The supported design is a single family with one known maternal
parent and unknown fathers (the open-pollinated / half-sib case). Pooling
several dams under one shared map and known-sire full-sib crosses are included
in the software as extensions of this core, and are documented in a companion
paper (forthcoming); they are not part of the released method's paper.

## Installation

```r
# install.packages("remotes")
remotes::install_github("mmollina/HSMap")
```

## Minimal example

A small **simulated** open-pollinated family ships with the package, so the
whole example runs without any external data:

```r
library(HSMap)

ped  <- system.file("extdata", "example_pedigree.csv",  package = "HSMap")
geno <- system.file("extdata", "example_genotypes.csv", package = "HSMap")
dat  <- read_HSMap_data(ped, geno)

## two-point analysis and filtering
tpt  <- pairwise_rf(dat)
tptf <- tpt_filter(tpt, diagnostic.plot = FALSE)

## group and order (the example is a single simulated chromosome)
grp <- group_markers(tptf, k = 1, inter = FALSE)
ord <- mds_order(grp, tptf, plot_each = FALSE)[[1]]

## dam phase along the order
ph <- phase_from_pairwise(tptf, order = ord, dam = "all")

## chromosome-wide maternal map with the source-aware HMM
map <- hmm_map_source_aware(dat, phased = ph, dam = 1, epsilon = 0.01)
map            # length, alpha/beta, sharing-mode occupancy, log-likelihood
```

## Main workflow

1. `read_HSMap_data()` — pedigree + genotype CSVs in, `HSMap.data` out.
2. `pairwise_rf()` / `tpt_filter()` — two-point recombination, phase LODs,
   marker filtering.
3. `group_markers()` / `mds_order()` — linkage-group assignment and
   within-group ordering.
4. `phase_from_pairwise()` — the dam's coupling/repulsion configuration along
   the order.
5. `hmm_map_source_aware()` — the production chromosome-wide estimator
   (maternal recombination + paternal sharing tracts, exact EM).
   `hmm_map()` / `hmm_map_blocks()` provide the population-mode special case
   and phase-block-safe fitting.
6. `get_block_map()`, `get_map()`, `plot_map_list()`, `plot_map_lmv()` —
   gap-aware map summaries and plots; `diagnose_map_intervals()`,
   `collapse_framework_blocks()`, `try_marker()`, `ripple_map()` — framework-map
   construction and local order diagnostics.
7. `calc_haploprob()` — per-offspring maternal inheritance probabilities for
   downstream use (e.g. QTL mapping, crossover analysis).

## Documentation

- `vignette("getting-started", package = "HSMap")` walks through the full
  workflow on simulated data.
- Every exported function has a help page (`?hmm_map_source_aware`, etc.).

## Citation

A methodological manuscript describing HSMap is in preparation. Until it is
available, please cite the package itself (`citation("HSMap")`).

## License

MIT © the HSMap authors.
