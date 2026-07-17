# Generate the small SIMULATED open-pollinated example dataset shipped in
# inst/extdata/ and used by README/examples/tests. Run once from the package root
# (the CSVs are committed; regenerating is only for provenance). Note:
# write_sim_genotypes() reseeds the RNG internally to vary REF/ALT alleles, so the
# REF/ALT letters are not reproducible run-to-run; the 0/1/2 genotype dosages that
# drive the analysis are fully determined by `seed` below.
#
#   Rscript -e "devtools::load_all()" data-raw/make_example_data.R
library(HSMap)

set.seed(2024)
# A single linkage group of 24 markers at ~5 cM (Haldane) spacing. Every dam is
# heterozygous at every marker, so the phase resolves cleanly and the example yields a
# connected map. (No-linkage intervals and unresolved phase ARE handled by the blockwise
# machinery -- reported as gaps with NA distance -- but this small quick-start dataset is
# deliberately clean; see the block/segment tests for the gap-handling paths.)
sim <- sim_multi_pop(
  T_markers          = 24,
  n_pops             = 3,
  n_ind_per_pop      = c(60, 50, 45),
  r_vec              = rep(0.05, 23),            # ~5 cM (Haldane) per interval
  phase_mode         = "random", repulsion_rate = 0.3,
  maternal_geno_mode = "all_het",                # every marker informative in every dam
  paternal_pA_base   = 0.4, paternal_pA_sd = 0.05,
  error_rate         = 0.01, seed = 2024
)

dir.create("inst/extdata", recursive = TRUE, showWarnings = FALSE)
write_sim_pedigree(sim,  file = "inst/extdata/example_pedigree.csv")
write_sim_genotypes(sim, file = "inst/extdata/example_genotypes.csv")
cat("wrote inst/extdata/example_pedigree.csv and example_genotypes.csv\n")
