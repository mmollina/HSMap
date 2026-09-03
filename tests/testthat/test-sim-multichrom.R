# Regression: sim_multi_chrom() must run WITHOUT the caller passing `miss_rate`
# (previously the formal default was self-referential, `miss_rate = miss_rate`, which
# recursed with "promise already under evaluation").

test_that("sim_multi_chrom() works with the default miss_rate (no self-referential default)", {
  mp <- make_map(n_chrom = 2, markers_per_chrom = c(5, 5))
  expect_silent(
    sim <- sim_multi_chrom(map = mp, n_pops = 1, n_ind_per_pop = 8,
                           maternal_geno_mode = "all_het", seed = 1)
  )
  expect_true(is.list(sim))
  # explicit miss_rate still honored
  expect_no_error(
    sim_multi_chrom(map = mp, n_pops = 1, n_ind_per_pop = 8,
                    maternal_geno_mode = "all_het", miss_rate = 0.1, seed = 1)
  )
})
