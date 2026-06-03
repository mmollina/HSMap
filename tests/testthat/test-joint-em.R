# Tests for the joint multi-dam EM (shared recombination map).
# Framing: the joint EM is the likelihood-based estimator implied by the model;
# the consensus (per-dam offspring-weighted average) is a diagnostic that behaves
# like shrinkage toward r_start.

make_dat <- function(sim) {
  structure(list(G_list = sim$G_list, M_list = sim$M_list), class = "HSMap.data")
}
oracle_phased <- function(markers, phase_vec, dam) {
  structure(list(dam = dam, order = markers,
                 clusters = integer(length(markers)),
                 phase_vec = as.integer(phase_vec)),
            class = "HSMap.phased")
}
oracle_multi <- function(sim, markers) {
  res <- lapply(seq_along(sim$G_list), function(g)
    oracle_phased(markers, 1L - sim$truth$v_true[[g]], names(sim$G_list)[g]))
  names(res) <- names(sim$G_list)
  class(res) <- "HSMap.phased.multi"
  res
}

test_that("joint EM with a single dam reproduces the single-dam path", {
  RcppParallel::setThreadOptions(numThreads = 1)
  set.seed(11); Tm <- 40
  r_true <- runif(Tm - 1, 0.01, 0.15)
  sim <- sim_multi_pop(T_markers = Tm, n_pops = 1, n_ind_per_pop = 150,
                       marker_intersection = 1, r_vec = r_true,
                       phase_mode = "random", repulsion_rate = 0.3,
                       maternal_geno_mode = "all_het",
                       paternal_pA_base = 0.4, error_rate = 0.01, seed = 11)
  dat <- make_dat(sim); mk <- sim$truth$markers_union
  oph <- oracle_phased(mk, 1L - sim$truth$v_true[[1]], "P1")

  m_single <- hmm_map(dat, phased = oph, dam = 1, epsilon = 0.01,
                      paternal_mode = "per_marker", tol = 1e-6, maxit = 500)
  m_joint  <- hmm_map_joint(dat, phased = oph, dam = "all", epsilon = 0.01,
                            paternal_mode = "per_marker", tol = 1e-6, maxit = 500)

  expect_equal(as.numeric(m_joint$fit$r), as.numeric(m_single$fit$r), tolerance = 1e-8)
})

test_that("joint multi-dam returns a shared map usable downstream", {
  RcppParallel::setThreadOptions(numThreads = 1)
  set.seed(22); Tm <- 40
  r_true <- runif(Tm - 1, 0.01, 0.15)
  sim <- sim_multi_pop(T_markers = Tm, n_pops = 3, n_ind_per_pop = c(80, 80, 80),
                       marker_intersection = 1, r_vec = r_true,
                       phase_mode = "random", repulsion_rate = 0.3,
                       maternal_geno_mode = "all_het",
                       paternal_pA_base = 0.4, paternal_pA_sd = 0.05,
                       error_rate = 0.01, seed = 22)
  dat <- make_dat(sim); mk <- sim$truth$markers_union
  oph <- oracle_multi(sim, mk)

  # hmm_map default method = "joint"
  m <- hmm_map(dat, phased = oph, dam = "all", epsilon = 0.01,
               paternal_mode = "per_marker", tol = 1e-6, maxit = 300)

  expect_s3_class(m, "HSMap.map")
  expect_true(inherits(m, "HSMap.map.joint"))
  expect_length(m$fit$r, Tm - 1)
  expect_true(all(m$fit$r >= 0 & m$fit$r <= 0.5))
  expect_equal(length(m$dams), 3L)

  # downstream consumer: cumulative cM positions, one per marker, monotone
  pos <- HSMap:::get_map(m, map.function = "haldane")
  expect_length(pos, Tm)
  expect_true(all(diff(pos) >= 0))

  # total length within 50% of truth (sanity, not an RMSE claim)
  expect_lt(abs(sum(inv_haldane(m$fit$r)) - sum(inv_haldane(r_true))),
            0.5 * sum(inv_haldane(r_true)))
})

test_that("joint is less sensitive to r_start than the consensus (informed intervals)", {
  RcppParallel::setThreadOptions(numThreads = 1)
  set.seed(33); Tm <- 50
  r_true <- runif(Tm - 1, 0.005, 0.12)
  sim <- sim_multi_pop(T_markers = Tm, n_pops = 3, n_ind_per_pop = c(60, 60, 60),
                       marker_intersection = 1, r_vec = r_true,
                       phase_mode = "random", repulsion_rate = 0.3,
                       maternal_geno_mode = "HWE", maternal_pA = 0.5,
                       paternal_pA_base = 0.4, paternal_pA_sd = 0.05,
                       error_rate = 0.01, seed = 33)
  dat <- make_dat(sim); mk <- sim$truth$markers_union
  oph <- oracle_multi(sim, mk)

  # intervals at least one dam can see (het at both flanking markers)
  Mhet <- sapply(sim$M_list, function(m) as.integer(m == 1))
  ndh  <- sapply(seq_len(Tm - 1), function(k) sum(Mhet[k, ] & Mhet[k + 1, ], na.rm = TRUE))
  inf  <- ndh >= 1
  skip_if(sum(inf) < 5, "too few informed intervals in this seed")

  fit_r <- function(method, rs) {
    res <- hmm_map(dat, phased = oph, dam = "all", method = method,
                   epsilon = 0.01, paternal_mode = "per_marker",
                   r_start = rs, tol = 1e-6, maxit = 300)
    if (method == "joint") as.numeric(res$fit$r) else as.numeric(res$consensus$r)
  }
  j_lo <- fit_r("joint", 0.05);     j_hi <- fit_r("joint", 0.25)
  c_lo <- fit_r("consensus", 0.05); c_hi <- fit_r("consensus", 0.25)

  joint_shift <- mean(abs(j_hi[inf] - j_lo[inf]))
  cons_shift  <- mean(abs(c_hi[inf] - c_lo[inf]))
  expect_lt(joint_shift, cons_shift)
})

test_that("calc_haploprob decodes the joint shared map", {
  RcppParallel::setThreadOptions(numThreads = 1)
  set.seed(44); Tm <- 30
  r_true <- runif(Tm - 1, 0.01, 0.15)
  sim <- sim_multi_pop(T_markers = Tm, n_pops = 2, n_ind_per_pop = c(60, 50),
                       marker_intersection = 1, r_vec = r_true,
                       phase_mode = "random", repulsion_rate = 0.3,
                       maternal_geno_mode = "all_het",
                       paternal_pA_base = 0.4, error_rate = 0.01, seed = 44)
  dat <- make_dat(sim); mk <- sim$truth$markers_union
  oph <- oracle_multi(sim, mk)
  m <- hmm_map(dat, phased = oph, dam = "all", epsilon = 0.01,
               paternal_mode = "per_marker", tol = 1e-6, maxit = 200)
  expect_true(inherits(m, "HSMap.map.joint"))

  gp <- calc_haploprob(dat, m)                # default: all dams
  expect_s3_class(gp, "HSMap.gamma.multi")
  expect_equal(length(gp), 2L)

  g1 <- gp[["P1"]]
  expect_s3_class(g1, "HSMap.gamma")
  arr <- g1$gamma                              # [2 haps, T markers, n ind]
  expect_equal(dim(arr), c(2L, Tm, nrow(sim$G_list[["P1"]])))
  expect_true(all(abs(apply(arr, c(2, 3), sum) - 1) < 1e-6))   # haplotype probs sum to 1

  g_one <- calc_haploprob(dat, m, dam = "P2") # single-dam selection
  expect_s3_class(g_one, "HSMap.gamma")
  expect_equal(g_one$dam, "P2")
})
