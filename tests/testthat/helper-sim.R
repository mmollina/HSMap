# Shared helpers for HSMap tests.

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

# HWE genotype-frequency columns (AA, Aa, aa) for allele-A frequency q.
hwe_cols <- function(q) { q <- pmin(pmax(q, 1e-9), 1 - 1e-9); rbind(q^2, 2*q*(1-q), (1-q)^2) }

# Pseudocount penalty added to the observed log-likelihood: sum_k[a log q_k + b log(1-q_k)],
# with a = lambda*q0, b = lambda*(1-q0) (q0 = pseudocount target / prior mode).
pen_q <- function(q, lambda, q0 = 0.5) {
  a <- lambda * q0; b <- lambda * (1 - q0)
  sum(a * log(pmin(pmax(q, 1e-12), 1)) + b * log(pmin(pmax(1 - q, 1e-12), 1)))
}

# One simulated single-dam data set with an oracle phase, returning the pieces the
# engine tests need (data list, oracle phase, and the raw G/M/phase/r_true).
one_dam <- function(seed, Tm = 60L, n = 300L, pA = 0.40, err = 0.01, rr = 0.3, mat = "all_het") {
  set.seed(seed); r_true <- runif(Tm - 1, 0.01, 0.15)
  sim <- sim_multi_pop(T_markers = Tm, n_pops = 1, n_ind_per_pop = n, marker_intersection = 1,
                       r_vec = r_true, phase_mode = "random", repulsion_rate = rr,
                       maternal_geno_mode = mat, maternal_pA = 0.5,
                       paternal_pA_base = pA, error_rate = err, seed = seed)
  mk <- sim$truth$markers_union
  dat <- make_dat(sim)
  oph <- oracle_phased(mk, 1L - sim$truth$v_true[[1]], "P1")
  G <- dat$G_list[[1]][, mk, drop = FALSE]; storage.mode(G) <- "integer"
  list(dat = dat, oph = oph, G = G, M = as.integer(dat$M_list[[1]][mk]),
       ph = as.integer(1L - sim$truth$v_true[[1]]), r_true = r_true, mk = mk)
}
