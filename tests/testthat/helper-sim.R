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
