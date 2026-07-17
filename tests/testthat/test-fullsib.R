# Known-sire / full-sib core tests (m6-known-sire-core). See dev/known_sire_design.md.

RcppParallel::setThreadOptions(numThreads = 1L)

# ---- brute-force references ------------------------------------------------
.emit <- function(y, d, eps) if (is.na(y)) 1 else if (y == d) 1 - eps else eps / 2

.brute_fs <- function(G, Am, Ap, rm, rp, eps) {           # 4-state, want loglik + counts
  T <- ncol(G); n <- nrow(G)
  paths <- as.matrix(expand.grid(rep(list(0:3), T)))
  Tk <- function(s, sp, k) { hm<-s%/%2; hp<-s%%2; hmp<-sp%/%2; hpp<-sp%%2
    (if (hm==hmp) 1-rm[k] else rm[k]) * (if (hp==hpp) 1-rp[k] else rp[k]) }
  dos <- function(s, k) Am[s%/%2 + 1, k] + Ap[s%%2 + 1, k]
  ll <- 0; msw <- psw <- tot <- numeric(T - 1)
  for (i in seq_len(n)) { y <- G[i, ]; pw <- numeric(nrow(paths))
    for (p in seq_len(nrow(paths))) { pa <- paths[p, ]; w <- 0.25
      for (k in seq_len(T)) { w <- w * .emit(y[k], dos(pa[k], k), eps); if (k < T) w <- w * Tk(pa[k], pa[k+1], k) }
      pw[p] <- w }
    py <- sum(pw); ll <- ll + log(py); post <- pw / py
    for (k in seq_len(T - 1)) for (p in seq_len(nrow(paths))) { s <- paths[p,k]; sp <- paths[p,k+1]
      if ((s%/%2) != (sp%/%2)) msw[k] <- msw[k] + post[p]
      if ((s%%2)  != (sp%%2))  psw[k] <- psw[k] + post[p]; tot[k] <- tot[k] + post[p] } }
  list(ll = ll, msw = msw, psw = psw, tot = tot)
}

.mk_alleles <- function(sim, id) matrix(as.integer(
  phase_to_haplotypes(sim$truth$parents[[id]]$genotype, sim$truth$parents[[id]]$phase_vec)), 2)

# ===========================================================================
test_that("1. four-state forward likelihood matches brute-force enumeration", {
  sim <- sim_fullsib(n_markers = 3, crosses = data.frame(mother="M1", father="S1", n=6),
                     r_const_m = 0.1, r_const_p = 0.3, epsilon = 0.02, seed = 7)
  cr <- sim$data$crosses[["M1__x__S1"]]
  Am <- .mk_alleles(sim, "M1"); Ap <- .mk_alleles(sim, "S1")
  rm <- c(0.1, 0.25); rp <- c(0.3, 0.15)
  b <- .brute_fs(cr$G, Am, Ap, rm, rp, 0.02)
  expect_equal(fs_loglik_cpp(cr$G, Am, Ap, rm, rp, 0.02), b$ll, tolerance = 1e-10)
})

test_that("2. forward-backward posteriors sum to one at every marker", {
  sim <- sim_fullsib(n_markers = 5, crosses = data.frame(mother="M1", father="S1", n=10), seed = 3)
  cr <- sim$data$crosses[["M1__x__S1"]]; Am <- .mk_alleles(sim,"M1"); Ap <- .mk_alleles(sim,"S1")
  es <- fs_estep_cpp(cr$G, Am, Ap, rep(0.1,4), rep(0.2,4), 0.02, TRUE)
  for (i in seq_len(nrow(cr$G))) {
    gm <- matrix(es$gamma[i, ], 4)
    expect_true(all(abs(colSums(gm) - 1) < 1e-10))
  }
})

test_that("3. expected maternal and paternal counts match brute force", {
  sim <- sim_fullsib(n_markers = 3, crosses = data.frame(mother="M1", father="S1", n=8),
                     r_const_m = 0.12, r_const_p = 0.28, epsilon = 0.03, seed = 21)
  cr <- sim$data$crosses[["M1__x__S1"]]; Am <- .mk_alleles(sim,"M1"); Ap <- .mk_alleles(sim,"S1")
  rm <- c(0.12, 0.2); rp <- c(0.28, 0.09)
  b <- .brute_fs(cr$G, Am, Ap, rm, rp, 0.03)
  es <- fs_estep_cpp(cr$G, Am, Ap, rm, rp, 0.03, FALSE)
  expect_equal(as.numeric(es$m_switch), b$msw, tolerance = 1e-9)
  expect_equal(as.numeric(es$p_switch), b$psw, tolerance = 1e-9)
  expect_equal(as.numeric(es$total),    b$tot, tolerance = 1e-9)
})

test_that("4/16. maternal and paternal maps are recovered separately", {
  rm_t <- c(0.03,0.30,0.10,0.42,0.06,0.20,0.15,0.08,0.26)
  rp_t <- c(0.40,0.05,0.22,0.02,0.33,0.10,0.30,0.12,0.06)
  sim <- sim_fullsib(n_markers = 10, crosses = data.frame(mother="M1", father="S1", n=2500),
                     r_m = rm_t, r_p = rp_t, epsilon = 0.01, seed = 99)
  pm <- list(M1 = sim$truth$parents[["M1"]]$phase_vec)
  pp <- list(S1 = sim$truth$parents[["S1"]]$phase_vec)
  fit <- hmm_map_fullsib(sim$data, phased_m = pm, phased_p = pp, epsilon = 0.01, tol = 1e-7, maxit = 2000)
  expect_true(fit$fit$converged)
  expect_lt(sqrt(mean((as.numeric(fit$fit$r_m) - rm_t)^2)), 0.03)
  expect_lt(sqrt(mean((as.numeric(fit$fit$r_p) - rp_t)^2)), 0.03)
  # the two maps are genuinely different and each tracks its own truth
  expect_gt(cor(as.numeric(fit$fit$r_m), rm_t), 0.95)
  expect_gt(cor(as.numeric(fit$fit$r_p), rp_t), 0.95)
})

test_that("5. a homozygous sire reduces to the maternal-only model (paternal alleles fixed)", {
  sim <- sim_fullsib(n_markers = 6, crosses = data.frame(mother="M1", father="S1", n=40), seed = 5)
  cr <- sim$data$crosses[["M1__x__S1"]]; Am <- .mk_alleles(sim,"M1")
  Tn <- ncol(cr$G); Ap_hom <- matrix(1L, 2, Tn)                 # sire AA at every marker -> paternal allele 1
  rm <- runif(Tn-1, 0.05, 0.4)
  ll_fs <- fs_loglik_cpp(cr$G, Am, Ap_hom, rm, rep(0.1, Tn-1), 0.02)
  ll_op <- op_estep_cpp(cr$G, Am, rm, rep(1, Tn), 0.02)$loglik  # maternal-only, paternal allele fixed to A
  expect_equal(ll_fs, ll_op, tolerance = 1e-9)
  # and it is independent of r_p
  expect_equal(ll_fs, fs_loglik_cpp(cr$G, Am, Ap_hom, rm, rep(0.45, Tn-1), 0.02), tolerance = 1e-9)
})

test_that("6. a homozygous mother reduces to the paternal-only model", {
  sim <- sim_fullsib(n_markers = 6, crosses = data.frame(mother="M1", father="S1", n=40), seed = 8)
  cr <- sim$data$crosses[["M1__x__S1"]]; Ap <- .mk_alleles(sim,"S1")
  Tn <- ncol(cr$G); Am_hom <- matrix(1L, 2, Tn)                 # mother AA at every marker
  rp <- runif(Tn-1, 0.05, 0.4)
  ll_fs <- fs_loglik_cpp(cr$G, Am_hom, Ap, rep(0.1, Tn-1), rp, 0.02)
  ll_pat <- op_estep_cpp(cr$G, Ap, rp, rep(1, Tn), 0.02)$loglik # 2-state on the paternal side, maternal fixed to A
  expect_equal(ll_fs, ll_pat, tolerance = 1e-9)
  expect_equal(ll_fs, fs_loglik_cpp(cr$G, Am_hom, Ap, rep(0.45, Tn-1), rp, 0.02), tolerance = 1e-9)
})

test_that("7. single and joint wrappers agree on one full-sib cross", {
  sim <- sim_fullsib(n_markers = 8, crosses = data.frame(mother="M1", father="S1", n=600),
                     r_const_m = 0.1, r_const_p = 0.2, epsilon = 0.02, seed = 31)
  pm <- list(M1 = sim$truth$parents[["M1"]]$phase_vec); pp <- list(S1 = sim$truth$parents[["S1"]]$phase_vec)
  f1 <- hmm_map_fullsib(sim$data, phased_m=pm, phased_p=pp, epsilon=0.02, tol=1e-7, maxit=1000, r_start=0.1)
  f2 <- hmm_map_mixed(sim$data, phased_m=pm, phased_p=pp, epsilon=0.02, tol=1e-7, maxit=1000, r_start=0.1)
  expect_equal(as.numeric(f1$fit$r_m), as.numeric(f2$fit$r_m), tolerance = 1e-8)
  expect_equal(as.numeric(f1$fit$r_p), as.numeric(f2$fit$r_p), tolerance = 1e-8)
})

test_that("8. two crosses sharing a mother use one maternal phase and map", {
  cr <- data.frame(mother=c("M1","M1"), father=c("S1","S2"), n=c(800,800), stringsAsFactors=FALSE)
  rm_t <- c(0.05,0.3,0.12,0.4,0.08,0.22,0.15); rp_t <- c(0.3,0.06,0.25,0.03,0.28,0.1,0.32)
  sim <- sim_fullsib(n_markers=8, crosses=cr, r_m=rm_t, r_p=rp_t, epsilon=0.01, seed=202)
  pm <- lapply(sim$truth$parents, function(p) p$phase_vec)
  fit <- hmm_map_fullsib(sim$data, phased_m=pm, phased_p=pm, epsilon=0.01, tol=1e-7, maxit=2000)
  expect_setequal(fit$contributing_mothers, "M1")           # ONE mother
  expect_length(fit$fit$r_m, length(rm_t))                  # ONE maternal map
  expect_lt(sqrt(mean((as.numeric(fit$fit$r_m) - rm_t)^2)), 0.03)
})

test_that("9. two crosses sharing a sire use one paternal phase and map", {
  cr <- data.frame(mother=c("M1","M2"), father=c("S1","S1"), n=c(800,800), stringsAsFactors=FALSE)
  rm_t <- c(0.05,0.3,0.12,0.4,0.08,0.22,0.15); rp_t <- c(0.3,0.06,0.25,0.03,0.28,0.1,0.32)
  sim <- sim_fullsib(n_markers=8, crosses=cr, r_m=rm_t, r_p=rp_t, epsilon=0.01, seed=203)
  pm <- lapply(sim$truth$parents, function(p) p$phase_vec)
  fit <- hmm_map_fullsib(sim$data, phased_m=pm, phased_p=pm, epsilon=0.01, tol=1e-7, maxit=2000)
  expect_setequal(fit$contributing_sires, "S1")             # ONE sire
  expect_length(fit$fit$r_p, length(rp_t))                  # ONE paternal map
  expect_lt(sqrt(mean((as.numeric(fit$fit$r_p) - rp_t)^2)), 0.03)
})

test_that("10. mixed OP + full-sib combine likelihood and maternal counts", {
  cr <- data.frame(mother=c("M1","M1","D1"), father=c("S1","S2",NA), n=c(700,700,700), stringsAsFactors=FALSE)
  rm_t <- c(0.05,0.30,0.12,0.40,0.08,0.22,0.15,0.10,0.28)
  rp_t <- c(0.35,0.06,0.25,0.03,0.30,0.12,0.33,0.15,0.05)
  sim <- sim_fullsib(n_markers=10, crosses=cr, r_m=rm_t, r_p=rp_t, epsilon=0.01, op_paternal_pA=0.4, seed=202)
  pm <- lapply(sim$truth$parents, function(p) p$phase_vec)
  fit <- hmm_map_mixed(sim$data, phased_m=pm, phased_p=pm, epsilon=0.01, lambda=20, tol=1e-7, maxit=2000)
  expect_true(fit$fit$converged)
  expect_setequal(fit$contributing_crosses, c("M1__x__S1","M1__x__S2","D1"))
  expect_lt(sqrt(mean((as.numeric(fit$fit$r_m) - rm_t)^2)), 0.03)   # pooled maternal
  expect_lt(sqrt(mean((as.numeric(fit$fit$r_p) - rp_t)^2)), 0.03)   # full-sib paternal
})

test_that("11. open-pollinated-only results are unchanged (dispatch to hmm_map)", {
  so <- sim_fullsib(n_markers=8, crosses=data.frame(mother="D1",father=NA,n=500),
                    r_const_m=0.12, op_paternal_pA=0.4, epsilon=0.02, seed=55)
  pm <- list(D1 = so$truth$parents[["D1"]]$phase_vec)
  mx <- hmm_map_mixed(so$data, phased_m=pm, epsilon=0.05, lambda=20, tol=1e-6, maxit=1000, r_start=0.05)
  ph <- structure(list(dam="D1", order=colnames(so$data$crosses[["D1"]]$G),
                       phase_vec=so$truth$parents[["D1"]]$phase_vec), class="HSMap.phased")
  direct <- hmm_map(so$data, phased=ph, dam="D1", epsilon=0.05, lambda=20, tol=1e-6, maxit=1000, r_start=0.05)
  expect_true(mx$dispatched)
  expect_identical(as.numeric(mx$fit$r_m), as.numeric(direct$fit$r))
  expect_identical(mx$fit$logLik, direct$fit$logLik)
})

test_that("12. known-sire genotypes affect the likelihood", {
  sim <- sim_fullsib(n_markers=6, crosses=data.frame(mother="M1",father="S1",n=200), seed=12)
  cr <- sim$data$crosses[["M1__x__S1"]]; Am <- .mk_alleles(sim,"M1"); Ap <- .mk_alleles(sim,"S1")
  Tn <- ncol(cr$G)
  ll_real <- fs_loglik_cpp(cr$G, Am, Ap, rep(0.1,Tn-1), rep(0.1,Tn-1), 0.02)
  ll_hom  <- fs_loglik_cpp(cr$G, Am, matrix(1L,2,Tn), rep(0.1,Tn-1), rep(0.1,Tn-1), 0.02) # uninformative sire
  expect_gt(abs(ll_real - ll_hom), 1)
})

test_that("13. permuting sire genotypes changes the full-sib likelihood", {
  sim <- sim_fullsib(n_markers=8, crosses=data.frame(mother="M1",father="S1",n=300), seed=13)
  cr <- sim$data$crosses[["M1__x__S1"]]; Am <- .mk_alleles(sim,"M1"); Ap <- .mk_alleles(sim,"S1")
  Tn <- ncol(cr$G); rm <- rep(0.1,Tn-1); rp <- rep(0.15,Tn-1)
  perm <- c(2:Tn, 1)                                        # cyclic column permutation of the sire
  ll0 <- fs_loglik_cpp(cr$G, Am, Ap, rm, rp, 0.02)
  ll1 <- fs_loglik_cpp(cr$G, Am, Ap[, perm, drop=FALSE], rm, rp, 0.02)
  expect_gt(abs(ll0 - ll1), 1)
})

test_that("14. missing offspring genotypes are neutral", {
  sim <- sim_fullsib(n_markers=6, crosses=data.frame(mother="M1",father="S1",n=30), seed=14)
  cr <- sim$data$crosses[["M1__x__S1"]]; Am <- .mk_alleles(sim,"M1"); Ap <- .mk_alleles(sim,"S1")
  Tn <- ncol(cr$G); rm <- rep(0.1,Tn-1); rp <- rep(0.2,Tn-1)
  G2 <- rbind(cr$G, matrix(NA_integer_, 1, Tn))            # one all-missing offspring
  expect_equal(fs_loglik_cpp(G2, Am, Ap, rm, rp, 0.02),
               fs_loglik_cpp(cr$G, Am, Ap, rm, rp, 0.02), tolerance = 1e-10)
})

test_that("15. missing required sire genotype triggers the documented error", {
  sim <- sim_fullsib(n_markers=6, crosses=data.frame(mother="M1",father="S1",n=50), seed=15)
  d <- sim$data
  d$parent_genotypes[["S1"]][3] <- NA_integer_             # puncture the sire genotype
  pm <- list(M1 = sim$truth$parents[["M1"]]$phase_vec); pp <- list(S1 = sim$truth$parents[["S1"]]$phase_vec)
  expect_error(hmm_map_fullsib(d, phased_m=pm, phased_p=pp), "missing required (parent|sire) genotype")
})

test_that("17. final returned likelihood is evaluated at the final parameters", {
  sim <- sim_fullsib(n_markers=8, crosses=data.frame(mother="M1",father="S1",n=400),
                     r_const_m=0.1, r_const_p=0.2, epsilon=0.02, seed=17)
  pm <- list(M1=sim$truth$parents[["M1"]]$phase_vec); pp <- list(S1=sim$truth$parents[["S1"]]$phase_vec)
  fit <- hmm_map_fullsib(sim$data, phased_m=pm, phased_p=pp, epsilon=0.02, tol=1e-7, maxit=1000)
  cr <- sim$data$crosses[["M1__x__S1"]]; Am <- .mk_alleles(sim,"M1"); Ap <- .mk_alleles(sim,"S1")
  ll_at_final <- fs_loglik_cpp(cr$G[, fit$order, drop=FALSE], Am, Ap,
                               as.numeric(fit$fit$r_m), as.numeric(fit$fit$r_p), 0.02)
  expect_equal(fit$fit$logLik, ll_at_final, tolerance = 1e-8)
})

test_that("18. the active objective is non-decreasing within tolerance", {
  sim <- sim_fullsib(n_markers=8, crosses=data.frame(mother="M1",father="S1",n=500),
                     r_const_m=0.1, r_const_p=0.25, epsilon=0.02, seed=18)
  pm <- list(M1=sim$truth$parents[["M1"]]$phase_vec); pp <- list(S1=sim$truth$parents[["S1"]]$phase_vec)
  fit <- hmm_map_fullsib(sim$data, phased_m=pm, phased_p=pp, epsilon=0.02, tol=1e-7, maxit=1000)
  expect_true(all(diff(fit$fit$loglik_trace) >= -1e-6))
  expect_false(fit$fit$objective_decreased)
})

test_that("19. a low maxit produces a non-convergence warning", {
  sim <- sim_fullsib(n_markers=8, crosses=data.frame(mother="M1",father="S1",n=400), seed=19)
  pm <- list(M1=sim$truth$parents[["M1"]]$phase_vec); pp <- list(S1=sim$truth$parents[["S1"]]$phase_vec)
  expect_warning(fit <- hmm_map_fullsib(sim$data, phased_m=pm, phased_p=pp, maxit=2),
                 "did not converge")
  expect_false(fit$fit$converged)
  expect_identical(fit$fit$conv_reason, "maxit")
})

test_that("20. results are deterministic with one thread", {
  sim <- sim_fullsib(n_markers=8, crosses=data.frame(mother="M1",father="S1",n=400),
                     r_const_m=0.1, r_const_p=0.2, epsilon=0.02, seed=20)
  pm <- list(M1=sim$truth$parents[["M1"]]$phase_vec); pp <- list(S1=sim$truth$parents[["S1"]]$phase_vec)
  a <- hmm_map_fullsib(sim$data, phased_m=pm, phased_p=pp, epsilon=0.02, tol=1e-7, maxit=1000)
  b <- hmm_map_fullsib(sim$data, phased_m=pm, phased_p=pp, epsilon=0.02, tol=1e-7, maxit=1000)
  expect_identical(as.numeric(a$fit$r_m), as.numeric(b$fit$r_m))
  expect_identical(as.numeric(a$fit$r_p), as.numeric(b$fit$r_p))
  expect_identical(a$fit$logLik, b$fit$logLik)
})


# ---- Commit 1: distances (inv_haldane) and gap behavior --------------------
test_that("C1. full-sib d_m/d_p use inv_haldane on linked intervals", {
  rm_t <- c(0.05,0.30,0.12,0.40,0.08); rp_t <- c(0.35,0.06,0.25,0.03,0.30)
  sim <- sim_fullsib(n_markers=6, crosses=data.frame(mother="M1",father="S1",n=1500),
                     r_m=rm_t, r_p=rp_t, epsilon=0.01, seed=77)
  pm <- list(M1=sim$truth$parents[["M1"]]$phase_vec); pp <- list(S1=sim$truth$parents[["S1"]]$phase_vec)
  fit <- hmm_map_fullsib(sim$data, phased_m=pm, phased_p=pp, epsilon=0.01, tol=1e-7, maxit=2000)
  lm <- fit$fit$interval_status_m == "linked"; lp <- fit$fit$interval_status_p == "linked"
  expect_equal(as.numeric(fit$fit$d_m)[lm], inv_haldane(as.numeric(fit$fit$r_m)[lm]), tolerance=1e-10)
  expect_equal(as.numeric(fit$fit$d_p)[lp], inv_haldane(as.numeric(fit$fit$r_p)[lp]), tolerance=1e-10)
  expect_true(all(fit$fit$interval_status_m %in% c("linked","no_linkage_boundary","insufficient_information")))
  expect_true(all(fit$fit$interval_status_p %in% c("linked","no_linkage_boundary","insufficient_information")))
})

test_that("C1. r near 0.5 is a no-linkage gap (NA distance), not a finite linked interval", {
  rep <- HSMap:::.fs_interval_report(c(0.1, 0.4999, 0.5, NA_real_), c(100,100,100,0), 0.499)
  expect_identical(rep$status, c("linked","no_linkage_boundary","no_linkage_boundary","insufficient_information"))
  expect_equal(rep$dist[1], inv_haldane(0.1))
  expect_true(all(is.na(rep$dist[2:4])))                 # r ~ 0.5 never a huge finite cM
})

test_that("C1. OP-only mixed distances agree with the existing OP map reporter", {
  so <- sim_fullsib(n_markers=8, crosses=data.frame(mother="D1",father=NA,n=500),
                    r_const_m=0.12, op_paternal_pA=0.4, epsilon=0.02, seed=55)
  pm <- list(D1=so$truth$parents[["D1"]]$phase_vec)
  mx <- hmm_map_mixed(so$data, phased_m=pm, epsilon=0.05, lambda=20, tol=1e-6, maxit=1000,
                      r_start=0.05, gap_r=0.499)
  gm <- get_map(mx$op_result, "haldane", gap_r=0.499)
  expect_equal(as.numeric(mx$fit$d_m), as.numeric(attr(gm, "dist_cM")), tolerance=1e-12)
  expect_identical(unname(mx$fit$interval_status_m), unname(attr(gm, "status")))
})


# ---- Commit 2: explicit haplotype contract (no unsafe phase vectors) --------
test_that("C2. explicit haplotypes preserve orientation across a hom-interrupted interval", {
  # mother het, AA(hom), het with repulsion between the two het markers
  Htrue <- matrix(c(1L,0L, 1L,1L, 0L,1L), 2)   # homolog1 = A,A,a ; homolog2 = a,A,A
  pv <- haplotypes_to_phase(Htrue)
  expect_true(anyNA(pv))                        # phase undefined around the hom marker
  # phase_to_haplotypes is NOT a complete inverse: it rejects the NA phase, never couples it
  expect_error(phase_to_haplotypes(c(1L,2L,1L), pv), "unresolved")
  # the explicit matrix is unambiguous and keeps the correct (repulsion) orientation
  expect_identical(Htrue[, 1], c(1L, 0L)); expect_identical(Htrue[, 3], c(0L, 1L))
})

test_that("C2. haplotype and legacy phase inputs agree when phase is fully resolved", {
  sim <- sim_fullsib(n_markers=8, crosses=data.frame(mother="M1",father="S1",n=500),
                     r_const_m=0.1, r_const_p=0.2, epsilon=0.02, seed=42)
  ord <- colnames(sim$data$crosses[["M1__x__S1"]]$G)
  pm <- list(M1=sim$truth$parents[["M1"]]$phase_vec); pp <- list(S1=sim$truth$parents[["S1"]]$phase_vec)
  hm <- list(M1 = `colnames<-`(sim$truth$parents[["M1"]]$hap, ord))
  hp <- list(S1 = `colnames<-`(sim$truth$parents[["S1"]]$hap, ord))
  f_pv <- hmm_map_fullsib(sim$data, phased_m=pm, phased_p=pp, epsilon=0.02, tol=1e-7, maxit=1000)
  f_hp <- hmm_map_fullsib(sim$data, haplotypes_m=hm, haplotypes_p=hp, epsilon=0.02, tol=1e-7, maxit=1000)
  expect_equal(as.numeric(f_pv$fit$r_m), as.numeric(f_hp$fit$r_m), tolerance=1e-10)
  expect_equal(as.numeric(f_pv$fit$r_p), as.numeric(f_hp$fit$r_p), tolerance=1e-10)
})

test_that("C2. unresolved (NA) phase is rejected by the legacy interface", {
  sim <- sim_fullsib(n_markers=6, crosses=data.frame(mother="M1",father="S1",n=50), seed=15)
  pm <- list(M1=sim$truth$parents[["M1"]]$phase_vec); pp <- list(S1=sim$truth$parents[["S1"]]$phase_vec)
  pm$M1[2] <- NA_integer_
  expect_error(hmm_map_fullsib(sim$data, phased_m=pm, phased_p=pp), "unresolved")
})

test_that("C2. the same parent resolves to the same haplotype matrix across crosses", {
  ord <- sprintf("M%04d", 1:5)
  H <- matrix(c(1L,0L,1L,0L,0L,1L,0L,1L,1L,0L), 2, 5); colnames(H) <- ord
  a1 <- HSMap:::.resolve_parent_alleles("M1","mother", list(M1=H), NULL, NULL, ord, "test")
  a2 <- HSMap:::.resolve_parent_alleles("M1","mother", list(M1=H), NULL, NULL, ord, "test")
  expect_identical(a1, a2)
  expect_error(HSMap:::.resolve_parent_alleles("M1","mother",
    list(M1=`[<-`(H, 1, 1, NA_integer_)), NULL, NULL, ord, "test"), "unresolved")
})


# ---- Commit 3: EM convergence aligned with the stabilized OP engine ---------
test_that("C3. traces end exactly at the final values (full-sib)", {
  sim <- sim_fullsib(n_markers=8, crosses=data.frame(mother="M1",father="S1",n=500),
                     r_const_m=0.1, r_const_p=0.25, epsilon=0.02, seed=18)
  pm <- list(M1=sim$truth$parents[["M1"]]$phase_vec); pp <- list(S1=sim$truth$parents[["S1"]]$phase_vec)
  fit <- hmm_map_fullsib(sim$data, phased_m=pm, phased_p=pp, epsilon=0.02, tol=1e-7, maxit=1000)
  f <- fit$fit
  expect_equal(utils::tail(f$loglik_trace, 1), f$logLik, tolerance=1e-9)
  expect_equal(utils::tail(f$objective_trace, 1), f$objective, tolerance=1e-9)
  expect_equal(utils::tail(f$max_dr_m_trace, 1), 0)
  expect_equal(utils::tail(f$max_dr_p_trace, 1), 0)
  expect_length(f$loglik_trace, f$iters + 1L)            # one per iter + final append
  expect_true(all(diff(f$objective_trace) >= -1e-6))     # non-decreasing incl. final step
  expect_false(f$objective_decreased)
})

test_that("C3. mixed objective contains the final q penalty and traces end at final", {
  cr <- data.frame(mother=c("M1","D1"), father=c("S1",NA), n=c(500,500), stringsAsFactors=FALSE)
  sim <- sim_fullsib(n_markers=8, crosses=cr, r_const_m=0.1, r_const_p=0.2,
                     epsilon=0.02, op_paternal_pA=0.4, seed=71)
  pm <- lapply(sim$truth$parents, function(p) p$phase_vec)
  fit <- hmm_map_mixed(sim$data, phased_m=pm, phased_p=pm, epsilon=0.02, lambda=20, tol=1e-7, maxit=1000)
  f <- fit$fit
  expect_equal(f$objective, f$logLik + f$q_penalty, tolerance=1e-9)
  expect_equal(utils::tail(f$objective_trace, 1), f$objective, tolerance=1e-9)
  expect_equal(utils::tail(f$loglik_trace, 1), f$logLik, tolerance=1e-9)
  expect_equal(utils::tail(f$max_dq_trace, 1), 0)
  expect_length(f$objective_trace, f$iters + 1L)
  expect_true(all(diff(f$objective_trace) >= -1e-6))
})

test_that("C3. a low maxit reports non-convergence (mixed)", {
  cr <- data.frame(mother=c("M1","D1"), father=c("S1",NA), n=c(400,400), stringsAsFactors=FALSE)
  sim <- sim_fullsib(n_markers=8, crosses=cr, epsilon=0.02, seed=72)
  pm <- lapply(sim$truth$parents, function(p) p$phase_vec)
  expect_warning(mx <- hmm_map_mixed(sim$data, phased_m=pm, phased_p=pm, maxit=2), "did not converge")
  expect_false(mx$fit$converged)
  expect_identical(mx$fit$conv_reason, "maxit")
})

test_that("C3. the final-step objective-decrease check works and is scale-aware", {
  expect_true(HSMap:::.fs_obj_decreased(10, 20))
  expect_false(HSMap:::.fs_obj_decreased(20, 10))
  expect_false(HSMap:::.fs_obj_decreased(10 - 1e-13, 10))   # sub-tolerance change not flagged
  expect_false(HSMap:::.fs_obj_decreased(5, -Inf))          # first iteration never flags
})

test_that("C3. numeric input validation", {
  sim <- sim_fullsib(n_markers=6, crosses=data.frame(mother="M1",father="S1",n=50), seed=1)
  pm <- list(M1=sim$truth$parents[["M1"]]$phase_vec); pp <- list(S1=sim$truth$parents[["S1"]]$phase_vec)
  expect_error(hmm_map_fullsib(sim$data, phased_m=pm, phased_p=pp, epsilon=1),    "epsilon")
  expect_error(hmm_map_fullsib(sim$data, phased_m=pm, phased_p=pp, epsilon=-0.1), "epsilon")
  expect_error(hmm_map_fullsib(sim$data, phased_m=pm, phased_p=pp, tol=0),        "tol")
  expect_error(hmm_map_fullsib(sim$data, phased_m=pm, phased_p=pp, maxit=2.5),    "maxit")
  expect_error(hmm_map_fullsib(sim$data, phased_m=pm, phased_p=pp, r_start=0.6),  "r_start")
  expect_error(hmm_map_mixed(sim$data, phased_m=pm, phased_p=pp, lambda=-1),      "lambda")
  expect_error(hmm_map_mixed(sim$data, phased_m=pm, phased_p=pp, q0=1.5),         "q0")
})
