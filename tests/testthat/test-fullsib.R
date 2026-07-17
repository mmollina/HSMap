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


# ---- Commit 4: sex-specific consensus map scope + pooling metadata ----------
test_that("C4. output declares sex-specific consensus scope and per-parent pooling", {
  cr <- data.frame(mother=c("M1","M2"), father=c("S1","S2"), n=c(400,400), stringsAsFactors=FALSE)
  sim <- sim_fullsib(n_markers=6, crosses=cr, r_const_m=0.1, r_const_p=0.2, epsilon=0.01, seed=301)
  pm <- lapply(sim$truth$parents, function(p) p$phase_vec)
  fit <- hmm_map_fullsib(sim$data, phased_m=pm, phased_p=pm, epsilon=0.01, tol=1e-6, maxit=1000)
  expect_identical(fit$map_scope, "sex_specific_consensus")
  expect_identical(fit$fit$map_scope, "sex_specific_consensus")
  # two mothers pool into ONE maternal map; two sires into ONE paternal map
  expect_setequal(colnames(fit$fit$maternal_meioses_by_mother), c("M1","M2"))
  expect_setequal(colnames(fit$fit$paternal_meioses_by_sire), c("S1","S2"))
  # each mother contributes ~n meioses per interval to the shared maternal map
  expect_true(all(abs(fit$fit$maternal_meioses_by_mother[, "M1"] - 400) < 1e-6))
  # row sums equal the total maternal meioses per interval
  expect_equal(rowSums(fit$fit$maternal_meioses_by_mother), as.numeric(fit$fit$meiosis_count),
               tolerance=1e-8, ignore_attr=TRUE)
})

test_that("C4. a repeated mother/sire pools into one shared map (no parent-specific map)", {
  # repeated mother M1 (two sires) and repeated sire S1 (two mothers)
  cr <- data.frame(mother=c("M1","M1","M2"), father=c("S1","S2","S1"),
                   n=c(300,300,300), stringsAsFactors=FALSE)
  sim <- sim_fullsib(n_markers=6, crosses=cr, r_const_m=0.1, r_const_p=0.2, epsilon=0.01, seed=302)
  pm <- lapply(sim$truth$parents, function(p) p$phase_vec)
  fit <- hmm_map_fullsib(sim$data, phased_m=pm, phased_p=pm, epsilon=0.01, tol=1e-6, maxit=1000)
  # ONE maternal map and ONE paternal map (single r_m / r_p vectors), not per parent
  expect_length(fit$fit$r_m, 5L); expect_length(fit$fit$r_p, 5L)
  # repeated mother M1 pools BOTH her crosses (600 meioses/interval); S1 pools its two
  expect_true(all(abs(fit$fit$maternal_meioses_by_mother[, "M1"] - 600) < 1e-6))
  expect_true(all(abs(fit$fit$paternal_meioses_by_sire[, "S1"] - 600) < 1e-6))
  expect_setequal(fit$maternal_crosses, c("M1__x__S1","M1__x__S2","M2__x__S1"))
})

test_that("C4. mixed pooling: maternal counts include OP, paternal only full-sib sires", {
  cr <- data.frame(mother=c("M1","D1"), father=c("S1",NA), n=c(400,400), stringsAsFactors=FALSE)
  sim <- sim_fullsib(n_markers=6, crosses=cr, r_const_m=0.1, r_const_p=0.2,
                     epsilon=0.01, op_paternal_pA=0.4, seed=303)
  pm <- lapply(sim$truth$parents, function(p) p$phase_vec)
  fit <- hmm_map_mixed(sim$data, phased_m=pm, phased_p=pm, epsilon=0.01, lambda=20, tol=1e-6, maxit=1000)
  expect_identical(fit$map_scope, "sex_specific_consensus")
  expect_setequal(colnames(fit$fit$maternal_meioses_by_mother), c("M1","D1"))  # OP dam included
  expect_setequal(colnames(fit$fit$paternal_meioses_by_sire), "S1")            # only full-sib sire
  expect_setequal(fit$maternal_crosses, c("M1__x__S1","D1"))
  expect_setequal(fit$paternal_crosses, "M1__x__S1")
})


# ---- Commit 5: identifiability and information diagnostics ------------------
# Build a one-cross HSMap.data directly from explicit parental haplotypes.
.mk_fs_data <- function(Hm, Hp, rm, rp, n, seed) {
  set.seed(seed); z <- ncol(Hm); ord <- sprintf("M%04d", 1:z)
  meio <- function(r) { P <- matrix(0L, n, z); P[, 1] <- sample(1:2, n, TRUE)
    for (k in 2:z) { sw <- rbinom(n, 1, r[k-1]) == 1; P[, k] <- ifelse(sw, 3L - P[, k-1], P[, k-1]) }; P }
  Pm <- meio(rm); Pp <- meio(rp)
  ma <- sapply(1:z, function(k) Hm[cbind(Pm[, k], k)]); pa <- sapply(1:z, function(k) Hp[cbind(Pp[, k], k)])
  G <- ma + pa; storage.mode(G) <- "integer"; rownames(G) <- sprintf("o%03d", 1:n); colnames(G) <- ord
  gm <- as.integer(Hm[1, ] + Hm[2, ]); gf <- as.integer(Hp[1, ] + Hp[2, ]); names(gm) <- names(gf) <- ord
  cr <- list("M1__x__S1" = list(cross_id="M1__x__S1", mother_id="M1", father_id="S1",
             family_type="known_sire_genotyped", offspring=rownames(G), M=gm, F=gf, G=G))
  dat <- list(G_list=list("M1__x__S1"=G), M_list=list("M1__x__S1"=gm),
    alleles=data.frame(marker_id=ord, REF="A", ALT="B", chrom=1L, position=1:z, stringsAsFactors=FALSE),
    pedigree=NULL,
    cross_table=data.frame(cross_id="M1__x__S1", mother_id="M1", father_id="S1",
      family_type="known_sire_genotyped", n_offspring=n, mother_genotyped=TRUE,
      father_genotyped=TRUE, stringsAsFactors=FALSE),
    crosses=cr, parent_genotypes=list(M1=gm, S1=gf), F_list=list("M1__x__S1"=gf))
  class(dat) <- "HSMap.data"
  list(data=dat, Hm=`colnames<-`(Hm, ord), Hp=`colnames<-`(Hp, ord), ord=ord)
}

test_that("C5. structural informativeness and exchangeability are detected", {
  Am <- matrix(c(1L,0L, 0L,1L, 1L,0L, 0L,1L), 2, 4)   # mother het at all markers
  built_matmaonly <- list(list(Am=Am, Ap=matrix(1L,2,4), G=matrix(0L,10,4)))  # sire AA -> uninformative
  i1 <- HSMap:::.fs_informativeness(built_matmaonly, built_matmaonly, 3L)
  expect_true(all(i1$mat_inf)); expect_true(all(!i1$pat_inf)); expect_false(i1$globally_exchangeable)
  built_exch <- list(list(Am=Am, Ap=Am, G=matrix(0L,10,4)))     # identical parents
  i2 <- HSMap:::.fs_informativeness(built_exch, built_exch, 3L)
  expect_true(i2$globally_exchangeable)
})

test_that("C5. only-maternal-informative: paternal intervals flagged, not forced", {
  Hm <- matrix(c(1L,0L, 1L,0L, 0L,1L, 1L,0L, 0L,1L, 1L,0L), 2, 6)  # mother het everywhere
  Hp <- matrix(1L, 2, 6)                                            # sire AA everywhere
  d <- .mk_fs_data(Hm, Hp, rep(0.1,5), rep(0.1,5), 1500, 11)
  fit <- hmm_map_fullsib(d$data, haplotypes_m=list(M1=d$Hm), haplotypes_p=list(S1=d$Hp),
                         epsilon=0.01, tol=1e-6, maxit=1000)
  expect_true(all(fit$fit$interval_status_p == "insufficient_information"))
  expect_true(all(is.na(fit$fit$d_p)))                             # not forced into a distance
  expect_false(any(fit$fit$identifiability$paternal_informative))
  expect_true(all(fit$fit$identifiability$maternal_informative))
})

test_that("C5. exchangeable parents: non-identifiable labels and symmetric likelihood", {
  H <- matrix(c(1L,0L, 0L,1L, 1L,0L, 0L,1L, 1L,0L), 2, 5)          # identical mother & sire
  d <- .mk_fs_data(H, H, rep(0.1,4), rep(0.3,4), 1200, 22)
  fit <- hmm_map_fullsib(d$data, haplotypes_m=list(M1=d$Hm), haplotypes_p=list(S1=d$Hp),
                         epsilon=0.01, tol=1e-6, maxit=1000)
  expect_false(fit$fit$identifiable_labels)
  expect_true(all(fit$fit$interval_status_m == "nonidentifiable_exchangeable"))
  # likelihood is symmetric under r_m <-> r_p (labels not identifiable)
  cr <- d$data$crosses[["M1__x__S1"]]
  rmv <- c(0.05,0.2,0.1,0.4); rpv <- c(0.3,0.15,0.25,0.05)
  ll1 <- fs_loglik_cpp(cr$G, matrix(as.integer(d$Hm),2), matrix(as.integer(d$Hp),2), rmv, rpv, 0.01)
  ll2 <- fs_loglik_cpp(cr$G, matrix(as.integer(d$Hm),2), matrix(as.integer(d$Hp),2), rpv, rmv, 0.01)
  expect_equal(ll1, ll2, tolerance=1e-9)
})

test_that("C5. distinguishable maps: unique optimum on the (r_m, r_p) grid + start-insensitive", {
  Hm <- matrix(c(1L,0L, 0L,1L, 1L,0L),2,3); Hp <- matrix(c(1L,0L, 1L,0L, 0L,1L),2,3)  # distinct phase
  d <- .mk_fs_data(Hm, Hp, c(0.08,0.35), c(0.40,0.08), 3000, 33)
  cr <- d$data$crosses[["M1__x__S1"]]; AM <- matrix(as.integer(d$Hm),2); AP <- matrix(as.integer(d$Hp),2)
  grid <- seq(0.01, 0.49, by=0.02)
  # first interval: profile ll over (rm, rp) with the 2nd interval at a fixed near-truth value
  best <- c(NA,NA); bestll <- -Inf
  for (a in grid) for (b in grid) {
    ll <- fs_loglik_cpp(cr$G, AM, AP, c(a,0.35), c(b,0.08), 0.01)
    if (ll > bestll) { bestll <- ll; best <- c(a,b) }
  }
  expect_lt(abs(best[1] - 0.08), 0.05)      # maternal r_1 optimum near truth 0.08
  expect_lt(abs(best[2] - 0.40), 0.06)      # paternal r_1 optimum near truth 0.40
  expect_gt(abs(best[1] - best[2]), 0.15)   # clearly distinct -> not exchangeable
  # start-insensitive convergence
  f1 <- hmm_map_fullsib(d$data, haplotypes_m=list(M1=d$Hm), haplotypes_p=list(S1=d$Hp), r_start=0.05, tol=1e-7)
  f2 <- hmm_map_fullsib(d$data, haplotypes_m=list(M1=d$Hm), haplotypes_p=list(S1=d$Hp), r_start=0.45, tol=1e-7)
  expect_equal(as.numeric(f1$fit$r_m), as.numeric(f2$fit$r_m), tolerance=1e-4)
  expect_equal(as.numeric(f1$fit$r_p), as.numeric(f2$fit$r_p), tolerance=1e-4)
  expect_true(f1$fit$identifiable_labels)
})


# ---- Commit 6: cross-aware OP dispatch safety ------------------------------
# Two OP-modelled (known-but-untyped-sire) crosses that SHARE one mother; cross_id !=
# mother_id, so the legacy G_list (keyed by mother) cannot represent them 1:1.
.mk_untyped_shared_mother <- function(z = 6, n = 400, seed = 1) {
  ord <- sprintf("M%04d", 1:z)
  Hm <- matrix(c(1L, 0L), 2, z)                      # mother het everywhere (coupling)
  gm <- as.integer(Hm[1, ] + Hm[2, ]); names(gm) <- ord
  rm <- rep(0.1, z - 1)
  mk_cross <- function(sire, sd) {
    set.seed(sd)
    P <- matrix(0L, n, z); P[, 1] <- sample(1:2, n, TRUE)
    for (k in 2:z) { sw <- rbinom(n, 1, rm[k-1]) == 1; P[, k] <- ifelse(sw, 3L - P[, k-1], P[, k-1]) }
    ma <- sapply(1:z, function(k) Hm[cbind(P[, k], k)])
    pa <- matrix(rbinom(n * z, 1, 0.4), n, z)        # independent paternal (OP)
    G <- ma + pa; storage.mode(G) <- "integer"; rownames(G) <- sprintf("%s_o%03d", sire, 1:n); colnames(G) <- ord
    list(cross_id = paste0("M1__x__", sire), mother_id = "M1", father_id = sire,
         family_type = "known_sire_untyped", offspring = rownames(G),
         M = gm, F = stats::setNames(rep(NA_integer_, z), ord), G = G)
  }
  c1 <- mk_cross("S1", seed + 1); c2 <- mk_cross("S2", seed + 2)
  cids <- c(c1$cross_id, c2$cross_id)
  dat <- list(G_list = stats::setNames(list(c1$G, c2$G), cids),
    M_list = stats::setNames(list(gm, gm), cids),
    alleles = data.frame(marker_id = ord, REF = "A", ALT = "B", chrom = 1L, position = 1:z, stringsAsFactors = FALSE),
    pedigree = NULL,
    cross_table = data.frame(cross_id = cids, mother_id = "M1", father_id = c("S1", "S2"),
      family_type = "known_sire_untyped", n_offspring = n, mother_genotyped = TRUE,
      father_genotyped = FALSE, stringsAsFactors = FALSE),
    crosses = stats::setNames(list(c1, c2), cids), parent_genotypes = list(M1 = gm),
    F_list = stats::setNames(list(NULL, NULL), cids))
  class(dat) <- "HSMap.data"
  list(data = dat, Hm = `colnames<-`(Hm, ord), ord = ord, rm = rm, n = n)
}

test_that("C6. an ordinary one-mother OP family dispatches to the legacy engine", {
  so <- sim_fullsib(n_markers=8, crosses=data.frame(mother="D1",father=NA,n=400),
                    r_const_m=0.12, op_paternal_pA=0.4, epsilon=0.02, seed=61)
  pm <- list(D1 = so$truth$parents[["D1"]]$phase_vec)
  mx <- hmm_map_mixed(so$data, phased_m=pm, epsilon=0.05, lambda=20, tol=1e-6, maxit=1000, r_start=0.05)
  expect_true(mx$dispatched)                          # cross_id == mother, unique -> safe
})

test_that("C6. untyped-sire crosses sharing a mother do NOT silently dispatch", {
  d <- .mk_untyped_shared_mother(z=6, n=400, seed=5)
  hm <- list(M1 = d$Hm)
  # default: an untyped sire errors
  expect_error(hmm_map_mixed(d$data, haplotypes_m=hm), "untyped")
  # OP fallback: cross-aware path (NOT legacy dispatch), one maternal map pooling BOTH
  mx <- hmm_map_mixed(d$data, haplotypes_m=hm, untyped_sire="open_pollinated",
                      epsilon=0.02, lambda=20, tol=1e-6, maxit=1000)
  expect_false(mx$dispatched)                          # cross identity preserved, not dispatched
  expect_setequal(mx$contributing_crosses, c("M1__x__S1","M1__x__S2"))
  expect_true(all(abs(mx$fit$maternal_meioses_by_mother[, "M1"] - 2 * d$n) < 1e-6))  # both pooled
  expect_lt(sqrt(mean((as.numeric(mx$fit$r_m) - d$rm)^2)), 0.03)  # maternal map recovered
})

test_that("C6. a single untyped-sire OP-fallback cross (cross_id != mother) is cross-aware", {
  d <- .mk_untyped_shared_mother(z=6, n=600, seed=9)
  # keep only one cross so mothers are unique but cross_id != mother_id
  d$data$crosses <- d$data$crosses["M1__x__S1"]
  d$data$cross_table <- d$data$cross_table[1, , drop=FALSE]
  d$data$G_list <- d$data$G_list["M1__x__S1"]; d$data$M_list <- d$data$M_list["M1__x__S1"]
  d$data$F_list <- d$data$F_list["M1__x__S1"]
  mx <- hmm_map_mixed(d$data, haplotypes_m=list(M1=d$Hm), untyped_sire="open_pollinated",
                      epsilon=0.02, lambda=20, tol=1e-6, maxit=1000)
  expect_false(mx$dispatched)                          # cross_id != mother -> not dispatched
  expect_lt(sqrt(mean((as.numeric(mx$fit$r_m) - d$rm)^2)), 0.035)
})
