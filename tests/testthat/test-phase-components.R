## Regression guard for .pf_components().
##
## The original implementation labelled vertices when POPPED and rebuilt the
## stack with `stack[-length(stack)]` on every pop. On a dense phase graph --
## the normal case for a linkage group when `min_phase_lod` admits most pairs --
## a vertex could be pushed once per visited neighbour, so the stack grew to
## O(n^2) and each pop copied it. Cost grew like n^4: a 400-marker clique took
## ~26 s inside phase_from_pairwise(), a 600-marker one ~65 s.
##
## These tests pin BOTH properties that matter: the labelling is unchanged, and
## a dense graph is handled in roughly quadratic time.

test_that(".pf_components labels components correctly", {
  # single vertex
  expect_equal(HSMap:::.pf_components(matrix(FALSE, 1, 1)), 1L)

  # all singletons -> one component each, in ascending vertex order
  expect_equal(HSMap:::.pf_components(matrix(FALSE, 5, 5)), 1:5)

  # clique -> one component
  S <- matrix(TRUE, 12, 12); diag(S) <- FALSE
  expect_equal(HSMap:::.pf_components(S), rep(1L, 12))

  # two disconnected blocks -> ids assigned in ascending vertex order
  S <- matrix(FALSE, 6, 6)
  S[1:3, 1:3] <- TRUE; S[4:6, 4:6] <- TRUE; diag(S) <- FALSE
  expect_equal(HSMap:::.pf_components(S), c(1L, 1L, 1L, 2L, 2L, 2L))

  # path graph: connected, and deep recursion is not used
  n <- 200L; S <- matrix(FALSE, n, n); i <- 1:(n - 1L)
  S[cbind(i, i + 1L)] <- TRUE; S[cbind(i + 1L, i)] <- TRUE
  expect_equal(HSMap:::.pf_components(S), rep(1L, n))
})

test_that(".pf_components agrees with an independent reference on random graphs", {
  # Reference: repeated boolean reachability closure -- a different algorithm,
  # so agreement is genuine corroboration rather than the same code twice.
  ref <- function(Sup) {
    n <- nrow(Sup); comp <- integer(n); cid <- 0L
    for (s in seq_len(n)) {
      if (comp[s] != 0L) next
      cid <- cid + 1L
      reach <- rep(FALSE, n); reach[s] <- TRUE
      repeat {
        nxt <- reach | apply(Sup[reach, , drop = FALSE], 2L, any)
        if (identical(nxt, reach)) break
        reach <- nxt
      }
      comp[reach & comp == 0L] <- cid
    }
    comp
  }
  set.seed(11)
  for (n in c(6L, 25L, 80L)) {
    for (p in c(0.01, 0.08, 0.5, 0.95)) {
      M <- matrix(runif(n * n) < p, n, n); M <- M | t(M); diag(M) <- FALSE
      expect_equal(HSMap:::.pf_components(M), ref(M),
                   info = sprintf("n=%d p=%.2f", n, p))
    }
  }
})

test_that(".pf_components stays near-quadratic on a dense graph", {
  # A 400-vertex clique. Quadratic: milliseconds. The old quartic behaviour took
  # tens of seconds, so this bound is generous while still catching a regression.
  S <- matrix(TRUE, 400L, 400L); diag(S) <- FALSE
  el <- system.time(cc <- HSMap:::.pf_components(S))[["elapsed"]]
  expect_equal(cc, rep(1L, 400L))
  expect_lt(el, 5)
})
