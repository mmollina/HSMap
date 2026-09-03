## plot_map_lmv(): LinkageMapView bridge (Suggests-only dependency)

.mk_map <- function(r, markers) {
  structure(list(fit = list(r = r), order = markers), class = "HSMap.map")
}

test_that("plot_map_lmv renders a PDF and returns the tidy frame", {
  skip_if_not_installed("LinkageMapView")
  m1 <- .mk_map(c(0.05, 0.10, 0.02), paste0("a", 1:4))
  m2 <- .mk_map(c(0.08, 0.03),       paste0("b", 1:3))
  out <- tempfile(fileext = ".pdf")
  on.exit(unlink(out), add = TRUE)

  df <- NULL
  invisible(capture.output(                      # lmv prints its pdf sizing
    df <- plot_map_lmv(list(chrA = m1, chrB = m2), outfile = out)
  ))
  expect_true(file.exists(out) && file.size(out) > 0)
  expect_identical(names(df), c("group", "position", "locus"))
  expect_identical(unique(df$group), c("chrA", "chrB"))
  expect_identical(df$locus, c(paste0("a", 1:4), paste0("b", 1:3)))
  ## positions are the gap-aware cumulative cM (Haldane by default)
  expect_equal(df$position[df$group == "chrA"],
               unname(cumsum(c(0, inv_haldane(c(0.05, 0.10, 0.02))))),
               tolerance = 1e-10)
  ## unnamed input gets LG1... labels; single map accepted directly
  out2 <- tempfile(fileext = ".pdf"); on.exit(unlink(out2), add = TRUE)
  invisible(capture.output(df2 <- plot_map_lmv(m1, outfile = out2)))
  expect_identical(unique(df2$group), "LG1")
})

test_that("maps with gaps split into per-segment bars (positions reset)", {
  skip_if_not_installed("LinkageMapView")
  ## interval 2 is a no-linkage gap (r = 0.5): segments a1-a2 and a3-a5
  m <- .mk_map(c(0.05, 0.5, 0.10, 0.02), paste0("a", 1:5))
  out <- tempfile(fileext = ".pdf")
  on.exit(unlink(out), add = TRUE)
  invisible(capture.output(
    df <- plot_map_lmv(list(chr7 = m), outfile = out)
  ))
  expect_true(file.exists(out) && file.size(out) > 0)
  expect_identical(unique(df$group), c("chr7.1", "chr7.2"))
  expect_equal(df$position[df$locus == "a3"], 0)          # reset after the gap
  expect_equal(df$position[df$locus == "a5"],
               unname(sum(inv_haldane(c(0.10, 0.02)))), tolerance = 1e-10)

  ## split_gaps = FALSE refuses, mirroring plot_map_list()
  expect_error(plot_map_lmv(list(chr7 = m), outfile = out, split_gaps = FALSE),
               "gap")
})

test_that("input validation", {
  skip_if_not_installed("LinkageMapView")
  m <- .mk_map(0.1, c("x1", "x2"))
  expect_error(plot_map_lmv(m), "outfile")
  expect_error(plot_map_lmv(list(m, 1), outfile = tempfile()), "HSMap.map")
  expect_error(plot_map_lmv(list(), outfile = tempfile()), "non-empty")
})
