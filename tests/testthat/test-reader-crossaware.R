# Cross-aware reader validation (m6-known-sire-core, Correct 7/7).

# Write pedigree + genotype data frames to temp CSVs and read them back.
.write_read <- function(ped, geno) {
  pf <- tempfile(fileext = ".csv"); gf <- tempfile(fileext = ".csv")
  utils::write.csv(ped, pf, row.names = FALSE)
  utils::write.csv(geno, gf, row.names = FALSE)
  read_HSMap_data(pf, gf)
}

# minimal genotype frame: markers x (meta + samples). `samples` is a named list of
# integer genotype vectors (length = n markers).
.geno_df <- function(markers, samples) {
  d <- data.frame(marker_id = markers, REF = "A", ALT = "B",
                  chrom = 1L, position = seq_along(markers), stringsAsFactors = FALSE)
  for (nm in names(samples)) d[[nm]] <- samples[[nm]]
  d
}

MK <- sprintf("M%02d", 1:4)

test_that("a known sire that is a pedigree founder is read as genotyped", {
  ped <- data.frame(
    id = c("M1","S1","o1","o2","o3"),
    mother = c(NA,NA,"M1","M1","M1"),
    father = c(NA,NA,"S1","S1","S1"),
    generation = c(1,1,2,2,2),
    family_id = c("M1","S1","c1","c1","c1"), stringsAsFactors = FALSE)
  geno <- .geno_df(MK, list(M1=c(1,1,1,1), S1=c(1,0,2,1),
                            o1=c(1,0,1,1), o2=c(2,1,2,1), o3=c(1,0,1,2)))
  d <- .write_read(ped, geno)
  expect_identical(d$cross_table$family_type, "known_sire_genotyped")
  expect_identical(d$cross_table$cross_id, "M1__x__S1")
  expect_false(is.na(d$parent_genotypes[["S1"]][1]))          # sire genotype from the file
})

test_that("a known sire absent from the pedigree but present in the genotype file is genotyped", {
  ped <- data.frame(                                          # NO S1 row
    id = c("M1","o1","o2"), mother = c(NA,"M1","M1"),
    father = c(NA,"S1","S1"), generation = c(1,2,2),
    family_id = c("M1","c1","c1"), stringsAsFactors = FALSE)
  geno <- .geno_df(MK, list(M1=c(1,1,1,1), S1=c(1,1,0,2), o1=c(1,1,1,2), o2=c(2,1,0,2)))
  d <- .write_read(ped, geno)
  expect_identical(d$cross_table$family_type, "known_sire_genotyped")  # parent row not required
  expect_true("S1" %in% names(d$parent_genotypes))
})

test_that("a known sire listed but not genotyped is known_sire_untyped (with a warning)", {
  ped <- data.frame(
    id = c("M1","o1","o2"), mother = c(NA,"M1","M1"),
    father = c(NA,"S1","S1"), generation = c(1,2,2),
    family_id = c("M1","c1","c1"), stringsAsFactors = FALSE)
  geno <- .geno_df(MK, list(M1=c(1,1,1,1), o1=c(1,1,1,2), o2=c(2,1,0,2)))  # no S1 column
  expect_warning(d <- .write_read(ped, geno), "untyped")
  expect_identical(d$cross_table$family_type, "known_sire_untyped")
  expect_true(is.na(d$parent_genotypes[["S1"]][1]))
})

test_that("conflicting parent records for one offspring are an error", {
  ped <- data.frame(
    id = c("M1","M2","o1","o1"), mother = c(NA,NA,"M1","M2"),
    father = c(NA,NA,"S1","S1"), generation = c(1,1,2,2),
    family_id = c("M1","M2","c1","c2"), stringsAsFactors = FALSE)
  geno <- .geno_df(MK, list(M1=c(1,1,1,1), M2=c(1,1,1,1), o1=c(1,1,1,2)))
  expect_error(.write_read(ped, geno), "conflicting parents")
})

test_that("a duplicate offspring record with identical parents is de-duplicated", {
  ped <- data.frame(
    id = c("M1","o1","o1","o2"), mother = c(NA,"M1","M1","M1"),
    father = c(NA,NA,NA,NA), generation = c(1,2,2,2),
    family_id = c("M1","M1","M1","M1"), stringsAsFactors = FALSE)
  geno <- .geno_df(MK, list(M1=c(1,1,1,1), o1=c(1,1,1,2), o2=c(2,1,0,2)))
  d <- .write_read(ped, geno)
  expect_identical(d$cross_table$n_offspring, 2L)             # o1 counted once
  expect_identical(nrow(d$crosses[["M1"]]$G), 2L)
})

test_that("repeated parents across crosses are stored once and shared", {
  ped <- data.frame(
    id = c("M1","o1","o2","o3","o4"),
    mother = c(NA,"M1","M1","M1","M1"),
    father = c(NA,"S1","S1","S2","S2"),
    generation = c(1,2,2,2,2),
    family_id = c("M1","c1","c1","c2","c2"), stringsAsFactors = FALSE)
  geno <- .geno_df(MK, list(M1=c(1,1,1,1), S1=c(1,0,1,1), S2=c(0,1,1,0),
                            o1=c(1,0,1,1), o2=c(2,1,2,1), o3=c(0,1,1,0), o4=c(1,1,1,1)))
  d <- .write_read(ped, geno)
  expect_setequal(d$cross_table$cross_id, c("M1__x__S1","M1__x__S2"))
  # M1 stored once; both crosses see the same maternal genotype object
  expect_identical(d$crosses[["M1__x__S1"]]$M, d$crosses[["M1__x__S2"]]$M)
  expect_identical(d$crosses[["M1__x__S1"]]$M, stats::setNames(d$parent_genotypes[["M1"]], MK))
})

test_that("open-pollinated reading is unchanged (legacy G_list keyed by mother)", {
  ped <- data.frame(
    id = c("M1","o1","o2"), mother = c(NA,"M1","M1"),
    father = c(NA,NA,NA), generation = c(1,2,2),
    family_id = c("M1","M1","M1"), stringsAsFactors = FALSE)
  geno <- .geno_df(MK, list(M1=c(1,1,1,1), o1=c(1,1,1,2), o2=c(2,1,0,2)))
  d <- .write_read(ped, geno)
  expect_identical(names(d$G_list), "M1")                     # legacy name preserved
  expect_identical(d$cross_table$family_type, "open_pollinated")
  expect_identical(d$cross_table$father_id, "__unknown_sire__")
})
