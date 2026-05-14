# test_getDataForKusterCompounds.R

# Unitary tests for getDataForKusterCompounds()

# to run just this file:
#   testthat::test_file("tests/R/test_getDataForKusterCompounds.R")

source("helpers.R")

# ==== TEST WORKBOOK ====

wb_simple <- make_workbook(
  fold_matrix = rbind(
    c("sh.index.sites", "FDR",  "MCF7.CompA.fold",    "MCF7.CompB.fold"),
    c("SITE1;",         "0.01", "-2.5",               "1.3"),
    c("SITE2;",         "0.50",  "0.8",               "-0.3")
  ),
  pvalue_matrix = rbind(
    c("sh.index.sites", NA,     "MCF7.CompA.p.value", "MCF7.CompB.p.value"),
    c("SITE1;",         "0.01", "0.010",              "0.300"),
    c("SITE2;",         "0.50", "0.800",              "0.900")
  )
)

result <- getDataForKusterCompounds(wb_simple, c("CompA", "CompB"), numCompounds = 2)
fValues         <- result$fValues
pValues         <- result$pValues
fdrValues       <- result$fdrValues
compoundsCL     <- result$compoundsCellLine
sitesCellLine   <- result$sitesCellLine


# ==== PSITES ====

test_that("extracts site list correctly", {
  expect_equal(sitesCellLine, c("SITE1;", "SITE2;"))
})

test_that("correct site number", {
  expect_equal(length(sitesCellLine), 2)
})


# ==== FOLD VALUES ====

test_that("obtain fold values correctly", {
  expect_equal(fValues[["SITE1;"]][["CompA"]], -2.5,  tolerance = 1e-9)
  expect_equal(fValues[["SITE1;"]][["CompB"]],  1.3,  tolerance = 1e-9)
  expect_equal(fValues[["SITE2;"]][["CompA"]],  0.8,  tolerance = 1e-9)
  expect_equal(fValues[["SITE2;"]][["CompB"]], -0.3,  tolerance = 1e-9)
})

test_that("clues of fValues match sitesCellLine", {
  expect_setequal(names(fValues), sitesCellLine)
})

test_that("values of fValues are numerical", {
  for (site in names(fValues)) {
    for (compound in names(fValues[[site]])) {
      expect_true(
        is.numeric(fValues[[site]][[compound]]),
        label = paste("fValues[[", site, "]][[", compound, "]] no es numérico")
      )
    }
  }
})


# ==== P-VALUES ====

test_that("obtain p-values correctly", {
  expect_equal(pValues[["SITE1;"]][["CompA"]], 0.010, tolerance = 1e-9)
  expect_equal(pValues[["SITE1;"]][["CompB"]], 0.300, tolerance = 1e-9)
  expect_equal(pValues[["SITE2;"]][["CompA"]], 0.800, tolerance = 1e-9)
})

test_that("clues of pValues match fValues", {
  expect_setequal(names(pValues), names(fValues))
  for (site in names(fValues)) {
    expect_setequal(names(pValues[[site]]), names(fValues[[site]]))
  }
})


# ==== FDR ====

test_that("obtain FDR correctly", {
  expect_equal(fdrValues[["SITE1;"]], 0.01, tolerance = 1e-9)
  expect_equal(fdrValues[["SITE2;"]], 0.50, tolerance = 1e-9)
})

test_that("FDR NA becomes 0", {
  # Note: sheet pvalue uses "" in col 2 (no NA) so openxlsx
  # does not discard the column when reading back with readWorkbook.
  wb_na_fdr <- make_workbook(
    fold_matrix = rbind(
      c("sh.index.sites", "FDR", "MCF7.CompA.fold"),
      c("SITE1;",          NA,    "-2.5")
    ),
    pvalue_matrix = rbind(
      c("sh.index.sites", "",    "MCF7.CompA.p.value"),
      c("SITE1;",          "",    "0.01")
    )
  )
  res <- getDataForKusterCompounds(wb_na_fdr, c("CompA"), numCompounds = 1)
  expect_equal(res$fdrValues[["SITE1;"]], 0.0, tolerance = 1e-9)
})


# ==== COMPOUNDSCELLLINE ====

test_that("clues of compoundsCellLine are short names", {
  expect_true("CompA" %in% names(compoundsCL))
  expect_true("CompB" %in% names(compoundsCL))
})

test_that("values of compoundsCellLine are full names of column", {
  expect_equal(compoundsCL[["CompA"]], "MCF7.CompA.fold")
  expect_equal(compoundsCL[["CompB"]], "MCF7.CompB.fold")
})


# ==== CREATED SHEETS ====

test_that("creates pvalue.select sheet", {
  expect_true("pvalue.select" %in% names(wb_simple))
})

test_that("creates fold.select sheet", {
  expect_true("fold.select" %in% names(wb_simple))
})

test_that("do not duplicates sheets if they already exist", {
  #second run
  getDataForKusterCompounds(wb_simple, c("CompA", "CompB"), numCompounds = 2)
  expect_equal(sum(names(wb_simple) == "pvalue.select"), 1)
  expect_equal(sum(names(wb_simple) == "fold.select"), 1)
})
