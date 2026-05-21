# test_makeKSlist.R
#
# Tests unitarios para makeKSlistOfKinaseDownstreamTargets()
#
# Para ejecutar:
#   testthat::test_file("tests/R/test_makeKSlist.R")

source("helpers.R")

# ==== TEST DATA ====

prob_kinases <- list(
  KinaseX = list("SITE1;" = 1.0, "SITE2;" = 0.0),
  KinaseY = list()
)

ratios <- list(
  KinaseX = list("SITE1;" = 1.0, "SITE2;" = 0.0),
  KinaseY = list()
)

fdr_values <- list(
  "SITE1;" = 0.010,   # <0.02 : FDR filter passes
  "SITE2;" = 0.500    # >0.02: not
)

p_values <- list(
  "SITE1;" = list(CompA = 0.01),
  "SITE2;" = list(CompA = 0.80)
)

#auxiliary function to create a clean workbook and run makeKSlist
run_makeks <- function(probs = prob_kinases, rats = ratios, fdr = fdr_values,
                       ratio_t = 0.5, prob_t = 0.5) {
  wb <- createWorkbook()
  addWorksheet(wb, "fold"); addWorksheet(wb, "pvalue")
  makeKSlistOfKinaseDownstreamTargets(wb, rats, probs, p_values, fdr,
                                       ratioThreshold       = ratio_t,
                                       probabilityThreshold = prob_t)
}


# ==== TESTS: TRIPLE FILTER ====

test_that("site included when overcomes prob, ratio y fdr", {
  result <- run_makeks()
  expect_true("SITE1" %in% result[["KinaseX"]])
})

test_that("site included at exact threshold boundary (>= VBA behaviour)", {
  # prob == probabilityThreshold and ratio == ratioThreshold must be INCLUDED
  probs_tmp <- list(KinaseX = list("SITE1;" = 0.5))   # exactly at threshold
  rats_tmp  <- list(KinaseX = list("SITE1;" = 0.5))   # exactly at threshold
  fdr_tmp   <- list("SITE1;" = 0.001)

  result <- run_makeks(probs = probs_tmp, rats = rats_tmp, fdr = fdr_tmp,
                       ratio_t = 0.5, prob_t = 0.5)
  expect_true("SITE1" %in% result[["KinaseX"]])
})

test_that("site excluded when FDR is high", {
  probs_tmp <- list(KinaseX = list("SITE2;" = 1.0))
  rats_tmp  <- list(KinaseX = list("SITE2;" = 1.0))
  fdr_tmp   <- list("SITE2;" = 0.5)  # >0.02: excluded

  result <- run_makeks(probs = probs_tmp, rats = rats_tmp, fdr = fdr_tmp)
  expect_false("SITE2" %in% result[["KinaseX"]])
})

test_that("excluded site when probability is low", {
  probs_tmp <- list(KinaseX = list("SITE1;" = 0.3))  # <0.5
  rats_tmp  <- list(KinaseX = list("SITE1;" = 1.0))
  fdr_tmp   <- list("SITE1;" = 0.001)

  result <- run_makeks(probs = probs_tmp, rats = rats_tmp, fdr = fdr_tmp)
  expect_false("SITE1" %in% result[["KinaseX"]])
})

test_that("excluded site when ratio is low", {
  probs_tmp <- list(KinaseX = list("SITE1;" = 1.0))
  rats_tmp  <- list(KinaseX = list("SITE1;" = 0.1))  # <0.5
  fdr_tmp   <- list("SITE1;" = 0.001)

  result <- run_makeks(probs = probs_tmp, rats = rats_tmp, fdr = fdr_tmp)
  expect_false("SITE1" %in% result[["KinaseX"]])
})


# ==== TESTS: RESIDUES EXCLUSION ====

test_that("exclude psites with residue (M", {
  sitio <- "PROT(M1);"
  probs_tmp <- list(KinaseX = setNames(list(1.0), sitio))
  rats_tmp  <- list(KinaseX = setNames(list(1.0), sitio))
  fdr_tmp   <- setNames(list(0.001), sitio)

  result <- run_makeks(probs = probs_tmp, rats = rats_tmp, fdr = fdr_tmp)
  expect_length(result[["KinaseX"]], 0)
})

test_that("exclude psites with residue (R", {
  sitio <- "PROT(R5);"
  probs_tmp <- list(KinaseX = setNames(list(1.0), sitio))
  rats_tmp  <- list(KinaseX = setNames(list(1.0), sitio))
  fdr_tmp   <- setNames(list(0.001), sitio)

  result <- run_makeks(probs = probs_tmp, rats = rats_tmp, fdr = fdr_tmp)
  expect_length(result[["KinaseX"]], 0)
})

test_that("exclude psites with residue (K", {
  sitio <- "PROT(K10);"
  probs_tmp <- list(KinaseX = setNames(list(1.0), sitio))
  rats_tmp  <- list(KinaseX = setNames(list(1.0), sitio))
  fdr_tmp   <- setNames(list(0.001), sitio)

  result <- run_makeks(probs = probs_tmp, rats = rats_tmp, fdr = fdr_tmp)
  expect_length(result[["KinaseX"]], 0)
})

test_that("exclude psites that are literally None", {
  sitio <- "None"
  probs_tmp <- list(KinaseX = setNames(list(1.0), sitio))
  rats_tmp  <- list(KinaseX = setNames(list(1.0), sitio))
  fdr_tmp   <- setNames(list(0.001), sitio)

  result <- run_makeks(probs = probs_tmp, rats = rats_tmp, fdr = fdr_tmp)
  expect_length(result[["KinaseX"]], 0)
})


# ==== TESTS: MULTI-SITE ====

test_that("multisite psite are divided into individual sites", {
  sitio <- "SITE1;SITE2;"
  probs_tmp <- list(KinaseX = setNames(list(1.0), sitio))
  rats_tmp  <- list(KinaseX = setNames(list(1.0), sitio))
  fdr_tmp   <- setNames(list(0.001), sitio)
  p_tmp     <- setNames(list(list(CompA = 0.01)), sitio)

  wb <- createWorkbook()
  addWorksheet(wb, "fold"); addWorksheet(wb, "pvalue")
  result <- makeKSlistOfKinaseDownstreamTargets(wb, rats_tmp, probs_tmp, p_tmp, fdr_tmp,
                                                 ratioThreshold = 0.5,
                                                 probabilityThreshold = 0.5)
  expect_true("SITE1" %in% result[["KinaseX"]])
  expect_true("SITE2" %in% result[["KinaseX"]])
})

test_that("multisite split does not include empty strings", {
  sitio <- "SITE1;"
  probs_tmp <- list(KinaseX = setNames(list(1.0), sitio))
  rats_tmp  <- list(KinaseX = setNames(list(1.0), sitio))
  fdr_tmp   <- setNames(list(0.001), sitio)
  p_tmp     <- setNames(list(list(CompA = 0.01)), sitio)

  wb <- createWorkbook()
  addWorksheet(wb, "fold"); addWorksheet(wb, "pvalue")
  result <- makeKSlistOfKinaseDownstreamTargets(wb, rats_tmp, probs_tmp, p_tmp, fdr_tmp,
                                                 ratioThreshold = 0.5,
                                                 probabilityThreshold = 0.5)
  expect_false("" %in% result[["KinaseX"]])
  expect_true("SITE1" %in% result[["KinaseX"]])
})


# ==== TESTS: CREATED SHEET ====

test_that("creates PutativeKinaseSubstrates sheet with the correct header", {
  wb <- createWorkbook()
  addWorksheet(wb, "fold"); addWorksheet(wb, "pvalue")
  makeKSlistOfKinaseDownstreamTargets(wb, ratios, prob_kinases, p_values, fdr_values,
                                       ratioThreshold = 0.5, probabilityThreshold = 0.5)

  expect_true("PutativeKinaseSubstrates" %in% names(wb))

  ws_data <- readWorkbook(wb, sheet = "PutativeKinaseSubstrates", colNames = FALSE)
  expect_equal(as.character(ws_data[1, 1]), "kinase")
  expect_equal(as.character(ws_data[1, 2]), "n")
  expect_equal(as.character(ws_data[1, 3]), "substrates")
})
