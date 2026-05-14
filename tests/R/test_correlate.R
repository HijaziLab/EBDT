# test_correlate.R
# Unitary test for correlatePhosphoPeptideWithInhibitorSpecificity()
#run:
#   testthat::test_file("tests/R/test_correlate.R")

source("helpers.R")

# ==== TEST DATA ====

# KinaseX:2 compounds : calculated correlation, calculated ratio
# KinaseY:1 compound : no correlation nor ratio (needs>=2)
kinases_inhibited <- list(
  KinaseX = list(c("CompA", "CompB"), c(0.01, 0.05)),
  KinaseY = list(c("CompA"),          c(0.02))
)

# SITE1: fold<-1 and p<0.025 for both compounds : ratio=1.0
# SITE2: fold>-1 : ratio=0
f_values <- list(
  "SITE1;" = list(CompA = -2.5, CompB = -1.5),
  "SITE2;" = list(CompA =  0.8, CompB =  0.5)
)

p_values <- list(
  "SITE1;" = list(CompA = 0.010, CompB = 0.010),
  "SITE2;" = list(CompA = 0.800, CompB = 0.900)
)

#empty workbook (function needs one to write intermediate sheets)
wb <- createWorkbook()
addWorksheet(wb, "fold")   # necessary dummy sheets
addWorksheet(wb, "pvalue")

#run function
ratios <- correlatePhosphoPeptideWithInhibitorSpecificity(
  wb,
  kinasesInhibitedCompounds = kinases_inhibited,
  fValues                   = f_values,
  pValues                   = p_values,
  compoundsKuster           = list(),
  compoundsCellLine         = list()
)


# ==== TESTS: RATIO ====

test_that("ratio es 1.0 cuando todos los compuestos inhiben SITE1", {
  # CompA and CompB have fold<-1 and p<0.025 : both count : 2/2=1.0
  expect_equal(ratios[["KinaseX"]][["SITE1;"]], 1.0, tolerance = 1e-9)
})

test_that("ratio is 0 when fold does not exceed the threshold of -1", {
  # SITE2 has fold=0.8 and 0.5, none<-1 : ratio=0
  expect_equal(ratios[["KinaseX"]][["SITE2;"]], 0.0, tolerance = 1e-9)
})

test_that("ratio is 0 when p-value is not significant", {
  wb2 <- createWorkbook()
  addWorksheet(wb2, "fold"); addWorksheet(wb2, "pvalue")

  kinases_tmp <- list(KinaseX = list(c("CompA", "CompB"), c(0.01, 0.05)))
  f_tmp <- list("SITE1;" = list(CompA = -3.0, CompB = -2.0))  # fold ok
  p_tmp <- list("SITE1;" = list(CompA =  0.1, CompB =  0.5))  # p-value not significant

  r <- correlatePhosphoPeptideWithInhibitorSpecificity(wb2, kinases_tmp, f_tmp, p_tmp, list(), list())
  expect_equal(r[["KinaseX"]][["SITE1;"]], 0.0, tolerance = 1e-9)
})

test_that("parcial ratio when only one compound meets criteria", {
  wb2 <- createWorkbook()
  addWorksheet(wb2, "fold"); addWorksheet(wb2, "pvalue")

  kinases_tmp <- list(KinaseX = list(c("CompA", "CompB"), c(0.01, 0.05)))
  f_tmp <- list("SITE1;" = list(CompA = -2.0, CompB =  0.5))  # only CompA<-1
  p_tmp <- list("SITE1;" = list(CompA =  0.01, CompB = 0.01)) # both significant

  r <- correlatePhosphoPeptideWithInhibitorSpecificity(wb2, kinases_tmp, f_tmp, p_tmp, list(), list())
  expect_equal(r[["KinaseX"]][["SITE1;"]], 0.5, tolerance = 1e-9)
})

test_that("KinaseY with just one compound has an empty dict ratio", {
  #function requieres>=2 compounds to calculate ratio
  expect_equal(length(ratios[["KinaseY"]]), 0)
})


# ==== TESTS: CORRELATION ====

test_that("correlation is +1 when inhibition and fold increase together", {
  # kinaseValues=[0.01, 0.05], foldValues=[-2.5, -1.5]: both increase : r=1.0
  ws_corr <- readWorkbook(wb, sheet = "corrPPwithKinases", colNames = FALSE)

  #search row of SITE1 and KinaseX column
  header_row  <- as.character(ws_corr[1, ])
  kinaseX_col <- which(header_row == "KinaseX")
  site1_row   <- which(as.character(ws_corr[, 1]) == "SITE1;")

  corr_val <- as.numeric(ws_corr[site1_row, kinaseX_col])
  expect_equal(corr_val, 1.0, tolerance = 1e-9)
})

test_that("correlation is -1 when inhibition rises and fold drops", {
  # kinaseValues=[0.01,0.05], foldValues=[0.8,0.5]: x rises, and drops : r=-1.0
  wb2 <- createWorkbook()
  addWorksheet(wb2, "fold"); addWorksheet(wb2, "pvalue")

  kinases_tmp <- list(KinaseX = list(c("CompA", "CompB"), c(0.01, 0.05)))
  f_tmp <- list("SITE2;" = list(CompA = 0.8, CompB = 0.5))
  p_tmp <- list("SITE2;" = list(CompA = 0.5, CompB = 0.5))

  correlatePhosphoPeptideWithInhibitorSpecificity(wb2, kinases_tmp, f_tmp, p_tmp, list(), list())

  ws_corr    <- readWorkbook(wb2, sheet = "corrPPwithKinases", colNames = FALSE)
  header_row <- as.character(ws_corr[1, ])
  kinaseX_col <- which(header_row == "KinaseX")
  site2_row   <- which(as.character(ws_corr[, 1]) == "SITE2;")

  corr_val <- as.numeric(ws_corr[site2_row, kinaseX_col])
  expect_equal(corr_val, -1.0, tolerance = 1e-9)
})

test_that("KinaseY does not appear in corrPPwithKinases sheet", {
  # KinaseY has <2 compounds : correlation is not calculated: it does not appear in the sheet
  ws_corr    <- readWorkbook(wb, sheet = "corrPPwithKinases", colNames = FALSE)
  header_row <- as.character(ws_corr[1, ])
  expect_false("KinaseY" %in% header_row)
})


# ==== TESTS: CREATED SHEETS ====

test_that("creates corrPPwithKinases sheet", {
  expect_true("corrPPwithKinases" %in% names(wb))
})

test_that("creates ratioSigPPoverSignInhSpeci sheet", {
  expect_true("ratioSigPPoverSignInhSpeci" %in% names(wb))
})
