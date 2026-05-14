# test_getProbability.R

#unitary tests for getProbabilityofBeingKinaseSubs()

#run:
#   testthat::test_file("tests/R/test_getProbability.R")

source("helpers.R")

# ==== TEST DATA ====

ratios <- list(
  KinaseX = list("SITE1;" = 1.0, "SITE2;" = 0.0),
  KinaseY = list()  # 1 compound : empty dict
)

list_of_kinases <- list(
  MCF7 = "KinaseX;PRKCA;CDK7"  #expressed KinaseX, KinaseY not
)

wb <- createWorkbook()
addWorksheet(wb, "fold")
addWorksheet(wb, "pvalue")

probs <- getProbabilityofBeingKinaseSubs(wb, list_of_kinases, ratios, "MCF7")


# ==== TESTS: EXPRESSED KINASE ====

test_that("probability is 1 when ratio is max", {
  # SITE1: ratio KinaseX=1.0, maxRatio=1.0 : prob=1.0/1.0=1.0
  expect_equal(probs[["KinaseX"]][["SITE1;"]], 1.0, tolerance = 1e-9)
})

test_that("probability is 0 when ratio is 0 and maxRatio is 0", {
  # SITE2: ratio=0.0, maxRatio=0.0 : protection against dividing by 0: 0
  expect_equal(probs[["KinaseX"]][["SITE2;"]], 0.0, tolerance = 1e-9)
})

test_that("all probabilities are between 0 and 1", {
  for (kinase in names(probs)) {
    for (site in names(probs[[kinase]])) {
      val <- probs[[kinase]][[site]]
      expect_true(
        val >= 0 & val <= 1,
        label = paste("prob[", kinase, "][", site, "] =", val, "fuera de [0,1]")
      )
    }
  }
})


# ==== TESTS: NOT EXPRESSED KINASE ====

test_that("not expressed kinase has an empty probability dict", {
  # KinaseY is not in the kinase list of MCF7
  # (its ratios dict is empty too, so probabilities is empty)
  expect_equal(length(probs[["KinaseY"]]), 0)
})

test_that("Kinase with ratios but not expressed has probability 0 at all its sites", {
  wb2 <- createWorkbook()
  addWorksheet(wb2, "fold"); addWorksheet(wb2, "pvalue")

  ratios_tmp <- list(
    KinaseNoExpresada = list("SITE1;" = 0.8, "SITE2;" = 0.5)
  )
  lk_tmp <- list(MCF7 = "OtraKinasa;PRKCA")  # KinaseNoExpresada missing

  p <- getProbabilityofBeingKinaseSubs(wb2, lk_tmp, ratios_tmp, "MCF7")

  expect_equal(p[["KinaseNoExpresada"]][["SITE1;"]], 0.0, tolerance = 1e-9)
  expect_equal(p[["KinaseNoExpresada"]][["SITE2;"]], 0.0, tolerance = 1e-9)
})


# ==== TESTS: BORDERLINE CASE: maxRatio=0 ====

test_that("maxRatio=0 does not produce division by 0", {
  wb2 <- createWorkbook()
  addWorksheet(wb2, "fold"); addWorksheet(wb2, "pvalue")

  ratios_cero <- list(
    KinaseX = list("SITE1;" = 0.0, "SITE2;" = 0.0)  #all 0
  )
  lk_tmp <- list(MCF7 = "KinaseX")

  #cannot rise error
  expect_no_error({
    p <- getProbabilityofBeingKinaseSubs(wb2, lk_tmp, ratios_cero, "MCF7")
  })

  p <- getProbabilityofBeingKinaseSubs(wb2, lk_tmp, ratios_cero, "MCF7")
  expect_equal(p[["KinaseX"]][["SITE1;"]], 0.0, tolerance = 1e-9)
})

test_that("normalization is correct with several ratios", {
  # maxRatio of SITE1=max(0.8, 0.4)=0.8
  # prob KinaseA=0.8/0.8=1.0
  # prob KinaseB=0.4/0.8=0.5
  wb2 <- createWorkbook()
  addWorksheet(wb2, "fold"); addWorksheet(wb2, "pvalue")

  ratios_tmp <- list(
    KinaseA = list("SITE1;" = 0.8),
    KinaseB = list("SITE1;" = 0.4)
  )
  lk_tmp <- list(MCF7 = "KinaseA;KinaseB")

  p <- getProbabilityofBeingKinaseSubs(wb2, lk_tmp, ratios_tmp, "MCF7")

  expect_equal(p[["KinaseA"]][["SITE1;"]], 1.0, tolerance = 1e-9)
  expect_equal(p[["KinaseB"]][["SITE1;"]], 0.5, tolerance = 1e-9)
})


# ==== TESTS: CREATED SHEET ====

test_that("creates ProbOfBeingKinaseSubs sheet", {
  expect_true("ProbOfBeingKinaseSubs" %in% names(wb))
})
