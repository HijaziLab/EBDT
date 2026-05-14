# run_tests.R: Script to run R tests

#All tests (including the slow tests, with real data):
#   source("tests/R/run_tests.R")
#
#    Only one file:
#   testthat::test_file("tests/R/test_correlate.R")
#
# Requirements: install.packages(c("testthat", "openxlsx"))

library(testthat)

#Change the directory to the tests so helpers.R can use relative routes.
old_wd <- getwd()
setwd(dirname(sys.frame(1)$ofile))

cat("\n========================================\n")
cat("  Unitary test and integration EBDT \n")
cat("========================================\n\n")

#run all the tests (except golden by default: they are really slow)
cat("--- Unitary tests + integration ---\n")
test_results <- test_dir(
  ".",
  filter    = "^test_(getDataForKusterCompounds|correlate|getProbability|makeKSlist|createNodeEdges|integration)",
  reporter  = "progress"
)

cat("\n\n--- Regresion tests golden (slow) ---\n")
cat("To run: testthat::test_file('tests/R/test_golden.R')\n\n")

setwd(old_wd)
