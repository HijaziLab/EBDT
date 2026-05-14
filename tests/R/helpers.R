# helpers.R : Auxiliary functions shared by all .R 
#
#loaded at the startup of each test with: source("helpers.R")
#It does not have tests, just utilities 

library(openxlsx)
library(testthat)

# ==== CONFIGURABLE ROUTES ====
# change if you move the files

TESTS_DIR    <- normalizePath(file.path(getwd(), ".."), mustWork = FALSE)
CODE_DIR     <- normalizePath(file.path(TESTS_DIR, ".."), mustWork = FALSE)
FIXTURES_DIR <- file.path(TESTS_DIR, "fixtures")
GOLDEN_DIR   <- file.path(TESTS_DIR, "golden")

# ==== SHIM: compatibility keepVBA ====
# openxlsx does not accept keepVBA in all its editions. This shim ignores it
# so ebdt.R (that uses keepVBA=TRUE) works properly
.loadWorkbook_orig <- openxlsx::loadWorkbook
loadWorkbook <- function(xlsxFile, ..., keepVBA = NULL, isUnzipped = FALSE) {
  .loadWorkbook_orig(xlsxFile, isUnzipped = isUnzipped)
}

#load the algorithm`s code
source(file.path(CODE_DIR, "ebdt.R"))


# ==== make_workbook() ====

# Creates an openxlsx workbook in memory with sheets 'fold' and 'pvalue'

# Parametres:
#   fold_matrix: character array Row1=header.
#                   Format: c("sh.index.sites", "FDR", "MCF7.CompA.fold", ...)
#   pvalue_matrix: same structure for p-values

# Returns: Workbook object


make_workbook <- function(fold_matrix, pvalue_matrix) {
  wb <- createWorkbook()
  addWorksheet(wb, "fold")
  addWorksheet(wb, "pvalue")

  writeData(wb, "fold",   fold_matrix,   startRow = 1, startCol = 1, colNames = FALSE)
  writeData(wb, "pvalue", pvalue_matrix, startRow = 1, startCol = 1, colNames = FALSE)

  #save temporary file and recharge for compatibility guarantee 
  #with readWorkbook() exactly the same to the production code 
  tmp <- tempfile(fileext = ".xlsx")
  on.exit(unlink(tmp), add = TRUE)
  saveWorkbook(wb, tmp, overwrite = TRUE)

  return(loadWorkbook(tmp))
}


# ==== make_run_dir() ====
#
#creates a temporary directory with all the neccessary to run the pipeline
# Returns the route to the created directory

make_run_dir <- function(xlsm_source = file.path(FIXTURES_DIR, "MCF7_mini.xlsm"),
                         cell_line   = "MCF7") {
  run_dir <- tempfile()
  dir.create(run_dir)

  # Copy the Excel named correctly
  file.copy(xlsm_source, file.path(run_dir, paste0(cell_line, ".xlsm")))

  # Copy requiered CSVss
  req_src <- file.path(FIXTURES_DIR, "requiredData")
  req_dst <- file.path(run_dir, "requiredData")
  dir.create(req_dst)
  file.copy(list.files(req_src, full.names = TRUE), req_dst)

  return(run_dir)
}
