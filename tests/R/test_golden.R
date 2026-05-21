# test_golden.R : Regression tests (golden files)

# Run the full pipeline with all the real data and verifies that
# results are identical to the reference files in tests/golden/.

# these tests are slow (few minutes), run:
#   testthat::test_file("tests/R/test_golden.R")
#
#To exclude when running all the tests:
#   testthat::test_dir("tests/R/", filter = "^(?!golden)")

source("helpers.R")

# ==== PIPELINE PARAMETRES ====

INHIBITION_THRESHOLD   <- 0.5
RATIO_THRESHOLD        <- 0.5
PROBABILITY_THRESHOLD  <- 0.75


# ==== HELPER: READ GOLDEN CSV ====

read_golden <- function(cell_line, sheet_name) {
  path <- file.path(GOLDEN_DIR, paste0(cell_line, "_", sheet_name, ".csv"))
  read.csv(path, header = FALSE, stringsAsFactors = FALSE)
}


# ==== HELPER: PREPARE AND EXECUTE A PIPELINE DIRECTORY ====

make_golden_run_dir <- function(cell_line) {
  run_dir <- tempfile(); dir.create(run_dir)
  file.copy(file.path(CODE_DIR, paste0(cell_line, ".xlsm")),
            file.path(run_dir, paste0(cell_line, ".xlsm")))
  req_dst <- file.path(run_dir, "requiredData"); dir.create(req_dst)
  file.copy(list.files(file.path(CODE_DIR, "requiredData"), full.names = TRUE), req_dst)
  run_dir
}

run_pipeline <- function(run_dir, cell_line) {
  old_wd <- getwd()
  tryCatch({
    setwd(run_dir)
    GetExpectancyOfBeingDownstreamTarget(
      inhibitionThreshold  = INHIBITION_THRESHOLD,
      ratioThreshold       = RATIO_THRESHOLD,
      probabilityThreshold = PROBABILITY_THRESHOLD,
      cellLinesFiles       = c(paste0(cell_line, ".xlsm"))
    )
  }, finally = setwd(old_wd))
}


# ==== RUNS PIPELINE WITH REAL DATA ====

run_dir_mcf7   <- make_golden_run_dir("MCF7");   run_pipeline(run_dir_mcf7,   "MCF7")
run_dir_hl60   <- make_golden_run_dir("HL60");   run_pipeline(run_dir_hl60,   "HL60")
run_dir_ntera2 <- make_golden_run_dir("NTERA2"); run_pipeline(run_dir_ntera2, "NTERA2")

# Load resulting workbooks
wb_mcf7   <- loadWorkbook(file.path(run_dir_mcf7,   "MCF7.xlsm"),   keepVBA = FALSE)
wb_hl60   <- loadWorkbook(file.path(run_dir_hl60,   "HL60.xlsm"),   keepVBA = FALSE)
wb_ntera2 <- loadWorkbook(file.path(run_dir_ntera2, "NTERA2.xlsm"), keepVBA = FALSE)


#==== TESTS: PutativeKinaseSubstrates ====

test_that("MCF7: number of Kinases with susbtrates that match golden", {
  golden <- read_golden("MCF7", "PutativeKinaseSubstrates")
  actual <- readWorkbook(wb_mcf7, sheet = "PutativeKinaseSubstrates", colNames = FALSE)

  golden_n <- sum(suppressWarnings(as.numeric(as.character(golden[-1, 2]))) > 0,
                  na.rm = TRUE)
  actual_n <- sum(suppressWarnings(as.numeric(as.character(actual[-1, 2]))) > 0,
                  na.rm = TRUE)

  expect_equal(actual_n, golden_n,
               label = paste("MCF7: kinases with substrates:", actual_n,
                             "vs", golden_n, "(golden)"))
})

test_that("HL60: number of kinases with substrates that match golden", {
  golden <- read_golden("HL60", "PutativeKinaseSubstrates")
  actual <- readWorkbook(wb_hl60, sheet = "PutativeKinaseSubstrates", colNames = FALSE)

  golden_n <- sum(suppressWarnings(as.numeric(as.character(golden[-1, 2]))) > 0,
                  na.rm = TRUE)
  actual_n <- sum(suppressWarnings(as.numeric(as.character(actual[-1, 2]))) > 0,
                  na.rm = TRUE)

  expect_equal(actual_n, golden_n)
})


# ==== TESTS: nodes.edges ====

test_that("MCF7: number of kinase pairs with edges that match golden", {
  golden <- read_golden("MCF7", "nodes.edges")
  actual <- readWorkbook(wb_mcf7, sheet = "nodes.edges", colNames = FALSE)

  golden_pares <- sum(!is.na(golden[-1, 1]) & golden[-1, 1] != "")
  actual_pares <- sum(!is.na(actual[-1, 1]) & actual[-1, 1] != "")

  expect_equal(actual_pares, golden_pares,
               label = paste("MCF7 pairs:", actual_pares, "vs", golden_pares, "(golden)"))
})

test_that("HL60: number of kinase pairs with edges that match golden", {
  golden <- read_golden("HL60", "nodes.edges")
  actual <- readWorkbook(wb_hl60, sheet = "nodes.edges", colNames = FALSE)

  golden_pares <- sum(!is.na(golden[-1, 1]) & golden[-1, 1] != "")
  actual_pares <- sum(!is.na(actual[-1, 1]) & actual[-1, 1] != "")

  expect_equal(actual_pares, golden_pares)
})


# ==== TESTS: CORRELATIONS (first row of MCF7) ====

test_that("MCF7: correlations first row that match golden (tolerance 1e-6)", {
  golden <- read_golden("MCF7", "corrPPwithKinases")
  actual <- readWorkbook(wb_mcf7, sheet = "corrPPwithKinases", colNames = FALSE)

  # Compare row 2 (first psite)) column to column
  golden_row <- suppressWarnings(as.numeric(as.character(golden[2, -1])))
  actual_row <- suppressWarnings(as.numeric(as.character(actual[2, -1])))

  expect_equal(length(actual_row), length(golden_row),
               label = "Column number different from corrPPwithKinases")

  for (j in seq_along(golden_row)) {
    if (!is.na(golden_row[j])) {
      expect_equal(actual_row[j], golden_row[j], tolerance = 1e-6,
                   label = paste("corrPPwithKinases col", j))
    }
  }
})


# ==== TESTS: TOP KINASE SUBSTRATES ====

test_that("MCF7: substrates of the kinase with more substrates that match golden", {
  golden <- read_golden("MCF7", "PutativeKinaseSubstrates")
  actual <- readWorkbook(wb_mcf7, sheet = "PutativeKinaseSubstrates", colNames = FALSE)

  # Find top kinase in golden
  golden_n    <- suppressWarnings(as.numeric(as.character(golden[-1, 2])))
  top_idx     <- which.max(golden_n)
  top_kinase  <- as.character(golden[top_idx + 1, 1])
  top_n       <- golden_n[top_idx]

  #Search in current
  actual_fila <- actual[as.character(actual[, 1]) == top_kinase, ]
  actual_n    <- suppressWarnings(as.numeric(as.character(actual_fila[1, 2])))

  expect_equal(actual_n, top_n,
               label = paste("Substrates of", top_kinase, ":", actual_n,
                             "vs", top_n, "(golden)"))
})


# ==== TESTS: NTERA2 ====

test_that("NTERA2: The number of kinases with substrates matches golden", {
  golden <- read_golden("NTERA2", "PutativeKinaseSubstrates")
  actual <- readWorkbook(wb_ntera2, sheet = "PutativeKinaseSubstrates", colNames = FALSE)

  golden_n <- sum(suppressWarnings(as.numeric(as.character(golden[-1, 2]))) > 0,
                  na.rm = TRUE)
  actual_n <- sum(suppressWarnings(as.numeric(as.character(actual[-1, 2]))) > 0,
                  na.rm = TRUE)

  expect_equal(actual_n, golden_n,
               label = paste("NTERA2: kinases with substrates:", actual_n,
                             "vs", golden_n, "(golden)"))
})

test_that("NTERA2: number of kinase pairs with edges that match golden", {
  golden <- read_golden("NTERA2", "nodes.edges")
  actual <- readWorkbook(wb_ntera2, sheet = "nodes.edges", colNames = FALSE)

  golden_pares <- sum(!is.na(golden[-1, 1]) & golden[-1, 1] != "")
  actual_pares <- sum(!is.na(actual[-1, 1]) & actual[-1, 1] != "")

  expect_equal(actual_pares, golden_pares,
               label = paste("NTERA2 pairs:", actual_pares, "vs", golden_pares, "(golden)"))
})

test_that("NTERA2: first row correlations that match golden (tolerancia 1e-6)", {
  golden <- read_golden("NTERA2", "corrPPwithKinases")
  actual <- readWorkbook(wb_ntera2, sheet = "corrPPwithKinases", colNames = FALSE)

  golden_row <- suppressWarnings(as.numeric(as.character(golden[2, -1])))
  actual_row <- suppressWarnings(as.numeric(as.character(actual[2, -1])))

  expect_equal(length(actual_row), length(golden_row),
               label = "different column number in NTERA2 corrPPwithKinases")

  for (j in seq_along(golden_row)) {
    if (!is.na(golden_row[j])) {
      expect_equal(actual_row[j], golden_row[j], tolerance = 1e-6,
                   label = paste("NTERA2 corrPPwithKinases col", j))
    }
  }
})

test_that("NTERA2: substrates of the kinase with more substrates that match golden", {
  golden <- read_golden("NTERA2", "PutativeKinaseSubstrates")
  actual <- readWorkbook(wb_ntera2, sheet = "PutativeKinaseSubstrates", colNames = FALSE)

  golden_n    <- suppressWarnings(as.numeric(as.character(golden[-1, 2])))
  top_idx     <- which.max(golden_n)
  top_kinase  <- as.character(golden[top_idx + 1, 1])
  top_n       <- golden_n[top_idx]

  actual_fila <- actual[as.character(actual[, 1]) == top_kinase, ]
  actual_n    <- suppressWarnings(as.numeric(as.character(actual_fila[1, 2])))

  expect_equal(actual_n, top_n,
               label = paste("NTERA2 substrates of", top_kinase, ":", actual_n,
                             "vs", top_n, "(golden)"))
})
