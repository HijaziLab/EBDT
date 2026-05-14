# test_integration.R

# Integration tests: run the full pipeline with mini data (fixtures)
# and verify that the result has the correct structure.

# Parametrized over MCF7, HL60 and NTERA2 via a helper function + loop.

# Run:
#   testthat::test_file("tests/R/test_integration.R")

source("helpers.R")

HOJAS_OUTPUT_ESPERADAS <- c(
  "pvalue.select", "fold.select",
  "corrPPwithKinases", "ratioSigPPoverSignInhSpeci",
  "ProbOfBeingKinaseSubs", "PutativeKinaseSubstrates",
  "tableEdgeSubs", "nodes.edges"
)


# ==== HELPER: EJECUTAR TESTS DE INTEGRACION PARA UN WORKBOOK ====

run_mini_tests <- function(cell_line, wb) {

  test_that(paste0(cell_line, ": All expected output sheets are generated"), {
    hojas_actuales <- names(wb)
    for (hoja in HOJAS_OUTPUT_ESPERADAS) {
      expect_true(hoja %in% hojas_actuales,
                  label = paste(cell_line, "-missing sheet:", hoja))
    }
  })

  test_that(paste0(cell_line, ": input sheets are preserved"), {
    expect_true("fold"   %in% names(wb))
    expect_true("pvalue" %in% names(wb))
  })

  test_that(paste0(cell_line, ": PutativeKinaseSubstrates has correct header"), {
    ws <- readWorkbook(wb, sheet = "PutativeKinaseSubstrates", colNames = FALSE)
    expect_equal(as.character(ws[1, 1]), "kinase")
    expect_equal(as.character(ws[1, 2]), "n")
    expect_equal(as.character(ws[1, 3]), "substrates")
  })

  test_that(paste0(cell_line, ": PutativeKinaseSubstrates has at least one row by kinase"), {
    ws <- readWorkbook(wb, sheet = "PutativeKinaseSubstrates", colNames = FALSE)
    filas_datos <- ws[-1, ]
    filas_con_datos <- filas_datos[!is.na(filas_datos[, 1]), ]
    expect_true(nrow(filas_con_datos) >= 1)
    expect_true(nrow(filas_con_datos) <= 163)
  })

  test_that(paste0(cell_line, ": column n of PutativeKinaseSubstrates contains no negative integers"), {
    ws <- readWorkbook(wb, sheet = "PutativeKinaseSubstrates", colNames = FALSE)
    valores_n <- suppressWarnings(as.numeric(as.character(ws[-1, 2])))
    valores_n <- valores_n[!is.na(valores_n)]
    expect_true(all(valores_n >= 0))
    expect_true(all(valores_n == floor(valores_n)))
  })

  test_that(paste0(cell_line, ": nodes.edges has a correct header"), {
    ws <- readWorkbook(wb, sheet = "nodes.edges", colNames = FALSE)
    expect_equal(as.character(ws[1, 1]), "edge")
    expect_equal(as.character(ws[1, 2]), "weight")
    expect_equal(as.character(ws[1, 3]), "subs")
  })

  test_that(paste0(cell_line, ": all correlations are in [-1, 1]"), {
    ws <- readWorkbook(wb, sheet = "corrPPwithKinases", colNames = FALSE)
    for (i in 2:min(6, nrow(ws))) {
      for (j in 2:ncol(ws)) {
        val <- suppressWarnings(as.numeric(as.character(ws[i, j])))
        if (!is.na(val)) {
          expect_true(val >= -1.01 & val <= 1.01,
                      label = paste(cell_line, "-out of range correlation:", val))
        }
      }
    }
  })

  test_that(paste0(cell_line, ":all ratios are between [0, 1]"), {
    ws <- readWorkbook(wb, sheet = "ratioSigPPoverSignInhSpeci", colNames = FALSE)
    for (i in 2:min(6, nrow(ws))) {
      for (j in 2:ncol(ws)) {
        val <- suppressWarnings(as.numeric(as.character(ws[i, j])))
        if (!is.na(val)) {
          expect_true(val >= 0 & val <= 1,
                      label = paste(cell_line, "-out of range ratio:", val))
        }
      }
    }
  })

  test_that(paste0(cell_line, ":all probabilities are between [0, 1]"), {
    ws <- readWorkbook(wb, sheet = "ProbOfBeingKinaseSubs", colNames = FALSE)
    for (i in 2:min(6, nrow(ws))) {
      for (j in 2:ncol(ws)) {
        val <- suppressWarnings(as.numeric(as.character(ws[i, j])))
        if (!is.na(val)) {
          expect_true(val >= 0 & val <= 1,
                      label = paste(cell_line, "-out of range probability:", val))
        }
      }
    }
  })
}


# ==== LOOP: RUN FOR EACH CELL LINE ====

old_wd <- getwd()

for (cell_line in c("MCF7", "HL60", "NTERA2")) {

  run_dir <- make_run_dir(
    xlsm_source = file.path(FIXTURES_DIR, paste0(cell_line, "_mini.xlsm")),
    cell_line   = cell_line
  )

  tryCatch({
    setwd(run_dir)
    suppressWarnings(
      GetExpectancyOfBeingDownstreamTarget(
        inhibitionThreshold  = 0.1,
        ratioThreshold       = 0.2,
        probabilityThreshold = 0.5,
        cellLinesFiles       = c(paste0(cell_line, ".xlsm"))
      )
    )
  }, finally = {
    setwd(old_wd)
  })

  wb <- loadWorkbook(file.path(run_dir, paste0(cell_line, ".xlsm")), keepVBA = FALSE)
  run_mini_tests(cell_line, wb)
}
