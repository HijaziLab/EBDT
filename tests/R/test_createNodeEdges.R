# test_createNodeEdges.R

# Unitary tests for createNodeEdges()

# Run:
#   testthat::test_file("tests/R/test_createNodeEdges.R")

source("helpers.R")

# ==== TEST DATA ====

# KinaseA and KinaseB share SITE1 and SITE3 : appear in nodes.edges with weigh 2
# KinaseC only has SITE4 : it shares anything
kinase_substrates <- list(
  KinaseA = c("SITE1", "SITE2", "SITE3"),
  KinaseB = c("SITE1", "SITE3"),
  KinaseC = c("SITE4")
)

wb <- createWorkbook()
addWorksheet(wb, "fold"); addWorksheet(wb, "pvalue")
createNodeEdges(wb, kinase_substrates)


# ==== TESTS: nodes.edges ====

test_that("pair with shared substrates appears in nodes.edges", {
  ws <- readWorkbook(wb, sheet = "nodes.edges", colNames = FALSE)
  edges <- as.character(ws[-1, 1])  #jump header
  expect_true("KinaseA.KinaseB" %in% edges)
})

test_that("pair without shared substrates does not appear in nodes.edges", {
  ws <- readWorkbook(wb, sheet = "nodes.edges", colNames = FALSE)
  edges <- as.character(ws[-1, 1])
  expect_false("KinaseA.KinaseC" %in% edges)
  expect_false("KinaseB.KinaseC" %in% edges)
})

test_that("weight of KinaseA-KinaseB is 2 (SITE1 and SITE3 shared)", {
  ws <- readWorkbook(wb, sheet = "nodes.edges", colNames = FALSE)
  fila <- ws[as.character(ws[, 1]) == "KinaseA.KinaseB", ]
  expect_equal(as.numeric(fila[[2]]), 2)
})

test_that("shared substrates appear in subs column", {
  ws <- readWorkbook(wb, sheet = "nodes.edges", colNames = FALSE)
  fila <- ws[as.character(ws[, 1]) == "KinaseA.KinaseB", ]
  subs <- strsplit(as.character(fila[[3]]), ";")[[1]]
  expect_true("SITE1" %in% subs)
  expect_true("SITE3" %in% subs)
})

test_that("nodes.edges has header: edge, weight, subs", {
  ws <- readWorkbook(wb, sheet = "nodes.edges", colNames = FALSE)
  expect_equal(as.character(ws[1, 1]), "edge")
  expect_equal(as.character(ws[1, 2]), "weight")
  expect_equal(as.character(ws[1, 3]), "subs")
})


# ==== TESTS: tableEdgeSubs ====

test_that("tableEdgeSubs has correct dimensions (n+1 x n+1)", {
  ws <- readWorkbook(wb, sheet = "tableEdgeSubs", colNames = FALSE)
  n <- length(kinase_substrates)
  expect_equal(nrow(ws), n + 1)
  expect_equal(ncol(ws), n + 1)
})

test_that("tableEdgeSubs is simetrycal", {
  ws <- readWorkbook(wb, sheet = "tableEdgeSubs", colNames = FALSE)
  n <- length(kinase_substrates)

  for (i in 2:(n + 1)) {
    for (j in 2:(n + 1)) {
      val_ij <- as.character(ws[i, j] %||% "")
      val_ji <- as.character(ws[j, i] %||% "")
      expect_equal(val_ij, val_ji,
                   label = paste0("asymmetry en [", i, ",", j, "]"))
    }
  }
})


# ==== TESTS: CREATED SHEETS ====

test_that("creates tableEdgeSubs sheet", {
  expect_true("tableEdgeSubs" %in% names(wb))
})

test_that("creates nodes.edges sheet", {
  expect_true("nodes.edges" %in% names(wb))
})


# ==== TEST: BORDERLINE CASE: no shared substrates ====

test_that("nodes.edges only has header if there is no shared substrates", {
  wb2 <- createWorkbook()
  addWorksheet(wb2, "fold"); addWorksheet(wb2, "pvalue")

  substrates_disjuntos <- list(
    KinaseA = c("SITE1"),
    KinaseB = c("SITE2")   #no overlap
  )
  createNodeEdges(wb2, substrates_disjuntos)

  ws <- readWorkbook(wb2, sheet = "nodes.edges", colNames = FALSE)
  #only header row, no data
  data_rows <- ws[-1, ]
  data_rows_con_datos <- data_rows[!is.na(data_rows[, 1]), ]
  expect_equal(nrow(data_rows_con_datos), 0)
})
