gold_path <- system.file("testdata", package = "countland", mustWork = TRUE)
gold.data <- Seurat::Read10X(data.dir = gold_path)

test_that("function SubsetGenes, verify works correctly", {
  C <- countland(gold.data)
  C <- SubsetGenes(C,seq_len(25))
  expect_equal(nrow(C@counts), 25)
  expect_equal(length(C@names_genes), 25)
})

test_that("function SubsetCells, verify works correctly", {
  C <- countland(gold.data)
  C <- SubsetCells(C,seq_len(25))
  expect_equal(ncol(C@counts), 25)
  expect_equal(length(C@names_cells), 25)
})

test_that("function RestoreCounts, verify works correctly", {
  C <- countland(gold.data)

  C@counts <- C@counts[1:25,1:25]

  expect_true(nrow(C@counts) != nrow(C@raw_counts))
  expect_true(ncol(C@counts) != ncol(C@raw_counts))

  C <- RestoreCounts(C)

  expect_equal(nrow(C@counts), nrow(C@raw_counts))
  expect_equal(ncol(C@counts), ncol(C@raw_counts))

  expect_equal(nrow(C@counts), length(C@names_genes))
  expect_equal(ncol(C@counts), length(C@names_cells))
})
