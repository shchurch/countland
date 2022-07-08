gold_path <- system.file("testdata", package = "countland", mustWork = TRUE)
gold.data <- Seurat::Read10X(data.dir = gold_path)

test_that("function countland, verify initialize correctly", {
  C <- countland(gold.data,remove_empty=F)
  expect_equal(nrow(C@counts), nrow(gold.data))
  expect_equal(ncol(C@counts), ncol(gold.data))
  expect_equal(nrow(C@counts), length(C@names_genes))
  expect_equal(ncol(C@counts), length(C@names_cells))

  expect_equal(C@raw_counts, C@counts)

})

test_that("function countland, verify remove_empty=TRUE works", {
  C <- countland(gold.data)

  cell_indices <- which(diff(gold.data@p)>0)
  gene_indices <- sort(unique(gold.data@i))+1
  gold.data <- gold.data[gene_indices,]
  gold.data <- gold.data[,cell_indices]

  expect_equal(nrow(C@counts), nrow(gold.data))
  expect_equal(ncol(C@counts), ncol(gold.data))

  expect_equal(nrow(C@counts), length(C@names_genes))
  expect_equal(ncol(C@counts), length(C@names_cells))

  expect_equal(C@raw_counts, C@counts)

})

