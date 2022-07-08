gold_path <- system.file("testdata", package = "countland", mustWork = TRUE)
gold.data <- Seurat::Read10X(data.dir = gold_path)

test_that("function Subsample, works correctly",{
  C <- countland(gold.data)

  C <- Subsample(C,gene_counts=10)
  expect_equal(max(apply(C@subsample,1,sum)),10)

  C <- Subsample(C,cell_counts=10)
  expect_equal(min(apply(C@subsample,2,sum)),10)

  C <- Subsample(C,cell_counts=10,gene_counts=10)
  expect_true(max(apply(C@subsample,1,sum)) <= 10)
  expect_equal(min(apply(C@subsample,2,sum)),10)
})
