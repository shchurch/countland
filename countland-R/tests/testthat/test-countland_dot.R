gold_path <- system.file("testdata", package = "countland", mustWork = TRUE)
gold.data <- Seurat::Read10X(data.dir = gold_path)

test_that("function Dot, returns correct values", {
  C <- countland(gold.data[1:10,1:3])
  C <- Dot(C)

  gd <- C@counts

  expect_equal(C@dots[1,1],as.numeric(gd[,1] %*% gd[,1]))
  expect_equal(C@dots[1,2],as.numeric(gd[,1] %*% gd[,2]))
})
