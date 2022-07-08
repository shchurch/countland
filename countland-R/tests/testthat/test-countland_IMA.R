gold_path <- system.file("testdata", package = "countland", mustWork = TRUE)
gold.data <- Seurat::Read10X(data.dir = gold_path)

test_that("function RunIMA, returns correct object", {
  C <- countland(gold.data)
  C <- RunIMA(C,features=10,u_bounds=c(10,10),subsample=F)
  expect_length(C@matrixU,1000)
  expect_s4_class(C@matrixU,"dgCMatrix")
  expect_s4_class(C@matrixV,"dgCMatrix")
  expect_s4_class(C@matrixLambda,"dgCMatrix")
})

test_that("function PlotIMA, returns plot object", {
  C <- countland(gold.data)
  C <- RunIMA(C,features=10,u_bounds=c(10,10),subsample=F)
  C@cluster_labels <- rep(c(1,2),length(C@names_cells))
  expect_true(ggplot2::is.ggplot(PlotIMA(C,subsample=F)))
})

test_that("function PlotEigengap, returns plot object", {
  C <- countland(gold.data)
  expect_true(ggplot2::is.ggplot(PlotIMAElbow(C,max_features=2,u_bounds=c(10,10),subsample=F)))
})
