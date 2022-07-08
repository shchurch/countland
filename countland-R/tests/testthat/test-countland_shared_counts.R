gold_path <- system.file("testdata", package = "countland", mustWork = TRUE)
gold.data <- Seurat::Read10X(data.dir = gold_path)

test_that("function SharedCounts, returns correct value", {
  C <- countland(gold.data[1:10,1:10])
  C <- SharedCounts(C,n_clusters=2,n_cells=10,subsample=F)
  expect_equal(C@sharedcounts[1,3],sum(pmin(C@counts[1,],C@counts[3,])))
  expect_equal(C@sharedcounts[2,3],sum(pmin(C@counts[2,],C@counts[3,])))
})

test_that("function PlotSharedCounts, returns plot object", {
  C <- countland(gold.data)
  C <- SharedCounts(C,n_clusters=2,n_cells=10,subsample=F)
  C@cluster_labels <- rep(c(1,2),length(C@names_cells))
  expect_true(ggplot2::is.ggplot(PlotSharedCounts(C)))
})
