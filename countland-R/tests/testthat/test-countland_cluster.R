mat <- rbind(c(0.01,0.01,0.01),
             c(0.01,0.02,0.01),
             c(100,100,100))

test_that("function Cluster, returns correct values",{
  C <- new("countland")
  C@embedding <- mat
  C@verbose=TRUE
  C <- Cluster(C,n_clusters=2)
  C@cluster_labels

  expect_true(C@cluster_labels[1] == C@cluster_labels[2])
  expect_true(C@cluster_labels[1] != C@cluster_labels[3])
})

test_that("function PlotEigengap, returns plot object", {
  C <- new("countland")
  C@embedding <- mat
  C@cluster_labels <- c(1,1,2)
  expect_true(ggplot2::is.ggplot(PlotEmbedding(C)))
})
