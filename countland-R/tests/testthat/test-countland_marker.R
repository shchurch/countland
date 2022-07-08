test_that("function RankMarkerGenes, option method='prop-zero', returns correct values", {
  counts <- t(matrix(rep(c(1,1,1,1,1,1),8),6))
  counts <- rbind(counts,c(1,1,1,0,0,0))
  counts <- rbind(counts,c(0,0,0,1,1,1))
  names <- paste0("gene",seq_len(10))
  clusters <- c(1,1,1,2,2,2)

  C <- new("countland")
  C@cluster_labels <- clusters
  C@counts <- as(counts,"dgCMatrix")
  C@names_genes <- names

  C <- RankMarkerGenes(C,method="prop-zero")

  expect_equal(C@marker_genes[1:10,]$names[1],"gene9")
  expect_equal(C@marker_genes[11:20,]$names[1],"gene10")
})

test_that("function RankMarkerGenes, option method='rank-sums', returns correct values", {
  counts <- t(matrix(rep(c(1,1,1,1,1,1),8),6))
  counts <- rbind(counts,c(10,10,10,1,1,1))
  counts <- rbind(counts,c(1,1,1,10,10,10))
  names <- paste0("gene",seq_len(10))
  clusters <- c(1,1,1,2,2,2)

  C <- new("countland")
  C@cluster_labels <- clusters
  C@counts <- as(counts,"dgCMatrix")
  C@names_genes <- names

  C <- RankMarkerGenes(C,method="rank-sums")

  expect_equal(C@marker_genes[1:10,]$names[1],"gene9")
  expect_equal(C@marker_genes[11:20,]$names[1],"gene10")
})

mat <- rbind(c(0.01,0.01,0.01),
             c(0.01,0.02,0.01),
             c(100,100,100))

test_that("function PlotMarker, returns two objects", {
  C <- new("countland")
  C@counts <- as(mat,"dgCMatrix")
  C@embedding <- mat
  C@cluster_labels <- c(1,1,2)
  expect_length(PlotMarker(C,1),2)
})
