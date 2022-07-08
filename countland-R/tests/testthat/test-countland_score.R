gold_path <- system.file("testdata", package = "countland", mustWork = TRUE)
gold.data <- Seurat::Read10X(data.dir = gold_path)

test_that("function ScoreCells, returns correct values", {
  C <- countland(gold.data)
  C <- ScoreCells(C)

  cell_index <- 13
  cell_scores <- C@cell_scores[cell_index,]
  cell_name <- C@cell_scores$names[cell_index]
  cell_vector <- gold.data[,cell_name]

  expect_equal(cell_scores$n_features,sum(cell_vector != 0))
  expect_equal(cell_scores$max_count_value,max(cell_vector))
  expect_equal(cell_scores$total_counts,sum(cell_vector))
  expect_equal(cell_scores$unique_count_values,sum(unique(cell_vector)!= 0))
  expect_equal(cell_scores$n_features_above1,sum(cell_vector > 1))
  expect_equal(cell_scores$n_features_above10,sum(cell_vector > 10))
})

test_that("function ScoreGenes, returns correct values", {
  C <- countland(gold.data)
  C <- ScoreGenes(C)

  gene_index <- 13
  gene_scores <- C@gene_scores[gene_index,]
  gene_name <- C@gene_scores$names[gene_index]
  gene_vector <- gold.data[gene_name,]

  expect_equal(gene_scores$n_cells,sum(gene_vector != 0))
  expect_equal(gene_scores$max_count_value,max(gene_vector))
  expect_equal(gene_scores$total_counts,sum(gene_vector))
  expect_equal(gene_scores$unique_count_values,sum(unique(gene_vector)!= 0))
  expect_equal(gene_scores$n_cells_above1,sum(gene_vector > 1))
  expect_equal(gene_scores$n_cells_above10,sum(gene_vector > 10))
})

test_that("function CountIndex, returns correct value", {
  expect_equal(countland:::CountIndex(c(0,0,0,0,5)),1)
  expect_equal(countland:::CountIndex(c(0,1,2,3,4,5)),3)
  expect_equal(countland:::CountIndex(c(0,4,4,4,4,5)),4)
  expect_equal(countland:::CountIndex(c(0,0,1,1,2,2,3,3,4,4,5,5)),4)
})

test_that("functionScoreCells, feature gene string works", {
  vec_test_new <- rep(0,ncol(gold.data))
  vec_test_new[10] <- 1
  new.gold.data <- rbind(gold.data,vec_test_new)
  C <- countland(new.gold.data)
  C <- ScoreCells(C,gene_string="test")
  expect_equal(sum(C@cell_scores$feature_match_counts),1)

  C <- ScoreCells(C,gene_string="siphonophore")
  expect_equal(sum(C@cell_scores$feature_match_counts),0)
})
