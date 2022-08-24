mat <- rbind(c(116,210,200),
             c(210,386,380),
             c(200,380,401))
dgC <- as(mat,"dgCMatrix")

test_that("function ScikitManifoldSpectralEmbedding, returns correct values", {
  embed <- countland:::ScikitManifoldSpectralEmbedding(dgC,n_components=2)

  ### FROM PYTHON SKLEARN ####

  # from sklearn.manifold import spectral_embedding
  # mat = np.array([[116,210,200],[210,386,380],[200,380,401]])
  # spectral_embedding(mat,n_components=2,drop_first=False)

  ## results
  ## array([[ 0.02515773,  0.04247028],
  ##       [ 0.02515773, -0.01382874],
  ##       [ 0.02515773, -0.01595493]])

  expect_equal(embed[[2]][1,1],0.02515773,tolerance = 0.000001)
  expect_equal(embed[[2]][1,2],0.04247028,tolerance = 0.000001)
  expect_equal(embed[[2]][3,2],-0.01595493,tolerance = 0.000001)
})

test_that("function Embed, returns embedding matrix of correct length", {
  C <- new("countland")
  C@dots <- dgC
  C@verbose=TRUE
  C <- Embed(C,n_components=2)

  expect_length(C@embedding,6)
})

test_that("function PlotEigengap, returns plot object", {
  C <- new("countland")
  C@eigenvals <- c(1,2,3,4,5)
  expect_true(ggplot2::is.ggplot(PlotEigengap(C)))
})

