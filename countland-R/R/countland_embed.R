#' Recapitulate scikit.manifold.spectral_embedding from python.
#'
#' @param A similarity matrix, dgCMatrix
#' @param n_components number of eigenvectors to retain, integer
#'
#' @return matrix of eigenvectors
ScikitManifoldSpectralEmbedding <- function(A,n_components){
	# calculate normalized graph laplacian using igraph
	Ai <- igraph::graph_from_adjacency_matrix(A,mode="undirected",weighted=T,diag=F)
	L <- igraph::laplacian_matrix(Ai,normalized=T)
	L <- L * -1 # flip sign to match scikit

	# calculate diagonal matrix
	w <- A
	diag(w) <- 0
	w <- apply(w,1,sum)
	dd <- sqrt(w)

	# find eigenvectors for graph laplacian with largest magnitude
	eigenv <- RSpectra::eigs(L,n_components,sigma=1,which="LM")
	diffusion_map <- eigenv$vectors
	embedding <- diffusion_map / dd

	# flip signs
	max_abs_rows <- apply(abs(embedding),2,which.max)
	signs <- sign(diag(embedding[max_abs_rows,]))
	embedding_sign <- Matrix::t(Matrix::t(embedding) * signs)


	embed <- embedding_sign

	return(list((-1 * eigenv$values),embed))
}

#' Perform spectral embedding on dot products.
#'
#' @param C countland object
#' @param n_components number of components, integer (default=10)
#'
#' @return countland object with slot `embedding`, `eigenvals`
#' @export
#' @examples
#' gold_path <- system.file("testdata", package = "countland", mustWork = TRUE)
#' gold.data <- Seurat::Read10X(data.dir = gold_path)
#' C <- countland(gold.data)
#' C <- Dot(C)
#' C <- Embed(C,n_components=5)
Embed <- function(C,n_components=10){

  if(C@verbose){message("Performing spectral embedding on dot products...")}

  stopifnot("dot product similarity matrix missing; run Dots() first"= length(C@dots) > 0)

  # set diagonal elements to zero
  A <- C@dots
  diag(A) <- 0

  embed <- ScikitManifoldSpectralEmbedding(A,n_components)
  C@eigenvals <- embed[[1]]
  C@embedding <- embed[[2]]

  if(C@verbose){message("    done.")}

  return(C)
}

#' Plots eigenvalues to investigate the optimal number of clusters
#'
#' @param C countland object
#'
#' @return generates plot of eigenvalues by number of components
#' @export
#' @examples
#' gold_path <- system.file("testdata", package = "countland", mustWork = TRUE)
#' gold.data <- Seurat::Read10X(data.dir = gold_path)
#' C <- countland(gold.data)
#' C <- Dot(C)
#' C <- Embed(C,n_components=5)
#' PlotEigengap(C)
PlotEigengap <- function(C){

  stopifnot("eigenvalues missing; run Embed() first"= length(C@eigenvals) > 0)

  e <- C@eigenvals
  edf <- data.frame(x = seq_len(length(e)), y = e)
  ggplot(edf,aes(x = .data$x, y = .data$y)) + geom_point() +
    xlab("index") +
    ylab("eigenvalues") +
    scale_x_continuous(breaks=seq_len(length(e)))
}

