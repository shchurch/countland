color_palette <- c("#8c564b", "#9467bd", "#2ca02c", "#e377c2", "#d62728", "#17becf", "#bcbd22", "#ff7f0e", "#7f7f7f", "#1f77b4")

#' Calculate pairwise dot products of counts between all cells.
#'
#' @param C countland object
#' @param subsample if TRUE, use subsampled counts, otherwise use counts (default=FALSE)
#'
#' @return countland object with slot `dots`
#' @export
Dot <- function(C,subsample=FALSE){
    if(subsample==FALSE){
      counts <- C@counts
    } else {
      if(length(C@subsample)!=0){
        counts <- C@subsample
      } else {
        stop("expecting array of subsampled counts, use subsample() or select subsample=False to use unsampled count matrix")
      }
    }

    print("Calculating dot products between rows...")
    C@dots <- Matrix::t(counts) %*% counts
    print("    done.")

    return(C)
}

#' Recapitulate scikit.manifold.spectral_embedding from python.
#'
#' @param A similarity matrix, dgCMatrix
#' @param n_components number of eigenvectors to retain, integer
#'
#' @return matrix of eigenvectors
ScikitManifoldSpectralEmbedding <- function(A,n_components){
	# calculate normalized graph laplacian using igraph
	Ai <- igraph::as.undirected(igraph::graph.adjacency(A,weighted=T))
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
Embed <- function(C,n_components=10){

  print("Performing spectral embedding on dot products...")

  stopifnot("dot product similarity matrix missing; run Dots() first"= length(C@dots) > 0)

  # set diagonal elements to zero
  A <- C@dots
  diag(A) <- 0

  embed <- ScikitManifoldSpectralEmbedding(A,n_components)
  C@eigenvals <- embed[[1]]
  C@embedding <- embed[[2]]

  print("    done.")

  return(C)
}

#' Plots eigenvalues to investigate the optimal number of clusters
#'
#' @param C countland object
#'
#' @export
PlotEigengap <- function(C){

  stopifnot("eigenvalues missing; run Embed() first"= length(C@eigenvals) > 0)

  e <- C@eigenvals
  edf <- data.frame(x = seq_len(length(e)), y = e)
  ggplot(edf,aes(x = .data$x, y = .data$y)) + geom_point() +
    xlab("index") +
    ylab("eigenvalues") +
    scale_x_continuous(breaks=seq_len(length(e)))

  #optimal_k <- which.max(diff(C@eigenvals))
  #if(optimal_k < min_clusters){optimal_k <- min_clusters}
  #return(optimal_k)
}

#optimal_k = np.argmax(np.diff(C.eigevals)) + 1
#if(optimal_k) < min_clusters:
#    optimal_k = min_clusters

#return(optimal_k)

#' Perform spectral clustering on dot products.
#'
#' @param C countland object
#' @param n_clusters number of clusters, integer
#' @param n_components number of compnonets from spectral embedding to use (default NULL, will be set to n_clusters), integer
#'
#' @return countland object with slot `cluster_labels`
#' @export
Cluster <- function(C,n_clusters,n_components=NULL){

  print("Performing spectral clustering on dot products...")

  stopifnot("embedding missing; run Embed() first"= length(C@embedding) > 0)

  E <- C@embedding
  if(is.null(n_components)){
    comp <- n_clusters
  } else {
    comp <- n_components
  }

  if(ncol(E) < comp){
    A <- C@dots
    diag(A) <- 0
    E <- ScikitManifoldSpectralEmbedding(A,comp)[[2]]
  } else {
    E <- E[,1:comp]
  }

  clust <- kmeans(E,n_clusters,nstart=10,iter.max=300,algorithm="Lloyd")
	C@cluster_labels <- clust$cluster
  print("    done.")

	return(C)
}

#' Plot cells using spectral embedding of dot products.
#'
#' @param C countland object
#' @param colors color palette for ggplot2, default=palette of 11 colors
#'
#' @export
PlotEmbedding <- function(C,colors=color_palette){

  stopifnot("embedding missing; run Embed() first"= length(C@embedding) > 0)

	embed <- C@embedding[,2:3]
	embed <- setNames(data.frame(embed),paste("component_",seq_len(2),sep=""))

	ggplot2::ggplot(embed,aes(x = .data$component_1,y = .data$component_2, color=as.character(C@cluster_labels))) +
	geom_point(size=1) +
	guides(color=guide_legend(title="cluster")) +
	scale_color_manual(values=colors)
}

#' Rank the top marker genes for each cluster from spectral clustering.
#'
#' @param C countland object
#' @param method `prop-zero` to rank by proportion of cells that are non-zero (default), or `rank-sums` to rank using Wilcoxon rank-sums test
#' @param subsample if TRUE, use subsampled counts, otherwise use counts (default=FALSE)
#'
#' @return countland object with slots `marker_genes` and `marker_full`
#' @export
RankMarkerGenes <- function(C,method='prop-zero',subsample=FALSE){
  if(subsample==FALSE){
    counts <- C@counts
  } else {
    if(length(C@subsample)!=0){
      counts <- C@subsample
    } else {
      stop("expecting array of subsampled counts, use subsample() or select subsample=False to use unsampled count matrix")
    }
  }

  cluster <- C@cluster_labels
	n <- length(C@names_genes)
	rankdfs <- list()
	for(c in seq_len(length(unique(cluster)))){
		cluster_cells <- which(cluster == c)
		noncluster_cells <- which(cluster != c)
		rankdf <- data.frame()
		if(method=="rank-sums"){
			cluster_counts <- as(Matrix::t(counts[,cluster_cells]),"matrix")
			noncluster_counts <- as(Matrix::t(counts[,noncluster_cells]),"matrix")
			res <- matrixTests::col_wilcoxon_twosample(cluster_counts,noncluster_counts,alternative="greater",exact=F)
			adjpvals <- p.adjust(res$pvalue,method="fdr")
			rankdf <- setNames(data.frame(C@names_genes,seq_len(n),res$statistic,res$pvalue,adjpvals),c("names","gene_index","statistic","pvalue","adj.pvalue"))
			rankdf$"significant" <- (rankdf$"adj.pvalue" < 0.05)
			rankdf <- rankdf[order(adjpvals),]
			rankdf$rank <- seq_len(n)
		} else if(method=="prop-zero"){
			cluster_positive <- apply(counts[,cluster_cells],1,function(x){sum(x!=0)})
			noncluster_positive <- apply(counts[,noncluster_cells],1,function(x){sum(x!=0)})
			prop_cluster <- cluster_positive / length(cluster_cells)
			prop_noncluster <- noncluster_positive / length(noncluster_cells)
			prop_a <- prop_cluster - prop_noncluster
			rankdf <- setNames(data.frame(C@names_genes,seq_len(n),prop_a),c("names","gene_index","diff.zeros"))
			rankdf <- rankdf[order(-prop_a),]
			rankdf$rank <- seq_len(n)
		}
		rankdf$cluster_label <- c
		rankdfs[[c]] <- rankdf
	}

	C@marker_genes <- do.call(rbind,lapply(rankdfs,function(x){x[1:10,]}))
	C@marker_full <- rankdfs

	return(C)
}

#' Plot cell using spectral embedding and display counts in a given gene.
#'
#' @param C countland object
#' @param gene_index index value for gene to visualize
#' @param colors color palette for ggplot2, default=palette of 11 colors
#'
#' @export
PlotMarker <- function(C,gene_index,colors=color_palette){

  embed <- C@embedding[,2:3]
  embed <- setNames(data.frame(embed),paste("component_",seq_len(2),sep=""))
  embed$counts <- C@counts[gene_index,]

  g1 <- ggplot2::ggplot(embed,aes(x = .data$component_1,y = .data$component_2, color=as.character(C@cluster_labels))) +
    geom_point(size=1) +
    guides(color=guide_legend(title="cluster")) +
    scale_color_manual(values=colors)
  g2 <- ggplot2::ggplot(embed[embed$counts != 0,],aes(x = .data$component_1,y = .data$component_2, color=.data$counts)) +
    geom_point(data = embed[embed$counts == 0,],size=0.5,color="gray") +
    geom_point(size=1) +
    guides(color=guide_legend(title="marker gene counts")) +
    viridis::scale_color_viridis()
  gridExtra::grid.arrange(g1,g2,ncol=2)
}
