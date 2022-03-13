#' @import ggplot2
#'
NULL

#' Calculate pairwise dot products of counts between all cells.
#'
#' @param C countland object
#'
#' @return countland object with slot `dots`
#' @export
Dots <- function(C){
    # logging
    C@dots <- t(C@counts) %*% C@counts

    return(C)
}

#' Recapitulate scikit.manifold.spectral_embedding from python.
#'
#' @param A similarity matrix, dgCMatrix
#' @param n_components number of eigenvectors to retain, integer
#' @param drop_first drop constant eigenvector, boolean
#'
#' @return matrix of eigenvectors
ScikitManifoldSpectralEmbedding <- function(A,n_components,drop_first=TRUE){
	if(drop_first==TRUE){
		n_components <- n_components+1
	}

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
	diffusion_map <- RSpectra::eigs(L,n_components,sigma=1,which="LM")$vectors
	embedding <- diffusion_map / dd

	# flip signs
	max_abs_rows <- apply(abs(embedding),2,which.max)
	signs <- sign(diag(embedding[max_abs_rows,]))
	embedding_sign <- t(t(embedding) * signs)

	if(drop_first==TRUE){
		# drop constant eigenvector
		embed <- embedding_sign[,(n_components-1):n_components]
	} else {
		embed <- embedding_sign
	}
	return(embed)
}

#' Perform spectral clustering on dot products.
#'
#' @param C countland object
#' @param n_clusters number of clusters, integer
#'
#' @return countland object with slot 'cluster_labels'
#' @export
Cluster <- function(C,n_clusters){
	if(length(C@embedding == 0) || length(C@dots != 0)){
		# spectral embedding of dot products, if you havent already
		# or if you have already clustered and want to recluster
		if(n_clusters < 3){
			n_components <- 3
		} else {
			n_components <- n_clusters
		}

		C@embedding <- ScikitManifoldSpectralEmbedding(C@dots,n_components,drop_first=FALSE)
	}

	clust <- kmeans(C@embedding,n_clusters,nstart=10,iter.max=300,algorithm="Lloyd")
	C@cluster_labels <- clust$cluster

	return(C)
}

#' Plot cells using spectral embedding of dot products.
#'
#' @param C countland object
#'
#' @export
PlotClusters <- function(C){

	embed <- C@embedding[,2:3]
	embed <- setNames(data.frame(embed),paste("component_",seq_len(2),sep=""))

	ggplot(embed,aes(x = component_1,y = component_2, color=as.character(C@cluster_labels))) +
	geom_point(size=1) +
	guides(color=guide_legend(title="cluster")) +
	scale_color_brewer(palette="Paired")

	# total counts
}

#' Rank the top marker genes for each cluster from spectral clustering.
#'
#' @param C countland object
#' @param method `prop-zero` to rank by proportion of cells that are non-zero (default), or `rank-sums` to rank using Wilcoxon rank-sums test
#'
#' @return countland object with slots `marker_genes` and `marker_full`
#' @export
RankMarkerGenes <- function(C,method='prop-zero'){
	cluster <- C@cluster_labels
	n <- length(C@names_genes)
	rankdfs <- list()
	for(c in seq_len(length(unique(cluster)))){
		cluster_cells <- which(cluster == c)
		noncluster_cells <- which(cluster != c)
		rankdf <- data.frame()
		if(method=="rank-sums"){
			cluster_counts <- as(t(C@counts[,cluster_cells]),"matrix")
			noncluster_counts <- as(t(C@counts[,noncluster_cells]),"matrix")
			res <- matrixTests::col_wilcoxon_twosample(cluster_counts,noncluster_counts,alternative="greater",exact=F)
			adjpvals <- p.adjust(res$pvalue,method="fdr")
			rankdf <- setNames(data.frame(C@names_genes,seq_len(n),res$statistic,res$pvalue,adjpvals),c("names","gene_index","statistic","pvalue","adj.pvalue"))
			rankdf$"significant" <- (rankdf$"adj.pvalue" < 0.05)
			rankdf <- rankdf[order(adjpvals),]
			rankdf$rank <- seq_len(n)
		} else if(method=="prop-zero"){
			cluster_positive <- apply(C@counts[,cluster_cells],1,function(x){sum(x!=0)})
			noncluster_positive <- apply(C@counts[,noncluster_cells],1,function(x){sum(x!=0)})
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
#'
#' @export
PlotMarker <- function(C,gene_index){

  embed <- C@embedding[,2:3]
  embed <- setNames(data.frame(embed),paste("component_",seq_len(2),sep=""))
  embed$counts <- C@counts[gene_index,]

  ggplot(embed[embed$counts != 0,],aes(x = component_1,y = component_2, color=counts)) +
    geom_point(data = embed[embed$counts == 0,],size=0.5,color="gray") +
    geom_point(size=1) +
    guides(color=guide_legend(title="marker gene counts")) +
    viridis::scale_color_viridis()

  # total counts
}
