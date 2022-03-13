Dots <- function(C){
    # logging
    C@dots <- t(C@counts) %*% C@counts

    return(C)
}


ScikitManifoldSpeCalEmbedding <- function(A,n_components,drop_first=TRUE){
	if(drop_first==TRUE){
		n_components <- n_components+1
	}

	# calculate normalized graph laplacian using igraph
	Ai <- as.undirected(graph.adjacency(A,weighted=T))
	L <- laplacian_matrix(Ai,normalized=T)
	L <- L * -1 # flip sign to match scikit

	# calculate diagonal matrix
	w <- A
	diag(w) <- 0
	w <- apply(w,1,sum)
	dd <- sqrt(w)

	# find eigenvectors for graph laplacian with largest magnitude
	diffusion_map <- eigs(L,n_components,sigma=1,which="LM")$vectors
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

Cluster <- function(C,n_clusters){
	if(length(C@embedding == 0) || length(C@dots != 0)){
		# speCal embedding of dot products, if you havent already
		# or if you have already clustered and want to recluster
		if(n_clusters < 3){
			n_components <- 3
		} else {
			n_components <- n_clusters
		}

		C@embedding <- ScikitManifoldSpeCalEmbedding(C@dots,n_components,drop_first=FALSE)
	}

	clust <- kmeans(C@embedding,n_clusters,nstart=10,iter.max=300,algorithm="Lloyd")
	C@cluster_labels <- clust$cluster

	return(C)
}

PlotClusters <- function(C){

	embed <- C@embedding[,2:3]
	embed <- setNames(data.frame(embed),paste("component_",seq_len(2),sep=""))

	ggplot(embed,aes(x = component_1,y = component_2, color=as.character(C@cluster_labels))) +
	geom_point(size=1) +
	guides(color=guide_legend(title="cluster")) +
	scale_color_brewer(palette="Paired")

	# total counts
}

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