color_palette <- c("#8c564b", "#9467bd", "#2ca02c", "#e377c2", "#d62728", "#17becf", "#bcbd22", "#ff7f0e", "#7f7f7f", "#1f77b4")

#' Rank the top marker genes for each cluster from spectral clustering.
#'
#' @param C countland object
#' @param method `prop-zero` to rank by proportion of cells that are non-zero (default), or `rank-sums` to rank using Wilcoxon rank-sums test
#' @param subsample if TRUE, use subsampled counts, otherwise use counts (default=FALSE)
#'
#' @return countland object with slots `marker_genes` and `marker_full`
#' @export
#' @examples
#' gold_path <- system.file("testdata", package = "countland", mustWork = TRUE)
#' gold.data <- Seurat::Read10X(data.dir = gold_path)
#' C <- countland(gold.data)
#' C <- Dot(C)
#' C <- Embed(C,n_components=5)
#' C <- Cluster(C,n_clusters=3)
#' C <- RankMarkerGenes(C,method='prop-zero',subsample=FALSE)
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
#' @return generates plot of cells with spectral embedding, colored by marker gene counts
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
