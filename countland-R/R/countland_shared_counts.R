color_palette <- c("#8c564b", "#9467bd", "#2ca02c", "#e377c2", "#d62728", "#17becf", "#bcbd22", "#ff7f0e", "#7f7f7f", "#1f77b4")

#' Combine groups of genes with similar counts by clustering and summing.
#'
#' @param C countland object
#' @param n_clusters number of clusters
#' @param n_cells number of cells to sample for gene clustering
#' @param subsample if TRUE, use subsampled counts (default), otherwise use counts
#'
#' @return countland object with slots `shared_counts`, `sum_sharedcounts`, `sum_sharedcounts_all`
#' @export
#' @examples
#' gold_path <- system.file("testdata", package = "countland", mustWork = TRUE)
#' gold.data <- Seurat::Read10X(data.dir = gold_path)
#' C <- countland(gold.data)
#' C <- SharedCounts(C,n_clusters=10,subsample=FALSE)
SharedCounts <- function(C,n_clusters,n_cells=100,subsample=TRUE){
  if(subsample==FALSE){
    sg <- C@counts
  } else {
    if(length(C@subsample)!=0){
      sg <- C@subsample
    } else {
      stop("expecting array of subsampled counts, use subsample() or select subsample=False to use unsampled count matrix")
    }
  }

  filtX <- sg[,sample(ncol(sg),n_cells,replace=F)]
  filt_gene_index <- which(apply(filtX,1,sum)>0)
  filtX <- filtX[filt_gene_index,]

  manhat <- as(rdist::rdist(filtX,metric="manhattan"),"matrix")
  sumgenes <- apply(filtX,1,sum)
  sumgenes_mat <- outer(sumgenes,sumgenes,FUN="+")
  C@sharedcounts <- as((sumgenes_mat - manhat) / 2,"dgCMatrix")

  spectral_embed <- ScikitManifoldSpectralEmbedding(C@sharedcounts,n_clusters)[[2]]
  spectral_cluster <- kmeans(spectral_embed,n_clusters,nstart=10,iter.max=300,algorithm="Lloyd")

  combX <- sg[filt_gene_index,]
  restX <- sg[-filt_gene_index,]

  C@sum_sharedcounts  <- as(do.call(rbind,lapply(seq_len(n_clusters),function(x){apply(combX[spectral_cluster$cluster ==x,],2,sum)})),"dgCMatrix")
  C@sum_sharedcounts_all <- as(rbind(C@sum_sharedcounts,restX),"dgCMatrix")

  return(C)
}

#' Plot cells using matrix of counts summed by clusters of genes.
#'
#' @param C countland object
#' @param x gene cluster to plot on x-axis, integer (default=1)
#' @param y gene cluster to plot on y-axis, integer (default=2)
#' @param colors color palette for ggplot2, default=palette of 11 colors
#'
#' @return generates plot of cells using shared counts
#' @export
PlotSharedCounts <- function(C,x = 1, y = 2,colors=color_palette){
  loading <- C@sum_sharedcounts

  ld <- setNames(data.frame(loading[x,],loading[y,],C@cluster_labels,row.names=NULL),c("f1","f2","cluster"))
  ggplot(ld,aes(x = .data$f1, y = .data$f2, color = as.character(.data$cluster))) +
    geom_point(size=1) +
    scale_color_manual(values = colors) +
    xlab(paste("feature: ",x,sep="")) +
    ylab(paste("feature: ",y,sep="")) +
    theme(legend.position = "None")
}
