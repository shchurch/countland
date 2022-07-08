color_palette <- c("#8c564b", "#9467bd", "#2ca02c", "#e377c2", "#d62728", "#17becf", "#bcbd22", "#ff7f0e", "#7f7f7f", "#1f77b4")

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
