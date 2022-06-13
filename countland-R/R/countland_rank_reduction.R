color_palette <- c("#8c564b", "#9467bd", "#2ca02c", "#e377c2", "#d62728", "#17becf", "#bcbd22", "#ff7f0e", "#7f7f7f", "#1f77b4")

#' Perofrm integer matrix approximatin on count matrix.
#'
#' @param C countland objecvt
#' @param features target number of features, integer
#' @param u_bounds upper bounds for U and V matrices, vector of length 2
#' @param l_bounds lower bounds for U and V matrices, vector of length 2 (default=c(0,0))
#' @param maxiter maximum number of iterations, integer (default=1000000)
#' @param stop_crit criterion for stopping based on difference between iterations, numeric (default=0.0001)
#' @param subsample if TRUE, use subsampled counts (default), otherwise use counts
#'
#' @return countland object with slots `matrixU`, `matrixV`, `matrixLambda`
#' @export
RunIMA <- function(C,features,u_bounds,l_bounds=c(0,0),maxiter=1000000,stop_crit=0.0001,subsample=TRUE){
	if(subsample==FALSE){
        sg <- as(Matrix::t(C@counts),"matrix")
    } else {
        if(length(C@subsample)!=0){
            sg <- as(Matrix::t(C@subsample),"matrix")
        } else {
            stop("expecting array of subsampled counts, use subsample() or select subsample=False to use unsampled count matrix")
        }
    }

    #logging
    params <- IMA::IMA_params(features,u_bounds,l_bounds,maxiter,stop_crit)
    res <- IMA::IMA(sg,params)

    C@matrixU <- as(res[[1]],"dgCMatrix")
    C@matrixV <- as(res[[2]],"dgCMatrix")
    C@matrixLambda <- as(res[[3]],"dgCMatrix")

    return(C)
}

#' Plot cells using integer matrix approximation
#'
#' @param C countland object
#' @param x feature on x-axis, integer (default=1)
#' @param y feature on y-axis, integer (default=2)
#' @param subsample if TRUE, use subsampled counts (default), otherwise use counts
#' @param colors color palette for ggplot2, default=palette of 11 colors
#'
#' @export
PlotIMA <- function(C,x = 1, y = 2, colors=color_palette, subsample=TRUE){
	if(subsample==FALSE){
        sg <- as(Matrix::t(C@counts),"matrix")
    } else {
        if(length(C@subsample)!=0){
            sg <- as(Matrix::t(C@subsample),"matrix")
        } else {
            stop("expecting array of subsampled counts, use subsample() or select subsample=False to use unsampled count matrix")
        }
    }

    loading <- sg %*% (C@matrixV %*% C@matrixLambda)

    ld <- setNames(data.frame(loading[,x],loading[,y],C@cluster_labels),c("f1","f2","cluster"))
    ggplot(ld,aes(x = .data$f1, y = .data$f2, color = as.character(.data$cluster))) +
    geom_point(size=1) +
    scale_color_manual(values = colors) +
    xlab(paste("feature: ",x,sep="")) +
    ylab(paste("feature: ",y,sep="")) +
    theme(legend.position = "None")
}

#' Plot the difference between the observe and reconstructed count matrix using integer matrix approximation and a series of total features.
#'
#' @param C countland object
#' @param max_features maximum number of features to assess, integer
#' @param u_bounds upper bounds for U and V matrices, vector of length 2
#' @param subsample if TRUE, use subsampled counts (default), otherwise use counts
#'
#' @export
PlotIMAElbow <- function(C,max_features,u_bounds,subsample=TRUE){
	if(subsample==FALSE){
        sg <- as(Matrix::t(C@counts),"matrix")
    } else {
        if(length(C@subsample)!=0){
            sg <- as(Matrix::t(C@subsample),"matrix")
        } else {
            stop("expecting array of subsampled counts, use subsample() or select subsample=False to use unsampled count matrix")
        }
    }

    obs_norm <- norm(sg,"F")

    norms <- sapply(seq(2,max_features),function(x){
    	params = IMA::IMA_params(x,u_bounds)
    	res = IMA::IMA(sg,params)
    	norm = norm(res[[1]]%*%res[[3]],"F")
    	norm_diff = obs_norm - norm
    	return(norm_diff)
    })

    ggplot(data.frame("features" = seq(2,max_features),"difference" = norms),aes(x=.data$features,y=.data$difference)) +
    geom_point()

}

#' Combine groups of genes with similar counts by clustering and summing.
#'
#' @param C countland object
#' @param n_clusters nubmer of clusters
#' @param n_cells number of cells to sample for gene clustering
#' @param subsample if TRUE, use subsampled counts (default), otherwise use counts
#'
#' @return countland object with slots `shared_counts`, `sum_sharedcounts`, `sum_sharedcounts_all`
#' @export
SharedCounts <- function(C,n_clusters,n_cells=100,subsample=T){
  if(subsample==FALSE){
      sg <- C@counts
  } else {
      if(length(C@subsample)!=0){
          sg <- C@subsample
      } else {
          stop("expecting array of subsampled counts, use subsample() or select subsample=False to use unsampled count matrix")
      }
  }

  filtX <- sg[,sample(col(sg),n_cells,replace=F)]
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
#' @export
PlotSharedCounts <- function(C,x = 1, y = 2,colors=color_palette){
    loading <- C@sum_sharedcounts

    ld <- setNames(data.frame(loading[x,],loading[y,],C@cluster_labels),c("f1","f2","cluster"))
    ggplot(ld,aes(x = .data$f1, y = .data$f2, color = as.character(.data$cluster))) +
    geom_point(size=1) +
    scale_color_manual(values = colors) +
    xlab(paste("feature: ",x,sep="")) +
    ylab(paste("feature: ",y,sep="")) +
    theme(legend.position = "None")
}
