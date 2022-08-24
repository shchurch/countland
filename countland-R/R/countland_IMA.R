color_palette <- c("#8c564b", "#9467bd", "#2ca02c", "#e377c2", "#d62728", "#17becf", "#bcbd22", "#ff7f0e", "#7f7f7f", "#1f77b4")

#' Perform integer matrix approximation on count matrix.
#'
#' @param C countland object
#' @param features target number of features, integer
#' @param u_bounds upper bounds for U and V matrices, vector of length 2
#' @param l_bounds lower bounds for U and V matrices, vector of length 2 (default=c(0,0))
#' @param maxiter maximum number of iterations, integer (default=1000000)
#' @param stop_crit criterion for stopping based on difference between iterations, numeric (default=0.0001)
#' @param subsample if TRUE, use subsampled counts (default), otherwise use counts
#'
#' @return countland object with slots `matrixU`, `matrixV`, `matrixLambda`
#' @export
#' @examples
#' gold_path <- system.file("testdata", package = "countland", mustWork = TRUE)
#' gold.data <- Seurat::Read10X(data.dir = gold_path)
#' C <- countland(gold.data)
#' C <- RunIMA(C,features=10,u_bounds=c(10,10),subsample=FALSE)
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
    params <- IMA_params(features,u_bounds,l_bounds,maxiter,stop_crit)
    res <- IMA(sg,params)

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
#' @return generates plot of cells using integer matrix approximation
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

    ld <- setNames(data.frame(loading[,x],loading[,y],C@cluster_labels,row.names=NULL),c("f1","f2","cluster"))
    ggplot(ld,aes(x = .data$f1, y = .data$f2, color = as.character(.data$cluster))) +
    geom_point(size=1) +
    scale_color_manual(values = colors) +
    xlab(paste("feature: ",x,sep="")) +
    ylab(paste("feature: ",y,sep="")) +
    theme(legend.position = "None")
}

#' Plot the difference between the observed and reconstructed count matrix using integer matrix approximation and a series of total features.
#'
#' @param C countland object
#' @param max_features maximum number of features to assess, integer
#' @param u_bounds upper bounds for U and V matrices, vector of length 2
#' @param subsample if TRUE, use subsampled counts (default), otherwise use counts
#'
#' @return generates elbow plot for the difference between observed and reconstructed matrices as number of features increases
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
    	params = IMA_params(x,u_bounds)
    	res = IMA(sg,params)
    	norm = norm(res[[1]]%*%res[[3]],"F")
    	norm_diff = obs_norm - norm
    	return(norm_diff)
    })

    ggplot(data.frame("features" = seq(2,max_features),"difference" = norms),aes(x=.data$features,y=.data$difference)) +
    geom_point()

}

