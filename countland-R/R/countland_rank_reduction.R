RunIMA <- function(C,features,u_bounds,l_bounds=c(0,0),maxiter=1000000.0,stop_crit=0.0001,subsample=TRUE){
	if(subsample==FALSE){
        sg <- as(t(C@counts),"matrix")
    } else {
        if(length(C@subsample)!=0){
            sg <- as(t(C@subsample),"matrix")
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

PlotIMA <- function(self,x = 1, y = 2,subsample=TRUE){
	if(subsample==FALSE){
        sg <- as(t(C@counts),"matrix")
    } else {
        if(length(C@subsample)!=0){
            sg <- as(t(C@subsample),"matrix")
        } else {
            stop("expecting array of subsampled counts, use subsample() or select subsample=False to use unsampled count matrix")
        }
    }

    loading <- sg %*% (C@matrixV %*% C@matrixLambda)

    ld <- setNames(data.frame(loading[,x],loading[,y],C@cluster_labels),c("f1","f2","cluster"))
    ggplot(ld,aes(x = f1, y = f2, color = as.character(cluster))) +
    geom_point(size=1) +
    xlab(paste("feature: ",x,sep="")) +
    ylab(paste("feature: ",y,sep="")) +
    theme(legend.position = "None")
}

PlotIMAElbow <- function(C,max_features,u_bounds,subsample=TRUE){
	if(subsample==FALSE){
        sg <- as(t(C@counts),"matrix")
    } else {
        if(length(C@subsample)!=0){
            sg <- as(t(C@subsample),"matrix")
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

    ggplot(data.frame("features" = seq(2,max_features),"difference" = norms),aes(x=features,y=difference)) +
    geom_point()

}

SharedCounts <- function(self,n_clusters,n_cells=100,subsample=T){
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

    speCal_embed <- ScikitManifoldSpeCalEmbedding(sharedcounts,n_clusters,drop_first=FALSE)
	speCal_cluster <- kmeans(speCal_embed,n_clusters,nstart=10,iter.max=300,algorithm="Lloyd")

	combX <- sg[filt_gene_index,]
	restX <- sg[-filt_gene_index,]

	C@sum_sharedcounts  <- do.call(rbind,lapply(seq_len(n_clusters),function(x){apply(combX[speCal_cluster$cluster ==x,],2,sum)}))
	C@sum_sharedcounts_all <- rbind(sum_sharedcounts,restX)

	return(C)
}

PlotSharedCounts <- function(self,x = 1, y = 2,subsample=TRUE){
    loading <- C@sum_sharedcounts

    ld <- setNames(data.frame(loading[x,],loading[y,],C@cluster_labels),c("f1","f2","cluster"))
    ggplot(ld,aes(x = f1, y = f2, color = as.character(cluster))) +
    geom_point(size=1) +
    xlab(paste("feature: ",x,sep="")) +
    ylab(paste("feature: ",y,sep="")) +
    theme(legend.position = "None")
}