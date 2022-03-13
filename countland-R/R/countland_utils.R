#' Recapitulate Seurat normalization
#'
#' @param C countland object
#'
#' @return countland object with slots `norm_factor`, `norm_counts`
Normalize <- function(C){
    # internal
    C@norm_factor <- 10000 / apply(C@counts,2,sum)
    C@norm_counts <- C@counts * C@norm_factor[col(C@counts)]
    return(C)
}

#' Recapitulate Seurat log transformation
#'
#' @param C countland object
#'
#' @return countland object with slots `log_counts`
Log <- function(C){
    # internal
    C@log_counts <- as(log(C@norm_counts + 1),"dgCMatrix")
    return(C)
}

#' Recapitulate Seurat scaling to unit variance
#'
#' @param C countland object
#'
#' @return countland object with slots `scaled_counts`
RescaleVariance <- function(C){
    # internal
    scaled <- scale(t(C@log_counts),center=F)
    scaled[is.na(scaled)] <- 0
    scaled <- as(t(scaled),"dgCMatrix")
    C@scaled_counts <- scaled
    return(C)
}

#' Recapitulate Seurat centering scaled and transformed data
#'
#' @param C countland object
#'
#' @return countland object with slots `centered_counts`
Center <- function(C){
    C@centered_counts <- as(t(scale(t(C@scaled_counts),scale=F)),"dgCMatrix")
    return(C)
}
