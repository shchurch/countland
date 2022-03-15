#' Recapitulate Seurat normalization
#'
#' @param C countland object
#'
#' @return countland object with slots `norm_factor`, `norm_counts`
Normalize <- function(C){
    # internal
    C@norm_factor <- 10000 / apply(C@counts,2,sum)
    C@norm_counts <- as(C@counts * C@norm_factor[col(C@counts)],"dgCMatrix")
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
    m <- as(Matrix::t(C@log_counts),"matrix")
    scaled <- scale(m,center=F)
    scaled[is.na(scaled)] <- 0
    scaled <- as(Matrix::t(scaled),"dgCMatrix")
    C@scaled_counts <- scaled
    return(C)
}

#' Recapitulate Seurat centering scaled and transformed data
#'
#' @param C countland object
#'
#' @return countland object with slots `centered_counts`
Center <- function(C){
  m <- as(Matrix::t(C@log_counts),"matrix")
  centered <- scale(m,center=T)
  centered[is.na(centered)] <- 0
  centered <- as(Matrix::t(centered),"dgCMatrix")
  C@centered_counts <- centered
  return(C)
}
