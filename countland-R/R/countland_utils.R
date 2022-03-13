Normalize <- function(C){
    # internal
    C@norm_factor <- 10000 / apply(C@counts,2,sum)
    C@norm_counts <- C@counts * C@norm_factor[col(C@counts)]
    return(C)
}

Log <- function(self){
    # internal
    C@log_counts <- as(log(C@norm_counts + 1),"dgCMatrix")
    return(C)
}

RescaleVariance <- function(self){
    # internal
    scaled <- scale(t(C@log_counts),center=F)
    scaled[is.na(scaled)] <- 0
    scaled <- as(t(scaled),"dgCMatrix")
    C@scaled_counts <- scaled
    return(C)
}

Center <- function(self){
    C@centered_counts <- as(t(scale(t(C@scaled_counts),scale=F)),"dgCMatrix")
    return(C)
}