#' Internal function for subsampling a column from a sparse matrix.
#'
#' @param lm column vector
#' @param li row positions
#' @param j column index
#' @param n_counts count to sample
#'
#' @return subsampled column as dgTMatrix components
SubsampleCol <- function(lm,li,j,n_counts){
    new_counts <- table(sample(rep(li,lm),n_counts,replace=F))
    new_x <- as.numeric(new_counts)
    new_i <- as.integer(names(new_counts))
    new_j <- as.integer(rep(j,length(new_x))-1)
    new_dgt <- list(new_x,new_i,new_j)
    names(new_dgt) <- c("x","i","j")
    return(new_dgt)
}

#' Split dgCMatrix into column vectors.
#'
#' @param m dgCMatrix
#'
#' @return list of column vectors, numeric
listCols<-function(m){
    #converts a sparse Matrix into a list of its columns
    res<-split(m@x, findInterval(seq_len(Matrix::nnzero(m)), m@p, left.open=TRUE))
    return(res)
}

#' Subsample cells to a standard number of counts by randomly sampling observations without replacement.
#'
#' @param C countland object
#' @param n_counts number of counts to subsample, must be larger than min total counts per cell
#'
#' @return countland object with slot `subsample`
#' @export
Subsample <- function(C,n_counts){
    lms <- listCols(C@counts)
    lis <- split(C@counts@i, findInterval(seq_len(Matrix::nnzero(C@counts)), C@counts@p, left.open=T))

    v <- lapply(seq_len(C@counts@Dim[[2]]),function(x){SubsampleCol(lms[[x]],lis[[x]],x,n_counts)})
    vd <- do.call(rbind,lapply(v,as.data.frame))
    vt <- as(C@counts,"dgTMatrix")
    vt@i <- vd$i
    vt@j <- vd$j
    vt@x <- vd$x
    vc <- as(vt,"dgCMatrix")

    C@subsample <- vc

    return(C)
}

#' Internal function for calculating count index.
#'
#' @param lm column vector
#'
#' @return count index = largest n where n cells have >= n counts
CountIndex<-function(lm){
    tb <- table(lm)
    namestb <- as.numeric(names(tb))
    if(any(tb >= as.numeric(names(tb)))) {
        return(max(namestb[tb >= as.numeric(names(tb))]))
    } else {
        return(0)
    }
}

#' Calculate several scores for count-based gene expression.
#'
#' @param C countland object
#' @param subsample if TRUE, use subsampled counts (default), otherwise use counts
#'
#' @return countland object with slot gene_scores
#' @export
ScoreGenes <- function(C,subsample=TRUE){

    if(subsample==FALSE){
        sg <- t(C@counts)
    } else {
        if(length(C@subsample)!=0){
            sg <- t(C@subsample)
        } else {
            stop("expecting array of subsampled counts, use subsample() or select subsample=False to use unsampled count matrix")
        }
    }

    fulldf <- setNames(data.frame(C@names_genes,diff(sg@p)),c("names","counts_above0"))

    mx <- vapply(listCols(sg),max,1)
    df <- setNames(data.frame(C@names_genes[as.numeric(names(mx))]),"names")
    df$max_count_value <- mx
    df$total_counts <- vapply(listCols(sg),sum,1)
    df$counts_above1 <- vapply(listCols(sg),function(x){sum(x>1)},1)
    df$counts_above10 <- vapply(listCols(sg),function(x){sum(x>10)},1)
    df$unique_count_values <- vapply(listCols(sg),function(x){length(unique(x))},1)
    df$count_index <- vapply(listCols(sg),CountIndex,1)

    gene_scores <- merge(fulldf,df,by="names",all.x=T)
    gene_scores[is.na(gene_scores)] <- 0

    C@gene_scores <- gene_scores

    return(C)
}
