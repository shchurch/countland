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
#' @param gene_counts maximum total counts for genes
#' @param cell_counts sequencing depth for all cells, or "min" to use the minimum cell total
#'
#' @return countland object with slot `subsample`
#' @export
#' @examples
#' gold_path <- system.file("testdata", package = "countland", mustWork = TRUE)
#' gold.data <- Seurat::Read10X(data.dir = gold_path)
#' C <- countland(gold.data)
#' C <- Subsample(C,gene_counts=250,cell_counts=100)
Subsample <- function(C,gene_counts=NA,cell_counts=NA){
    stopifnot("must choose either gene_counts or cell_counts"= !(is.na(gene_counts) && is.na(cell_counts)))

    if(!is.na(gene_counts)) {
        # invert the matrix to subsample gene vectors
        t_counts <- Matrix::t(C@counts)

        # get gene vectors from dgC matrix
        gene_vectors <- listCols(t_counts)
        gene_row_positions <- split(t_counts@i, findInterval(seq_len(Matrix::nnzero(t_counts)), t_counts@p, left.open=T))

        # set number to subsample to
        gene_total_counts <- sapply(gene_vectors,sum)
        above_threshold <- gene_total_counts > gene_counts
        gene_total_counts[above_threshold] <- gene_counts
        if(C@verbose){message(paste0("subsampling ",sum(above_threshold)," genes to a max total counts of ",gene_counts))}

        # subsample genes
        new_gene_vectors <- lapply(seq_len(t_counts@Dim[[2]]),
                                   function(x){SubsampleCol(gene_vectors[[x]],
                                                            gene_row_positions[[x]],
                                                            x,
                                                            gene_total_counts[[x]])})

        # build dgT matrix, then convert to dgC
        new_gene_df <- do.call(rbind,lapply(new_gene_vectors,as.data.frame))
        new_gene_mat <- as(t_counts,"dgTMatrix")
        new_gene_mat@i <- new_gene_df$i; new_gene_mat@j <- new_gene_df$j; new_gene_mat@x <- new_gene_df$x
        subsample <- as(Matrix::t(new_gene_mat),"dgCMatrix")
    }

    if(!is.na(cell_counts)) {
        # get the matrix
        if(is.na(gene_counts)){
            counts <- C@counts
        } else {
            counts <- subsample
        }

        # set number to subsample to
        if(cell_counts == "min") {
            cell_counts <- min(apply(counts,2,sum))
        }
        if(C@verbose){message(paste0("subsampling all cells to a standard sequencing depth of ",cell_counts))}

        # subsample cells
        cell_vectors <- listCols(counts)
        cell_row_positions <- split(counts@i, findInterval(seq_len(Matrix::nnzero(counts)), counts@p, left.open=T))
        new_cell_vectors <- lapply(seq_len(counts@Dim[[2]]),
                                   function(x){SubsampleCol(cell_vectors[[x]],
                                                            cell_row_positions[[x]],
                                                            x,
                                                            cell_counts)})

        # build dgT matrix, then convert to dgC
        new_cell_df <- do.call(rbind,lapply(new_cell_vectors,as.data.frame))
        new_cell_mat <- as(counts,"dgTMatrix")
        new_cell_mat@i <- new_cell_df$i; new_cell_mat@j <- new_cell_df$j; new_cell_mat@x <- new_cell_df$x
        subsample <- as(new_cell_mat,"dgCMatrix")
    }

    C@subsample <- subsample
    return(C)
}
