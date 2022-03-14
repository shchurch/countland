#' @importClassesFrom Matrix dgCMatrix
#' @import methods
#' @import stats
NULL

#' An S4 class to represent a countland object
#'
#' @slot counts A dgCMatrix with rows as cells, columns as genes.
#' @slot names_genes A character vector of column names.
#' @slot names_cells A character vector of row names.
#' @slot raw_counts The count dgCMatrix as originally loaded.
#' @slot raw_names_genes The gene name character vector as originally loaded.
#' @slot raw_names_cells The cell name characgter vector as originally loaded.
#' @slot subsample A dgCMatrix with row sums equal.
#' @slot gene_scores A data.frame of gene expression measures.
#' @slot dots A similarity dgCMatrix of dot products.
#' @slot embedding An array of two columns (spectral embeddings).
#' @slot cluster_labels A numeric vector of cluster assignemnts of length n cells.
#' @slot marker_full A list of data.frames with genes ranked for each cluster.
#' @slot marker_genes A data.frame of top ten marker genes per cluster.
#' @slot matrixU A dgCMatrix of dimensions cells x features.
#' @slot matrixV A dgCMatrix of dimensions genes x features.
#' @slot matrixLambda A diagonal dgCMatrix of scaling factors.
#' @slot sharedcounts A similarity dgCMatrix of shared counts between genes.
#' @slot sum_sharedcounts A dgCMatrix with counts summed within gene clusters.
#' @slot sum_sharedcounts_all A dgCMatrix with counts summed and including all genes not present in any cluster.
#' @slot norm_factor A numeric vector of cell normalization factors.
#' @slot norm_counts A dgCMatrix of normalized counts.
#' @slot log_counts A dgCMatrix of log trasnformed counts.
#' @slot scaled_counts A dgCMatrix of counts scaled by gene unit variance.
#' @slot centered_counts A dgCMatrix of counts centered at zero.
#'
#' @export
setClass("countland", slots=list(counts="dgCMatrix",
                                 names_genes="character",
                                 names_cells="character",
                                 raw_counts="dgCMatrix",
                                 raw_names_genes="character",
                                 raw_names_cells="character",
                                 subsample="dgCMatrix",
                                 gene_scores="data.frame",
                                 dots="dgCMatrix",
                                 embedding="array",
                                 cluster_labels="numeric",
                                 marker_full="list",
                                 marker_genes="data.frame",
                                 matrixU="dgCMatrix",
                                 matrixV="dgCMatrix",
                                 matrixLambda="dgCMatrix",
                                 sharedcounts="dgCMatrix",
                                 sum_sharedcounts="dgCMatrix",
                                 sum_sharedcounts_all="dgCMatrix",
                                 norm_factor="numeric",
                                 norm_counts="dgCMatrix",
                                 log_counts="dgCMatrix",
                                 scaled_counts="dgCMatrix",
                                 centered_counts="dgCMatrix")
)

#' Initialize a countland object
#'
#' @param m A matrix of counts (dense or sparse)
#' @param verbose print statements (default=TRUE)
#'
#' @return countland object
#' @export
#'
#' @examples
#' m <- matrix(sample(seq(0,4),prob=c(0.95,0.3,0.1,0.05,0.05),2000,replace=TRUE),ncol=50)
#' rownames(m) <- paste0("gene",seq_len(nrow(m)))
#' colnames(m) <- paste0("cell",seq_len(ncol(m)))
#' C <- countland(m)
countland <- function(m,verbose=TRUE){
    # assertions
    if(class(m)[1] != "dgCMatrix"){
        m <- as(m,"dgCMatrix")
    }

    C <- new("countland",counts=m)
    C@names_genes <- C@counts@Dimnames[[1]]
    C@names_cells <- C@counts@Dimnames[[2]]

    print("countland object")
    print(paste0("Count matrix has ",length(C@names_genes)," genes (rows)"))
    print(paste0("    and ",length(C@names_cells)," cells (columns)"))
    print(paste0("The fraction of entries that are nonzero is ",
                 round(Matrix::nnzero(C@counts)/length(C@counts),4)))

    C@raw_counts <- C@counts
    C@raw_names_genes <- C@names_genes
    C@raw_names_cells <- C@names_cells

    return(C)
}

#' Restore count matrix to original state
#'
#' @param C countland object
#'
#' @return countland object
#' @export
RestoreCounts <- function(C){

    C@counts <- C@raw_counts
    C@names_genes <- C@counts@Dimnames[[1]]
    C@names_cells <- C@counts@Dimnames[[2]]

    return(C)
}

#' Subsets genes using a vector of gene indices
#'
#' @param C  countland object
#' @param gene_indices vector of gene index values
#'
#' @return countland object, count matrix updated
#' @export
SubsetGenes <- function(C,gene_indices){
    C@counts <- C@counts[gene_indices,]
    C@names_genes <- C@names_genes[gene_indices]

    print(paste0("New number of genes: ",length(C@names_genes)))

    return(C)
}

#' Subsets cells using a vecvtor of cell indices
#'
#' @param C countland object
#' @param cell_indices vector of cell index values
#'
#' @return countland object, count matrix updated
#' @export
SubsetCells <- function(C,cell_indices){
    C@counts <- C@counts[,cell_indices]
    C@names_cells <- C@names_cells[cell_indices]

    print(paste0("New number of cells: ",length(C@names_cells)))

    return(C)
}
