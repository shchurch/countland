#' @importClassesFrom Matrix dgCMatrix
#' @importFrom ggplot2 ggplot aes geom_jitter position_jitter geom_point guides guide_legend
#' @importFrom ggplot2 scale_color_manual scale_x_continuous theme xlab ylab
#' @importFrom rlang .data
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
#' @slot raw_names_cells The cell name character vector as originally loaded.
#' @slot subsample A dgCMatrix with row sums equal.
#' @slot cell_scores A data.frame of cell count measures.
#' @slot gene_scores A data.frame of gene expression measures.
#' @slot dots A similarity dgCMatrix of dot products.
#' @slot eigenvals An vector of eigenvalues from spectral embedding
#' @slot embedding An array of two columns (spectral embeddings).
#' @slot cluster_labels A numeric vector of cluster assignments of length n cells.
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
#' @slot log_counts A dgCMatrix of log transformed counts.
#' @slot scaled_counts A dgCMatrix of counts scaled by gene unit variance.
#' @slot centered_counts A dgCMatrix of counts centered at zero.
#' @slot verbose A T/F object for suppressing messages
#'
#' @export
setClass("countland", slots=list(counts="dgCMatrix",
                                 names_genes="character",
                                 names_cells="character",
                                 raw_counts="dgCMatrix",
                                 raw_names_genes="character",
                                 raw_names_cells="character",
                                 subsample="dgCMatrix",
                                 cell_scores="data.frame",
                                 gene_scores="data.frame",
                                 dots="dgCMatrix",
                                 eigenvals="numeric",
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
                                 centered_counts="dgCMatrix",
                                 verbose="logical")
)

#' Initialize a countland object from a dgCMatrix
#'
#' @param m A matrix of counts (dense or sparse)
#' @param remove_empty filter out cells and genes with no observed counts (default=TRUE)
#' @param verbose show stderr message statements (default=TRUE)
#'
#' @return countland object
#' @export
#'
#' @examples
#' gold_path <- system.file("testdata", package = "countland", mustWork = TRUE)
#' gold.data <- Seurat::Read10X(data.dir = gold_path)
#' C <- countland(gold.data)
countland <- function(m,remove_empty=TRUE,verbose=TRUE){
    # assertions
    if(class(m)[1] != "dgCMatrix"){
        m <- as(m,"dgCMatrix")
    }

    C <- new("countland",counts=m)
    C@names_genes <- C@counts@Dimnames[[1]]
    C@names_cells <- C@counts@Dimnames[[2]]

    C@verbose=verbose
    if(C@verbose){message("countland object")}
    if(remove_empty==TRUE){
        C <- RemoveEmpty(C)
    }
    if(C@verbose){
        message(paste0("the count matrix has ",nrow(C@counts)," genes (rows)"))
        message(paste0("    and ",ncol(C@counts)," cells (columns)"))
        message(paste0("the fraction of entries that are nonzero is ",
            round(Matrix::nnzero(C@counts)/length(C@counts),4)))
    }

    C@raw_counts <- C@counts
    C@raw_names_genes <- C@names_genes
    C@raw_names_cells <- C@names_cells

    return(C)
}

#' Internal function to remove empty columns and rows
#'
#' @param C  countland object
#'
#' @return countland object, count matrix updated
RemoveEmpty <- function(C){

    cell_indices <- which(diff(C@counts@p)>0)
    gene_indices <- sort(unique(C@counts@i))+1
    C@counts <- C@counts[gene_indices,]
    C@names_genes <- C@names_genes[gene_indices]
    C@counts <- C@counts[,cell_indices]
    C@names_cells <- C@names_cells[cell_indices]
    if(C@verbose){message("after removing empty cells and genes")}

    return(C)
}
