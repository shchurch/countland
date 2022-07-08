#' Subsets genes using a vector of gene indices
#'
#' @param C  countland object
#' @param gene_indices vector of gene index values
#' @param remove_empty filter out cells and genes with no observed counts (default=TRUE)
#'
#' @return countland object, count matrix updated
#' @export
SubsetGenes <- function(C,gene_indices,remove_empty=TRUE){
  C@counts <- C@counts[gene_indices,]
  C@names_genes <- C@names_genes[gene_indices]

  if(remove_empty==TRUE){
    print("after subsetting and removing empty cells and genes,")
    C <- RemoveEmpty(C)
    LogGeneNumber(C)
  } else {
    LogGeneNumber(C)
  }

  return(C)
}

#' Subsets cells using a vecvtor of cell indices
#'
#' @param C countland object
#' @param cell_indices vector of cell index values
#' @param remove_empty filter out cells and genes with no observed counts (default=TRUE)
#'
#' @return countland object, count matrix updated
#' @export
SubsetCells <- function(C,cell_indices,remove_empty=TRUE){
  C@counts <- C@counts[,cell_indices]
  C@names_cells <- C@names_cells[cell_indices]
  if(remove_empty==TRUE){
    print("after subsetting and removing empty cells and genes,")
    C <- RemoveEmpty(C)
    LogGeneNumber(C)
  } else {
    LogGeneNumber(C)
  }

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
  LogGeneNumber(C)

  return(C)
}
