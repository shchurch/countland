#' Restore count matrix to original state
#'
#' @param C countland object
#'
#' @return countland object
PrintGeneNumber <- function(C){
  message(paste0("new number of genes: ",nrow(C@counts)))
  message(paste0("new number of cells: ",ncol(C@counts)))
}

#' Subsets genes using a vector of gene indices
#'
#' @param C  countland object
#' @param gene_indices vector of gene index values
#' @param remove_empty filter out cells and genes with no observed counts (default=TRUE)
#'
#' @return countland object, count matrix updated
#' @export
#' @examples
#' gold_path <- system.file("testdata", package = "countland", mustWork = TRUE)
#' gold.data <- Seurat::Read10X(data.dir = gold_path)
#' C <- countland(gold.data)
#' C <- SubsetGenes(C,gene_indices=1:200)
SubsetGenes <- function(C,gene_indices,remove_empty=TRUE){
  C@counts <- C@counts[gene_indices,]
  C@names_genes <- C@names_genes[gene_indices]

  if(remove_empty==TRUE){
    if(C@verbose){message("after subsetting and removing empty cells and genes,")}
    C <- RemoveEmpty(C)
    if(C@verbose){PrintGeneNumber(C)}
  } else {
    if(C@verbose){PrintGeneNumber(C)}
  }

  return(C)
}

#' Subsets cells using a vector of cell indices
#'
#' @param C countland object
#' @param cell_indices vector of cell index values
#' @param remove_empty filter out cells and genes with no observed counts (default=TRUE)
#'
#' @return countland object, count matrix updated
#' @export
#' @examples
#' gold_path <- system.file("testdata", package = "countland", mustWork = TRUE)
#' gold.data <- Seurat::Read10X(data.dir = gold_path)
#' C <- countland(gold.data)
#' C <- SubsetCells(C,cell_indices=1:50)
SubsetCells <- function(C,cell_indices,remove_empty=TRUE){
  C@counts <- C@counts[,cell_indices]
  C@names_cells <- C@names_cells[cell_indices]
  if(remove_empty==TRUE){
    if(C@verbose){message("after subsetting and removing empty cells and genes,")}
    C <- RemoveEmpty(C)
    if(C@verbose){PrintGeneNumber(C)}
  } else {
    if(C@verbose){PrintGeneNumber(C)}
  }

  return(C)
}

#' Restore count matrix to original state
#'
#' @param C countland object
#'
#' @return countland object
#' @export
#' @examples
#' gold_path <- system.file("testdata", package = "countland", mustWork = TRUE)
#' gold.data <- Seurat::Read10X(data.dir = gold_path)
#' C <- countland(gold.data)
#' C <- SubsetGenes(C,gene_indices=1:200)
#' C <- SubsetCells(C,cell_indices=1:50)
#' C <- RestoreCounts(C)
RestoreCounts <- function(C){

  C@counts <- C@raw_counts
  C@names_genes <- C@counts@Dimnames[[1]]
  C@names_cells <- C@counts@Dimnames[[2]]
  if(C@verbose){PrintGeneNumber(C)}

  return(C)
}
