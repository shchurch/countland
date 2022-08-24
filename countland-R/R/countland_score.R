color_palette <- c("#8c564b", "#9467bd", "#2ca02c", "#e377c2", "#d62728", "#17becf", "#bcbd22", "#ff7f0e", "#7f7f7f", "#1f77b4")

#' Internal function for calculating count index.
#'
#' @param lm column vector
#'
#' @return count index = largest n where n cells have >= n counts
CountIndex<-function(lm){
  values <- seq(min(lm),max(lm))
  vals_geq <- sapply(values,function(x){sum(lm >= x)})

  if(any(vals_geq >= values)) {
    return(max(values[vals_geq >= values]))
  } else {
    return(0)
  }
}

#' Calculate several scores for counts across cells
#'
#' @param C countland object
#' @param gene_string string with regular expression expression matching gene names of interest (default=NULL)
#'
#' @return countland object with slot cell_scores
#' @export
#' @examples
#' gold_path <- system.file("testdata", package = "countland", mustWork = TRUE)
#' gold.data <- Seurat::Read10X(data.dir = gold_path)
#' C <- countland(gold.data)
#' C <- ScoreCells(C,gene_string="*149932$")
ScoreCells <- function(C,gene_string=NULL){
  cts <- C@counts

  fulldf <- setNames(data.frame(C@names_cells,diff(cts@p)),c("names","n_features"))

  cts_cols <- listCols(cts)
  mx <- vapply(cts_cols,max,1)
  df <- setNames(data.frame(C@names_cells[as.numeric(names(mx))]),"names")
  df$max_count_value <- mx
  df$total_counts <- vapply(cts_cols,sum,1)
  df$n_features_above1 <- vapply(cts_cols,function(x){sum(x>1)},1)
  df$n_features_above10 <- vapply(cts_cols,function(x){sum(x>10)},1)
  df$unique_count_values <- vapply(cts_cols,function(x){length(unique(x))},1)
  df$count_index <- vapply(cts_cols,CountIndex,1)

  if(!is.null(gene_string)){
    gene_string_match <- grep(gene_string,C@names_genes)
    if(length(gene_string_match) == 1){
      df$feature_match_counts <- C@counts[gene_string_match,]
    } else {
      df$feature_match_counts <- apply(C@counts[gene_string_match,],2,sum)
    }
  }

  cell_scores <- merge(fulldf,df,by="names",all.x=T)
  cell_scores[is.na(cell_scores)] <- 0

  C@cell_scores <- cell_scores

  return(C)
}

#' Calculate several scores for count-based gene expression.
#'
#' @param C countland object
#' @param subsample if TRUE, use subsampled counts, otherwise use counts (default=FALSE)
#'
#' @return countland object with slot gene_scores
#' @export
#' @examples
#' gold_path <- system.file("testdata", package = "countland", mustWork = TRUE)
#' gold.data <- Seurat::Read10X(data.dir = gold_path)
#' C <- countland(gold.data)
#' C <- ScoreGenes(C)
ScoreGenes <- function(C,subsample=FALSE){

  if(subsample==FALSE){
    sg <- Matrix::t(C@counts)
  } else {
    if(length(C@subsample)!=0){
      sg <- Matrix::t(C@subsample)
    } else {
      stop("expecting array of subsampled counts, use subsample() or select subsample=False to use unsampled count matrix")
    }
  }

  fulldf <- setNames(data.frame(C@names_genes,diff(sg@p)),c("names","n_cells"))

  mx <- vapply(listCols(sg),max,1)
  df <- setNames(data.frame(C@names_genes[as.numeric(names(mx))]),"names")
  df$max_count_value <- mx
  df$total_counts <- vapply(listCols(sg),sum,1)
  df$n_cells_above1 <- vapply(listCols(sg),function(x){sum(x>1)},1)
  df$n_cells_above10 <- vapply(listCols(sg),function(x){sum(x>10)},1)
  df$unique_count_values <- vapply(listCols(sg),function(x){length(unique(x))},1)
  df$count_index <- vapply(listCols(sg),CountIndex,1)

  gene_scores <- merge(fulldf,df,by="names",all.x=T)
  gene_scores[is.na(gene_scores)] <- 0

  C@gene_scores <- gene_scores

  return(C)
}

#' Generate a strip plot for counts across selected genes
#'
#' @param C countland object
#' @param gene_indices vector of gene index values
#' @param colors color palette for ggplot2, default=palette of 11 colors
#'
#' @return generates plot of gene count distributions
#' @export
#' @examples
#' gold_path <- system.file("testdata", package = "countland", mustWork = TRUE)
#' gold.data <- Seurat::Read10X(data.dir = gold_path)
#' C <- countland(gold.data)
#' PlotGeneCounts(C,gene_indices=1:10)
PlotGeneCounts <- function(C,gene_indices,colors=color_palette){
  counts <- t(as(C@counts[gene_indices,],"matrix"))
  new_counts <- do.call(rbind,lapply(seq_len(ncol(counts)),function(x){data.frame("name" = colnames(counts)[x], "counts" = counts[,x])}))
  if(length(gene_indices) < 10){
    pal <- colors[1:length(gene_indices)]
  } else {
    pal <- rep("black",length(gene_indices))
  }
  ggplot(new_counts,aes(y=.data$name,x=as.numeric(.data$counts),color=.data$name)) +
    geom_jitter(size=0.5,position = position_jitter(width=0)) +
    scale_color_manual(values=pal) +
    xlab("counts") +
    ylab("gene name") +
    theme(legend.position = "none")
}


