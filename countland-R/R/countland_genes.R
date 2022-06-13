color_palette <- c("#8c564b", "#9467bd", "#2ca02c", "#e377c2", "#d62728", "#17becf", "#bcbd22", "#ff7f0e", "#7f7f7f", "#1f77b4")


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
        print(paste0("subsampling ",sum(above_threshold)," genes to a max total counts of ",gene_counts))

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
        print(paste0("subsampling all cells to a standard sequencing depth of ",cell_counts))

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

#' Calculate several scores for count-based gene expression.
#'
#' @param C countland object
#' @param subsample if TRUE, use subsampled counts, otherwise use counts (default=FALSE)
#'
#' @return countland object with slot gene_scores
#' @export
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
#' @export
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
