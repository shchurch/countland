SubsampleCol <- function(lm,li,j,n_counts){
    new_counts <- table(sample(rep(li,lm),n_counts,replace=F))
    new_x <- as.numeric(new_counts)
    new_i <- as.integer(names(new_counts))
    new_j <- as.integer(rep(j,length(new_x))-1)
    new_dgt <- list(new_x,new_i,new_j)
    names(new_dgt) <- c("x","i","j")
    return(new_dgt)
}

listCols<-function(m){
    #converts a sparse Matrix into a list of its columns
    res<-split(m@x, findInterval(seq_len(nnzero(m)), m@p, left.open=TRUE))
    return(res)
}

Subsample <- function(C,n_counts){
    lms <- listCols(C@counts)
    lis <- split(C@counts@i, findInterval(seq_len(nnzero(C@counts)), C@counts@p, left.open=T))

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

CountIndex<-function(lm){
    tb <- table(lm)
    namestb <- as.numeric(names(tb))
    if(any(tb >= as.numeric(names(tb)))) {
        return(max(namestb[tb >= as.numeric(names(tb))]))
    } else {
        return(0)
    }
}

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

    # fix columns

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

    C@scores_genes <- gene_scores

    return(C)
}