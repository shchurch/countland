#' Calculate pairwise dot products of counts between all cells.
#'
#' @param C countland object
#' @param subsample if TRUE, use subsampled counts, otherwise use counts (default=FALSE)
#'
#' @return countland object with slot `dots`
#' @export
#' @examples
#' gold_path <- system.file("testdata", package = "countland", mustWork = TRUE)
#' gold.data <- Seurat::Read10X(data.dir = gold_path)
#' C <- countland(gold.data)
#' C <- Dot(C)
Dot <- function(C,subsample=FALSE){
  if(subsample==FALSE){
    counts <- C@counts
  } else {
    if(length(C@subsample)!=0){
      counts <- C@subsample
    } else {
      stop("expecting array of subsampled counts, use subsample() or select subsample=False to use unsampled count matrix")
    }
  }

  if(C@verbose){message("Calculating dot products between rows...")}
  C@dots <- Matrix::t(counts) %*% counts
  if(C@verbose){message("    done.")}

  return(C)
}
