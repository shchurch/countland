#' Calculate pairwise dot products of counts between all cells.
#'
#' @param C countland object
#' @param subsample if TRUE, use subsampled counts, otherwise use counts (default=FALSE)
#'
#' @return countland object with slot `dots`
#' @export
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

  print("Calculating dot products between rows...")
  C@dots <- Matrix::t(counts) %*% counts
  print("    done.")

  return(C)
}
