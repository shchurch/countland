# This code was copied from the repository at
# github.com/shchurch/integer-matrix-approximation
# on July 11, 2022

#' Parameter class for IMA
#'
#' @param rank target number of features in final matrices
#' @param u_bounds upper bounds on integers
#' @param l_bounds lower bounds on integers (default = c(0,0))
#' @param maxiter maximum number of iterations (default = 1000000)
#' @param stop_crit criterion of difference at which to stop (default = 0.0001)
#'
#' @return parameter object
IMA_params <- function(rank,u_bounds,l_bounds=c(0,0),maxiter=1000000,stop_crit=0.0001){
  max_inner_iter <- 1
  alpha <- 0.01

  params <- list(rank,u_bounds,l_bounds,maxiter,stop_crit,max_inner_iter,alpha)
  names(params) <- c("rank","u_bounds","l_bounds","maxiter","stop_crit","max_inner_iter","alpha")
  return(params)
}

#' Update factor matrix - see SUSTain code
#'
#' @param M matrix to be updated (either U or V)
#' @param coeff matrix used in updating algorithm
#' @param mkrp matrix used in updating algorithm
#' @param mode whether update U or V
#' @param lambda_ scaling matrix
#' @param params parameter object
#'
#' @return updated matrix and scaling factors
IMA_Update_Factor <- function(M,coeff,mkrp,mode,lambda_,params){

  r <- ncol(M)
  for(k in seq(1:r)){
    core <- M %*% (lambda_ * coeff[,k])

    ## will be subtracted to update based on new lambda[k]
    core_k <- M[,k] * (lambda_[k]*coeff[k,k])

    ## update scalar weight of k-th component
    delta_lambda_k <- (t(M[,k]) %*% (mkrp[,k]-core)) / (coeff[k,k]*(norm(M[,k],"2")^2))
    lambda_[k] <- max(1,round(lambda_[k]+delta_lambda_k))

    ## adjust core with new lambda[k]
    core <- core - core_k + (M[,k] * (lambda_[k]*coeff[k,k]))

    ## get k-th column of factor matrix
    Mk <- M[,k] + ((mkrp[,k]-core) / (lambda_[k]*coeff[k,k]))

    ## project to integer constrained set
    Mk <- round(Mk)
    if (!is.infinite(params$l_bounds[mode])){
      Mk[Mk < params$l_bounds[mode]] <- params$l_bounds[mode]
    }
    if (!is.infinite(params$u_bounds[mode])){
      Mk[Mk > params$u_bounds[mode]] <- params$u_bounds[mode]
    }

    ## update k-th column of factor matrix
    M[,k] <- Mk

    ## avoid zero lock
    if(all(M[,k]==0)){
      stopifnot(params$u_bounds[mode] >= 1 && params$l_bounds[mode] >= 0)
      M[sample(seq(1:nrow(M)),1),k] <- 1 # add 1 to random value in k-th column
    }
  }

  return(list(M,lambda_))
}

#' rescale if max val is above upper bound
#'
#' @param h matrix to be rescaled
#' @param l_bound lower bound
#' @param u_bound upper bound
#'
#' @return rescaled matrix
IMA_Compute_Init_Scaled <- function(h,l_bound,u_bound){
  for(j in seq(1:ncol(h))){
    maxval <- max(h[,j])

    # if(maxval==0){
    #    stop("Caution, will divide by zero!")
    #}
    if (maxval>u_bound){
      mult_j <- u_bound / maxval
      h[,j] <- round(mult_j * h[,j])
    }
  }

  h[h < l_bound] <- l_bound
  return(h)
}

#' function to initialize U, V, and Lambda
#'
#' @param X observed data matrix
#' @param params parameter object
#'
#' @return initialized U, V, and Lambda matrices
IMA_init <- function(X,params){
  initU <- matrix(rep(0,nrow(X)*params$rank),ncol=params$rank)

  for(r in seq(1:params$rank)){
    p <- sample(nrow(X),nrow(X))
    nbig <- round((1 / params$rank) * nrow(X))
    initU[p[1:nbig],r] <- sample(seq(params$l_bounds[1]:params$u_bounds[1]),nbig,replace=T)
  }

  if(Matrix::rankMatrix(initU)!=params$rank){
    stop("initializing U failed to achieve target rank, try again")
  }

  rank_flag <- F
  max_rank_flag <- 0
  while(rank_flag == F){
    pat_picked <- sample(seq(1:nrow(X)),params$rank,replace=T)
    initV <- t(X[pat_picked,])
    initV <- IMA_Compute_Init_Scaled(initV,params$l_bounds[2],params$u_bounds[2])

    rank_flag <- (Matrix::rankMatrix(initV)==params$rank)
    max_rank_flag <- max_rank_flag + 1
    if(max_rank_flag == 1000){
      stop("cannot initialize V because rank is not achieved")
    }
  }

  initLambda <- rep(1,params$rank)
  return(list(initU,initV,initLambda))
}

#' run integer matrix approximation
#'
#' @param X observed data matrix
#' @param params parameter object
#'
#' @return U, V, and Lambda matrix factors
IMA <- function(X,params){
  # IMA init

  init_matrices <- IMA_init(X,params)
  U <- init_matrices[[1]]
  V <- init_matrices[[2]]
  lambda_ <- init_matrices[[3]]

  stopifnot(sum(lambda_ >= 1) == length(lambda_))

  nX <- norm(X,"F")^2

  iter_ <- 0

  ## convergence arrays
  e <- list(0)
  fit <- list(0)

  ## main loop
  while(iter_ <= params$maxiter){

    ## Update U and lambda
    A <- X %*% V
    B <- t(V) %*% V
    res1 <- IMA_Update_Factor(U,B,A,1,lambda_,params)
    U <- res1[[1]]
    lambda_ <- res1[[2]]

    ## Update V and lambda
    A <- t(X) %*% U
    B <- t(U) %*% U
    res2 <- IMA_Update_Factor(V,B,A,2,lambda_,params)
    V <- res2[[1]]
    lambda_ <- res2[[2]]

    ## Check convergence
    V_ <- V %*% diag(lambda_)
    e <- append(e,sqrt(max(nX - (2*sum(V_ * A)) + sum(sum( B * (t(V_) %*% V_) )),0)))
    fit <- append(fit,1-(e[[length(e)]]/sqrt(nX)))

    if(iter_ > 0){
      if(abs(fit[[length(fit)-1]] - fit[[length(fit)]]) < params$stop_crit){
        #print("SUSTain_M fit for each iteration")
        #print(fit)
        break
      }
    }

    iter_ <- iter_ + 1
  }

  return(list(U,V,diag(lambda_)))
}
