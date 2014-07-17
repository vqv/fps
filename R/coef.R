#' Extract basis vector estimates
#'
#' Returns estimated basis vectors
#'
#' @param       object fps object
#' @param       lambda lambda value to extract
#' @param ...   further arguments passed to or from other methods
#' @export
coef.fps <- function(object, lambda, ...) {
    i <- which.min(abs(object$lambda - lambda))
    v <- svd(object$projection[[i]], nu = 0, nv = ceiling(object$ndim))$v
    rownames(v) <- rownames(object$projection[[i]])
    class(v) <- c('fps_coef', 'matrix')
    return(v)
}

#' Extract projection matrix estimate
#' 
#' \code{projection} is a generic function for extracting a fitted projection 
#' matrix from a dimension reduction object
#'
#' @param object  object
#' @param ...     other arguments
#' @export
projection <- function(object, ...) UseMethod("projection")

#' Extract projection matrix estimates
#'
#' Returns a fitted projection matrix
#'
#' @param object        fps object
#' @param lambda        lambda value to extract
#' @param fixrank       should the rank of the projection matrix be fixed?
#' @param ...           other arguments
#' @export
#' @examples
#' data(wine)
#' out <- fps(cor(wine), ndim = 2)
#' projection(out, lambda = 0.5)
#' 
projection.fps <- function(object, lambda, fixrank = FALSE, ...) {
    
    if(fixrank) {
      return(tcrossprod(coef(object, lambda)))
    } else {
      i <- which.min(abs(object$lambda - lambda))
      object$projection[[i]] 
    }
}
