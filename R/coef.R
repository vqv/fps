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
    v <- Matrix(svd(object$projection[[i]], 
                    nu = 0, nv = ceiling(object$ndim))$v)
    rownames(v) <- rownames(object$projection[[i]])
    return(v)
}

#' Extract basis vector estimates
#'
#' Returns estimated basis vectors
#'
#' @param       object fps object
#' @param       lambda lambda value to extract
#' @param ...   further arguments passed to or from other methods.
#' @export
coef.svps <- function(object, lambda, ...) {
    i <- which.min(abs(object$lambda - lambda))
    s <- svd(object$projection[[i]], 
             nu = ceiling(object$ndim), nv = ceiling(object$ndim))

    out <- list(u = Matrix(s$u), v = Matrix(s$v))
    rownames(out$u) <- rownames(object$projection[[i]])
    rownames(out$v) <- colnames(object$projection[[i]])

    return(out)
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

#' Extract projection estimates
#'
#' Returns estimated projection matrix
#'
#' @param object  fps object
#' @param lambda  lambda value to extract
#' @param type    row-, column-, or bi- projection?
#' @param ...   further arguments passed to or from other methods.
#' @export
projection.svps <- function(object, lambda, 
                            type = c('bi', 'row', 'column'), ...) {
    type <- match.arg(type)
    i <- which.min(abs(object$lambda - lambda))
    if(type == 'row') {
      tcrossprod(object$projection[[i]])
    } else if(type == 'column') {
      crossprod(object$projection[[i]])
    } else {
      object$projection[[i]]
    }
}
