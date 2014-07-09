# 
# predict.R
# 
# Copyright 2014 Vincent Q. Vu. All rights reserved
# 

#' Extract basis vector estimates
#'
#' Returns estimated basis vectors
#'
#' @param       object fps object
#' @param       lambda lambda value to extract
#' @param ...   further arguments passed to or from other methods.
#' @export
coef.fps <- function(object, lambda, ...) {
    i <- which.min(abs(object$lambda - lambda))
    v <- svd(object$projection[[i]], nu = 0, nv = object$ndim)$v
    rownames(v) <- rownames(object$projection[[i]])
    class(v) <- 'fps_coef'
    return(v)
}

projection <- function(object, ...) UseMethod("projection")

#' Extract projection estimates
#'
#' Returns estimated projection matrix
#'
#' @param object  fps object
#' @param lambda  lambda value to extract
#' @export
projection.fps <- function(object, lambda) {
    i <- which.min(abs(object$lambda - lambda))
    return(object$projection[[i]])
}

#' Project data onto estimated subspace
#'
#' Returns the projection of x onto the estimated subspace
#'
#' @param object    fps object
#' @param x         data to project
#' @param lambda    lambda value to extract
#' @param ...       further arguments passed to or from other methods.
#' @export
predict.fps <- function(object, x, lambda, ...) {
    v <- coef.fps(object, lambda)
    return(x %*% v)
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

    out <- list(u = s$u, v = s$v)
    class(out) <- 'svps_coef'
    rownames(out$u) <- rownames(object$projection[[i]])
    rownames(out$v) <- colnames(object$projection[[i]])

    return(out)
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
