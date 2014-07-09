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
    v <- svd(object$projection[[i]], nu = 0, nv = ceiling(object$ndim))$v
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
#' @param ...   further arguments passed to or from other methods.
#' @export
projection.fps <- function(object, lambda, ...) {
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
