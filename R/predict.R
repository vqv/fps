# 
# predict.R
# 
# Copyright 2014 Vincent Q. Vu. All rights reserved
# 

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
