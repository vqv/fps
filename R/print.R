# 
# print.R
# 
# Copyright 2014 Vincent Q. Vu. All rights reserved
# 

#' Print FPS object
#'
#' This function prints
#'
#' @param x       fps object
#' @param digits  number of significant digits
#' @param ...     further arguments passed to or from other methods.
#' @method print fps
#' @export
print.fps <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  print(cbind(
    "sqrt(L0)"  = sapply(x$projection, 
                         function(y) ceiling(sqrt(sum(y != 0)))), 
    "%Var"      = signif(x$var.explained / x$var.total * 100, digits), 
    "Lambda"    = signif(x$lambda, digits))
  )
  invisible(x)
}

#' Print SVPS object
#'
#' This function prints
#'
#' @param x       svps object
#' @param digits  number of significant digits
#' @param ...     further arguments passed to or from other methods.
#' @method print svps
#' @export
print.svps <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  print(cbind(
    "sqrt(L0)"  = sapply(x$projection, 
                         function(y) ceiling(sqrt(sum(y != 0)))), 
    "%Var(row)" = signif(x$var.row / x$var.total * 100, digits), 
    "%Var(col)" = signif(x$var.col / x$var.total * 100, digits), 
    "Lambda"    = signif(x$lambda, digits))
  )
  invisible(x)
}
