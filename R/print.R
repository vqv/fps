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
#' @export
print.fps <- function(x, digits = max(3, getOption("digits") - 3), ...) {
    print(cbind(
      "sqrt(L0)"  = sapply(x$projection, 
                           function(y) ceiling(sqrt(sum(y != 0)))), 
      "%Var"      = signif(x$var.explained / x$var.total * 100, digits), 
      "Lambda"    = signif(x$lambda, digits))
    )
}

#' Print FPS basis coefficients
#'
#' This function prints
#'
#' @param x       fps_coef object
#' @param digits  number of significant digits
#' @param ...     further arguments passed to or from other methods.
#' @export
print.fps_coef <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  print.table(x, digits = digits, zero.print = ".", ...)
}

