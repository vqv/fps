#' Convert to prcomp object
#'
#' Converts an fps object to a prcomp object
#'
#' @param object fps object
#' @param x data to project
#' @param lambda lambda value to extract
#' @export
as.prcomp <- function(object, x, lambda) {
    rotation <- as.matrix(coef(object, lambda))
    colnames(rotation) <- paste('PC', 1:ncol(rotation), sep = '')
    scores <- as.matrix(x) %*% rotation
    sdev <- as.vector(apply(scores, 2, sd))

    pcobj <- list(
        sdev = sdev,
        sdev.total = sqrt(sum(apply(x, 2, var))),
        rotation = rotation,
        x = as.matrix(scores),
        center = object$center,
        scale = object$scale
    )
    class(pcobj) <- 'prcomp'

    return(pcobj)
}
