#' Plot FPS x
#'
#' This function plots the regularization path
#'
#' @param x          fps x
#' @param type       type of plot
#' @param ...        further arguments passed to or from other methods.
#' @export
#' @examples
#' data(wine)
#' out <- fps(cor(wine), ndim = 2)
#' plot(out)
plot.fps <- function(x, type = c('leverage', 'variance', 'coherence', 'crossleverage', 'sparsity'), ...) {
    type <- match.arg(type)
    if(type == 'leverage') {
        df <- melt(data.frame(lambda = x$lambda, t(x$leverage)), id.vars = "lambda")
        p <- qplot(x = lambda, y = value, geom = 'line', color = variable, data = df)
        p <- p + xlab(expression(lambda)) + ylab("leverage")
        if(is.null(rownames(x$leverage))) {
            p <- p + theme(legend.position = 'none')
        }
    } else if(type == 'variance') {
        df <- data.frame(lambda = x$lambda, variance = x$var.explained / x$var.total)
        p <- qplot(x = lambda, y = 100 * variance, geom = 'line', data = df, xlab = expression(lambda), ylab = 'percent of variance explained') +
                expand_limits(y=c(0,100))
    } else if(type == 'coherence') {
        df <- data.frame(lambda = x$lambda, coherence = apply(x$leverage, 2, max))
        p <- qplot(x = lambda, y = coherence, geom = 'line', data = df, xlab = expression(lambda), )        
    } else if(type == 'crossleverage') {
        df <- data.frame(lambda = x$lambda, active = sapply(x$x, function(x) sum(x[upper.tri(x)] != 0)))
        p <- qplot(x = lambda, y = active, geom = 'line', data = df, 
                   xlab = expression(lambda), ylab = 'nonzero cross-leverages')
    } else if(type == 'sparsity') {
        df <- data.frame(lambda = x$lambda, active = apply(x$leverage != 0, 2, sum))
        p <- qplot(x = lambda, y = active, geom = 'line', data = df, 
                   xlab = expression(lambda), ylab = 'number of selected variables')
    }
    return(p)
}
