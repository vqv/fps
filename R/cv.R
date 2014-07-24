# 
# cv.R
# 
# Copyright 2014 Vincent Q. Vu. All rights reserved
# 

#' Cross-validate
#' 
#' \code{cv} is a generic function for cross-validating
#'
#' @param object  object
#' @param ...     other arguments
#' @export
cv <- function(object, ...) UseMethod("cv")

#' Estimate projection score by cross-validation
#'
#' @param object        fps object
#' @param x             Data
#' @param K             Number of cross-validation folds
#' @param FUN           Function to compute input matrix to fps() from subsets of x
#' @param ...           Additional arguments to fps()
#' @export
cv.fps <- function(object, x, K = 5, FUN = var, ...) {
    cv.score <- matrix(nrow = K, ncol = length(object$lambda))
    cv.total <- numeric(K)
    nobs <- nrow(x)
    xi <- sample(1:K, size = nobs, replace = TRUE)
    ni <- table(xi)

    for(i in 1:K) {
        S.train <- FUN(x[xi != i, ])
        S.test <- FUN(x[xi == i, ])
        out.cv <- fps(S.train, ndim = object$ndim, lambda = object$lambda, ...)
        cv.score[i, ] <- sapply(out.cv$projection, 
                                function(H) { sum(S.test * H) })
        cv.total[i] <- sum(diag(S.test))
    }

    cv.null <- (sum(cv.total * ni) / nobs) * (object$ndim / ncol(x))
    cv.m <- apply(sweep(cv.score, 1, ni, FUN = '*'), 2, sum) / nobs
    cv.se <- sqrt( apply(sweep(cv.score^2, 1, ni, FUN = '*'), 2, sum) / nobs - cv.m^2 )

    out <- list(
        cv = cv.m,
        cv.gap = cv.m - cv.null, 
        cv.se = cv.se, 
        cv.null = cv.null,
        lambda = object$lambda,
        lambda.cv = object$lambda[which.max(cv.m)],
        lambda.1se = max(object$lambda[cv.m >= (max(cv.m) - cv.se[which.max(cv.m)])])
    )
    class(out) <- "fps_cv"
    return(out)
}

#' Plot result of cross-validation
#'
#' @param x     fps_cv object
#' @param ...   further arguments passed to or from other methods.
#' @export
plot.fps_cv <- function(x, ...) {
    df <- data.frame(lambda = x$lambda, 
                     cv = x$cv.gap, cv.se = x$cv.se)

    p <- qplot(x = lambda, y = cv, data = df, geom = 'point', 
               xlab = expression(lambda), ylab = 'cross-validation score')
    p <- p + geom_errorbar(aes(ymin = cv - cv.se, ymax = cv + cv.se), alpha = 1/2)
    p <- p + geom_vline(aes(xintercept = x$lambda.cv,
                            color = 'best lambda'), linetype = 'dashed')
    p <- p + geom_vline(aes(xintercept = x$lambda.1se,
                            color = '1SE lambda'), linetype = 'dashed')
    p <- p + geom_hline(aes(yintercept = 0), color = muted('red'))
    p <- p + scale_x_continuous()

    return(p)
}
