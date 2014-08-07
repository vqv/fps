library(microbenchmark)

set.seed(1)
S <- var(matrix(rnorm(40*100), ncol = 100))
microbenchmark(out <- fps(S, ndim = 2, verbose = 1), times = 1)
plot(out)

set.seed(1)
S <- var(matrix(rnorm(400*1000), ncol = 1000))
microbenchmark(out <- fps(S, ndim = 10, maxnvar = 100, verbose = 1), times = 1)
plot(out)

library(elasticnet)
data(pitprops)
pitprops <- as.matrix(pitprops)
microbenchmark(out <- fps(pitprops, ndim = 6, lambdamin = 0, verbose = 1), times = 10)
plot(out)

data(wine)
S <- cor(wine)
microbenchmark(out <- fps(S, ndim = 2, lambdamin = 0, verbose = 1), times = 10)
plot(out)
