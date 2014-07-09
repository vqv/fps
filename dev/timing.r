library(microbenchmark)

set.seed(1)
S <- var(matrix(rnorm(40*100), ncol = 100))
microbenchmark(out <- fps(S, ndim = 10, verbose = 1), times = 10)
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
