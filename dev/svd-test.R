n <- 100
m <- 100

u <- matrix(0, nrow = n, ncol = 2)
u[1:5, 1] <- 1
u[11:15, 2] <- 1
u <- qr.Q(qr(u))

v <- matrix(0, nrow = m, ncol = 2)
v[1:5, 1] <- 1
v[11:15, 2] <- 1
v <- qr.Q(qr(v))

mu <- tcrossprod(u,v)
image(mu)

set.seed(1)
x <- 16 * mu + matrix(rnorm(n*m), nrow = n)
image(x)

s <- svd(x)
image(tcrossprod(s$u[,1:2], s$v[,1:2]), main = 'SVD')

out <- svps(x, ndim = 2, verbose = 1)
print(out)
image(projection(out, 1), main = 'SVPS')
