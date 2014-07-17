Fantope Projection & Selection
==============================

This is a private branch of the fps R package.  The plan is to develop 
new functionality for regularized singular value decompositions here.

Installation
------------

Use [devtools](https://github.com/hadley/devtools) to install directly from GitHub:

```R
library(devtools)
install_bitbucket("svd", "vqv", ref = "svd", 
                  auth_user = "yourbitbucketusername", 
                  password = "yourbitbucketpassword")
```

Example Usage
-------------

```R
library(fps)

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
```

Issues
------

1. This package is currently under development and has not been thoroughly tested.
2. The documentation is under development and minimal.
3. The author has only built the package on Mac OS X 10.9.
