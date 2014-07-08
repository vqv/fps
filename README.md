Fantope Projection & Selection
==============================

`fps` is an R package that provides an implementation of an ADMM algorithm for computing the 
Fantope projection and selection estimator.  Most of the package is written in 
C++ using Rcpp and the Armadillo C++ library. The estimator is based on 
a convex relaxation of the sparse PCA problem based on the convex hull 
projection matrices (the Fantope).  This circumvents orthogonality and 
deflation issues that plague previous approaches to sparse PCA.  A 
preliminary report describing this estimator and some of its near-optimal 
statistical properties estimator can found in the 
[NIPS conference proceedings](http://papers.nips.cc/paper/5136-fantope-projection-and-selection-a-near-optimal-convex-relaxation-of-sparse-pca). A longer report with more details and new 
results is forthcoming and will be posted to [arXiv](http://arxiv.org).

Dependencies
------------

1. methods,
2. ggplot2,
3. grid,
4. reshape2,
5. Rcpp (>= 0.11.0),
6. RcppArmadillo (>= 0.4.0),
7. scales

Installation
------------

Use the [devtools](https://github.com/hadley/devtools) R package to install directly from GitHub:

    library(devtools)
    install_github("fps", "vqv")

Example Usage
-------------

    library(fps)
    data(wine)
    out <- fps(cor(wine), ndim = 2)
    plot(out)

Issues
------

1. This package is currently under development and has not been thoroughly tested.
2. The documentation is under development and minimal.
3. The author has only built the package on Mac OS X 10.9.
