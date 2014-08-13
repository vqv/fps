#ifndef __SVPS_H
#define __SVPS_H

#include <RcppArmadillo.h>

Rcpp::List svps(Rcpp::NumericMatrix x, double ndim,
                unsigned int nsol = 50, int maxnvar = -1,
                double lambdaminratio = -1, double lambdamin = -1,
                Rcpp::NumericVector lambda = Rcpp::NumericVector::create(),
                int maxiter = 100, double tolerance = 1e-3, int verbose = 0);

#endif