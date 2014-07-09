//
// svd.cpp
// 
// Created by Vincent Q. Vu on 2014-07-09
// Copyright 2014 Vincent Q. Vu. All rights reserved
// 

#include <RcppArmadillo.h>
#include <algorithm>
#include <vector>
#include <utility>
#include <cmath>

#include "admm.h"
#include "projection.h"
#include "softthreshold.h"
#include "utility.h"

using namespace Rcpp;
using namespace arma;

//' Singular Value Projection and Selection
//'
//' This function computes a solution path of the Singular Value Projection 
//' and Selection (SVPS) estimator.  It takes a data matrix \code{x} as input 
//' and returns an object containing a list of projection matrices 
//' estimated by SVPS over a sequence of regularization parameter values.
//'
//' By default, a sensible minimum value of the regularization parameter 
//' is automatically chosen so that the minimally regularized solution
//' is different from the ordinary SVD.
//'
//' @param x              Input matrix
//' @param ndim           Target subspace dimension (can be fractional)
//' @param nsol           Number of solutions to compute
//' @param lambda         Vector of regularization parameter values
//' @param lambdamin      Minimum value of lambda (set automatically if 
//'                       \code{lambdamin < 0})
//' @param maxnvar        Suggested maximum number of variables to include 
//'                       (ignored if \code{maxnvar <= 0})
//' @param maxiter        Maximum number of iterations for each solution
//' @param tolerance      Convergence threshold
//' @param verbose        Level of verbosity; Silent if \code{verbose = 0}, otherwise 
//'                       display more messages and progress indicators as \code{verbose} 
//'                       increases
//'
//' @return An S3 object of class \code{svps} which is a list with the 
//'         following components:
//'   \item{ndim}{sum of squares (dimension) of the estimate}
//'   \item{lambda}{a vector containing the regularization parameters of each 
//'                 estimate}
//'   \item{projection}{a list containing the the (bi-)projection matrix 
//'                     estimates}
//'   \item{leverage.row}{a matrix whose columns are the row leverages}
//'   \item{leverage.col}{a matrix whose columns are the column leverages}
//'   \item{L1}{a vector of the sum of absolute values of each estimate}
//'   \item{var.explained}{}
//'   \item{var.total}{}
//'   \item{niter}{a vector containing the number of ADMM iterations for each 
//'                estimate}
//'
//' @export
//'
//' @examples
//' # Apply FPS to the standardized wine data from the UCI ML repository
//' data(wine)
//' out <- fps(cor(wine), ndim = 2)
//' print(out)
//' plot(out)
//'
// [[Rcpp::export]]
List svps(NumericMatrix x, double ndim,
          NumericVector lambda = NumericVector::create(), 
          int nsol = 50, double lambdamin = -1, int maxnvar = -1, 
          int maxiter = 100, double tolerance = 1e-3, int verbose = 0) {

  if(x.ncol() < 2 || x.nrow() < 2) {
    stop("Expected x to be a matrix");
  }
  if(ndim < 0.0 || ndim >= std::min(x.nrow(), x.ncol())) {
    stop("Expected ndim to be between 0 and the number of rows/columns of S");
  }

  // Map x to an arma::mat
  const mat _x(x.begin(), x.nrow(), x.ncol(), false);

  // Generate lambda sequence if necessary
  vec _lambda;
  if(lambda.size() > 0) {
    _lambda = vec(lambda.begin(), lambda.size(), false);
    nsol = lambda.size();
  } else {
    vec max_row = sort(vectorise(max(abs(_x), 1)), "descend"),
        max_col = sort(vectorise(max(abs(_x), 0)), "descend");

    if(maxnvar > 0 && (uword) maxnvar < _x.n_rows 
                   && (uword) maxnvar < _x.n_cols) {
      lambdamin = std::max(max_row[maxnvar], 
                           max_col[maxnvar]);
    } else if(lambdamin < 0) {
      lambdamin = std::min(max_row[_x.n_rows - 1], 
                           max_col[_x.n_cols - 1]) * 0.1;
    }
    double lambdamax = std::min(max_row[std::ceil(ndim)], 
                                max_col[std::ceil(ndim)]);

    loglinearseq(_lambda, lambdamin, lambdamax, nsol);
  }

  // Placeholders for solutions
  List            projection(nsol);
  NumericVector   niter(nsol);
  NumericVector   L1(nsol);
  NumericVector   varexplained(nsol);
  NumericMatrix   leverage_row(_x.n_rows, nsol);
  NumericMatrix   leverage_col(_x.n_cols, nsol);
  mat _leverage_row(leverage_row.begin(), _x.n_rows, nsol, false),
      _leverage_col(leverage_col.begin(), _x.n_cols, nsol, false);

  // Set dimnames attribute of leverages array
  List dimnames(x.attr("dimnames"));
  if(dimnames.size() > 0) {
    leverage_row.attr("dimnames") = List::create(_["variable"] = dimnames[0],
                                                 _["lambda"] = R_NilValue);
    leverage_col.attr("dimnames") = List::create(_["variable"] = dimnames[1],
                                                 _["lambda"] = R_NilValue);
  }

  // ADMM variables
  mat     h = zeros<mat>(_x.n_rows, _x.n_cols), 
          z = zeros<mat>(_x.n_rows, _x.n_cols),
          u = zeros<mat>(_x.n_rows, _x.n_cols),
          z_old = zeros<mat>(_x.n_rows, _x.n_cols);

  // ADMM parameters
  double tolerance_abs = std::sqrt(ndim) * tolerance,
         admm_penalty = norm(vectorise(_x), "inf"),
         admm_adjust = 2.0;

  // Outer loop to compute solution path
  for(int i = 0; i < nsol; i++) {
    if(verbose > 0) Rcout << ".";

    // ADMM
    niter[i] = admm(SingularValueProjection(ndim), EntrywiseSoftThreshold(), 
                    _x, _lambda[i], 
                    h, z, u, z_old, 
                    admm_penalty, admm_adjust,
                    maxiter, tolerance_abs);

    // Store solution
    NumericMatrix p(_x.n_rows, _x.n_cols);
    p.attr("dimnames") = x.attr("dimnames");
    mat _p(p.begin(), p.nrow(), p.ncol(), false);
    _p = z;
    projection(i) = p;

    L1(i) = norm(vectorise(z), 1);
    if(_x.n_rows < _x.n_cols) {
      varexplained(i) = accu(square(_x * z.t()));
    } else {
      varexplained(i) = accu(square(_x.t() * z));
    }
    _leverage_row.col(i) = vectorise(sum(square(z), 1));
    _leverage_col.col(i) = vectorise(sum(square(z), 0));

    if(verbose > 1) Rcout << niter[i];
    if(verbose > 2) Rcout << "(" << admm_penalty << ")";
  }

  if(verbose > 0) Rcout << std::endl;

  // Return
  List out = List::create(
    Named("ndim") = ndim,
    Named("lambda") = wrap(_lambda.memptr(), 
                           _lambda.memptr() + _lambda.n_elem),
    Named("projection") = projection,
    Named("leverage.row") = leverage_row,
    Named("leverage.col") = leverage_col,
    Named("L1") = L1,
    Named("var.explained") = varexplained,
    Named("var.total") = accu(square(_x)),
    Named("niter") = niter
  );
  out.attr("class") = "svps";

  return out;
}
