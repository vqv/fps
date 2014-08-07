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

// Computes a set of indices that contains the active variables of an SVPS
// solution.
void find_active(uvec& active_row, uvec& active_col, 
                 const vec& max_row, const vec& max_col, 
                 const double lambda, const double ndim) {

  active_row = find(max_row > lambda);
  active_col = find(max_col > lambda);

  // Make sure there is at least one active row/col
  if (active_row.n_elem == 0) find(max_row == max(max_row));
  if (active_col.n_elem == 0) find(max_col == max(max_col));
}

// Computes minimum and maximum values for lambda based on theory 
// and heuristic.
void compute_lambdarange(double& lambdamin, double& lambdamax, 
                         const double lambdaminratio, 
                         const vec& max_row, const vec& max_col,
                         const int maxnrow, const int maxncol, 
                         const double ndim) {

  lambdamax = std::max(max_row.max(), max_col.max());

  if (lambdamin < 0) {
    if (lambdaminratio < 0) {
      lambdamin = lambdamax * lambdaminratio;
    } else {
      lambdamin = std::min(max_row.min(), max_col.min());
    }
  }

  if (maxncol > 0 && (uword) maxncol < max_col.n_elem) {
    vec m = sort(max_row, "descend");
    lambdamin = std::max(lambdamin, m[maxncol]);
  }

  if (maxnrow > 0 && (uword) maxnrow < max_row.n_elem) {
    vec m = sort(max_row, "descend");
    lambdamin = std::max(lambdamin, m[maxnrow]);
  }
}

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
//' @param maxnrow        Suggested maximum number of rows to include 
//'                       (ignored if \code{maxnrow <= 0})
//' @param maxncol        Suggested maximum number of cols to include 
//'                       (ignored if \code{maxnrow <= 0})
//' @param lambdaminratio Minimum value of lambda as a fraction of 
//'                       the automatically determined maximum value of 
//'                       lambda (ignored if \code{lambdaminratio < 0})
//' @param lambdamin      Minimum value of lambda (set automatically if 
//'                       \code{lambdamin < 0})
//' @param lambda         Vector of regularization parameter values
//' @param maxiter        Maximum number of iterations for each solution
//' @param tolerance      Convergence threshold
//' @param verbose        Level of verbosity (silent if \code{verbose = 0}; 
//'                       otherwise display more messages and progress 
//'                       indicators as \code{verbose} increases)
//'
//' @details
//'
//' The default automatic choice of lambdamin ensures that the least 
//' regularized estimate omits at least one row or column.
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
//'   \item{var.row}{}
//'   \item{var.col}{}
//'   \item{var.total}{}
//'   \item{niter}{a vector containing the number of ADMM iterations for each 
//'                estimate}
//'
//' @export
//'
//' @examples
//' # Apply FPS to the standardized wine data from the UCI ML repository
//' data(wine)
//' out <- svps(scale(wine), ndim = 3)
//' print(out)
//' plot(out)
//'
// [[Rcpp::export]]
List svps(NumericMatrix x, double ndim,
          int nsol = 50, int maxnrow = -1, int maxncol = -1, 
          double lambdaminratio = -1, double lambdamin = -1, 
          NumericVector lambda = NumericVector::create(), 
          int maxiter = 100, double tolerance = 1e-3, int verbose = 0) {

  // Sanity checks
  if(x.ncol() < 2 || x.nrow() < 2) stop("Expected x to be a matrix");
  if(ndim <= 0.0 || ndim >= std::min(x.nrow(), x.ncol())) {
    stop("Expected ndim to be between 0 and the number of rows/columns of S");
  }
  if(nsol < 1) stop("Expected nsol > 0");
  if(maxiter < 1) stop("Expected maxiter > 0");
  if(tolerance <= 0.0) stop("Expected tolerance > 0");

  // Map x to an arma::mat
  const mat _x(x.begin(), x.nrow(), x.ncol(), false);

  // Compute row- and column-wise infinity norms.
  vec max_row = vectorise(max(abs(_x), 1)), 
      max_col = vectorise(max(abs(_x), 0));

  // Generate lambda sequence if necessary
  vec _lambda;
  if (lambda.size() > 0) {
    _lambda = arma::sort(vec(lambda.begin(), lambda.size()), "descend");
    nsol = _lambda.n_elem;
    lambdamin = _lambda(nsol - 1);
  } else {
    double lambdamax;
    compute_lambdarange(lambdamin, lambdamax, lambdaminratio, 
                        max_row, max_col, maxnrow, maxncol, ndim);
    loglinearseq(_lambda, lambdamin, lambdamax, nsol);
  }

  // Placeholders for solutions
  List            projection(nsol);
  NumericVector   niter(nsol);
  NumericVector   L1(nsol);
  NumericVector   var_row(nsol);
  NumericVector   var_col(nsol);
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

  // Find active variables
  uvec active_row, active_col;
  find_active(active_row, active_col, max_row, max_col, lambdamin, ndim);
  mat x_working = _x.submat(active_row, active_col);
  if (verbose > 0 && (active_row.n_elem < _x.n_rows || 
                      active_col.n_elem < _x.n_cols)) {
    Rcout << "Reduced active set size to " 
          << active_row.n_elem << " by " << active_col.n_elem << std::endl;
  }

  // ADMM variables
  mat z = zeros<mat>(x_working.n_rows, x_working.n_cols),
      u = zeros<mat>(x_working.n_rows, x_working.n_cols);

  // ADMM parameters
  double tolerance_abs = std::sqrt(ndim) * tolerance,
         admm_penalty = norm(vectorise(_x), "inf"),
         admm_adjust = 2.0;

  // Outer loop to compute solution path
  for(int i = 0; i < nsol; i++) {
    if(verbose > 0) Rcout << ".";

    // ADMM
    niter[i] = admm(SingularValueProjection(ndim), 
                    EntrywiseSoftThreshold(_lambda[i]), 
                    x_working, z, u, 
                    admm_penalty, admm_adjust,
                    maxiter, tolerance_abs);

    // Store solution
    NumericMatrix p(_x.n_rows, _x.n_cols);
    p.attr("dimnames") = x.attr("dimnames");
    mat _p(p.begin(), p.nrow(), p.ncol(), false);
    _p.submat(active_row, active_col) = z;
    projection(i) = p;

    L1(i) = norm(vectorise(z), 1);
    var_row(i) = accu(square(x_working.t() * z)); // trace(xx' zz')
    var_col(i) = accu(square(x_working * z.t())); // trace(x'x z'z)
    _leverage_row.col(i) = vectorise(sum(square(_p), 1));
    _leverage_col.col(i) = vectorise(sum(square(_p), 0));

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
    Named("var.row") = var_row,
    Named("var.col") = var_col,
    Named("var.total") = accu(square(_x)),
    Named("niter") = niter
  );
  out.attr("class") = "svps";

  return out;
}
