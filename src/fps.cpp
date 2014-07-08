//
// fps.cpp
// 
// Created by Vincent Q. Vu on 2014-03-14
// Copyright 2014 Vincent Q. Vu. All rights reserved
// 

#include <RcppArmadillo.h>
#include <algorithm>
#include <vector>
#include <utility>
#include <cmath>

#include "projection.h"
#include "softthreshold.h"
#include "utility.h"

using namespace Rcpp;
using namespace arma;

int fps_admm(const mat& S, const double& ndim, const double& lambda,
             mat& x, mat& z, mat& u, mat& z_old, 
             double& admm_penalty, const double& admm_adjust,
             int maxiter, const double& tolerance)
{
  int niter;
  double rr, ss;

  for(niter = 0; niter < maxiter; niter++) {
    // Store previous value of z
    z_old = z;

    // Fantope projection
    x = z - u + (S / admm_penalty);
    fantope_projection(x, ndim);

    // Elementwise soft thresholding
    z = x + u;
    z.transform( SoftThresholdOp(lambda / admm_penalty) );

    // Dual variable update
    u = u + x - z;

    // Compute primal and dual residual norms
    rr = norm(x - z, "fro");
    ss = admm_penalty * norm(z - z_old, "fro");

    // Check convergence criterion
    if(rr < tolerance && ss < tolerance) {
      niter++;
      break;
    }

    // Penalty adjustment (Boyd, et al. 2010)
    if(rr > 10.0 * ss) {
        admm_penalty = admm_penalty * admm_adjust;
        u = u / admm_adjust;
    } else if(ss > 10.0 * rr) {
        admm_penalty = admm_penalty / admm_adjust;
        u = u * admm_adjust;
    }
  }

  return niter;
}

//' Fantope Projection and Selection
//'
//' This function computes a solution path of the Fantope Projection and
//' Selection (FPS) estimator.  It takes a symmetric matrix \code{S} as input 
//' and returns an object containing a list of projection matrices 
//' estimated by FPS over a sequence of regularization parameter values.
//'
//' By default, a sensible minimum value of the regularization parameter 
//' is automatically chosen so that the minimally regularized solution
//' is different from PCA.
//'
//' @param S              Input matrix (assumed to be symmetric)
//' @param ndim           Target subspace dimension (can be fractional)
//' @param nsol           Number of solutions to compute
//' @param lambda         Vector of regularization parameter values
//' @param lambdamin      Minimum value of lambda (automatic if < 0)
//' @param maxnvar        Suggested maximum number of variables to include 
//'                       (ignored if <= 0)
//' @param maxiter        Maximum number of iterations for each solution
//' @param tolerance      Convergence threshold
//' @param verbose        Level of verbosity (0 = silent, >0 = increasing 
//'                       verbosity)
//'
//' @return Object of class 'fps'
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
//' # Create a synthetic dataset by adding noise variables to the wine data
//' j <- sample(1:ncol(wine), size = 500 - ncol(wine), replace = TRUE)
//' noise <- apply(wine[, j], 2, sample, replace = TRUE)
//' colnames(noise) <- rep('noise', ncol(noise))
//' x <- cbind(wine, noise)
//' out <- fps(cor(x), ndim = 2, verbose = 1)
//'
//' # Choose lambda by cross-validation (this may take a few minutes)
//' \dontrun{
//' cvout <- cv(out, x, FUN = cor, verbose = 1)
//' plot(cvout)
//' v <- coef(out, lambda = cvout$lambda.cv)
//' print(v)
//' }
//'
// [[Rcpp::export]]
List fps(NumericMatrix S, double ndim,
         NumericVector lambda = NumericVector::create(), 
         int nsol = 50, double lambdamin = -1, int maxnvar = -1, 
         int maxiter = 100, double tolerance = 1e-3, int verbose = 0) {

  const mat _S(S.begin(), S.nrow(), S.ncol(), false);

  if(S.nrow() < 2) stop("Expected S to be a matrix");
  if(ndim < 0.0 || ndim > S.nrow()) stop("Expected ndim to be between 0.0 and the number of rows/columns of S");

  // Compute the maximum absolute off-diagonal value for each column of S
  vec maxoffdiag(_S.n_cols), maxoffdiag_sorted;
  compute_maxoffdiag(maxoffdiag, _S);
  maxoffdiag_sorted = sort(maxoffdiag, "descend");

  // Add a constant to the diagonals of S to ensure that
  // S(i,i) >= |S(i,j)| for all i and j
  vec diag_adjusted = _S.diag() - min(_S.diag() - maxoffdiag);

  // Generate lambda sequence if necessary
  vec _lambda;
  if(lambda.size() > 0) {
    _lambda = vec(lambda.begin(), lambda.size(), false);
    nsol = lambda.size();
  } else {
    if(maxnvar > 0 && (uword) maxnvar < _S.n_cols) {
      lambdamin = maxoffdiag_sorted[maxnvar];
    } else if(lambdamin < 0) {
      // lambdamin automatically set to the minimum maximum off-diagonal
      lambdamin = maxoffdiag_sorted[maxoffdiag_sorted.n_elem - 1];
    }
    double lambdamax = std::max(lambdamin, maxoffdiag_sorted[0]);

    _lambda = vec(nsol);
    loglinearseq(_lambda, lambdamin, lambdamax, nsol);
  }

  // Placeholders for solutions
  List            projection(nsol);
  NumericVector   niter(nsol);
  NumericVector   L1(nsol);
  NumericVector   varexplained(nsol);
  NumericMatrix   leverage(_S.n_rows, nsol);
  mat _leverage(leverage.begin(), leverage.rows(), leverage.cols(), false);

  // Set dimnames attribute of leverages array
  List dimnames(S.attr("dimnames"));
  if(dimnames.size() > 0) {
    leverage.attr("dimnames") = List::create(_["variable"] = dimnames[0],
                                             _["lambda"] = R_NilValue);
  }

  // ADMM variables
  mat     x = zeros<mat>(S.nrow(), S.ncol()), 
          z = zeros<mat>(S.nrow(), S.ncol()),
          u = zeros<mat>(S.nrow(), S.ncol()),
          z_old = zeros<mat>(S.nrow(), S.ncol());

  // ADMM parameters
  double tolerance_abs = sqrt(ndim) * tolerance,
         admm_penalty = norm(vectorise(_S), "inf"),
         admm_adjust = 2.0;

  // Outer loop to compute solution path
  for(int i = 0; i < nsol; i++) {
    if(verbose > 0) Rcout << ".";

    // ADMM
    niter[i] = fps_admm(_S, ndim, _lambda[i], 
                        x, z, u, z_old, 
                        admm_penalty, admm_adjust,
                        maxiter, tolerance_abs);

    // Store solution
    NumericMatrix p(_S.n_rows, _S.n_cols);
    p.attr("dimnames") = S.attr("dimnames");
    mat _p(p.begin(), p.nrow(), p.ncol(), false);
    _p = z;
    projection(i) = p;

    L1(i) = norm(vectorise(z), 1);
    varexplained(i) = dot(_S, z);
    _leverage.col(i) = z.diag();

    if(verbose > 1) Rcout << niter[i];
    if(verbose > 2) Rcout << "(" << admm_penalty << ")";
  }

  if(verbose > 0) Rcout << std::endl;

  // Return
  List out = List::create(
    Named("ndim") = ndim,
    Named("lambda") = _lambda,
    Named("projection") = projection,
    Named("leverage") = leverage,
    Named("L1") = L1,
    Named("var.explained") = varexplained,
    Named("var.total") = trace(_S),
    Named("niter") = niter
  );
  out.attr("class") = "fps";

  return out;
}
