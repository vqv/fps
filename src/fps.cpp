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

#include "admm.h"
#include "projection.h"
#include "softthreshold.h"
#include "graphsequence.h"
#include "utility.h"

using namespace Rcpp;
using namespace arma;


// Computes minimum and maximum values for lambda based on theory 
// and heuristic.  The maximum value of lambda is equal to the 
// largest maxoffdiag.  The minimum value of lambda is set 
// equal to an order statistic of maxoffdiag.
void compute_lambdarange(const GraphSeq& gs, 
                         double& lambdamin, double& lambdamax, 
                         const double lambdaminratio, 
                         const int maxnvar, const double ndim) {

  for (const auto& i : gs) {
    if (std::isfinite(i.first)) { 
      lambdamax = i.first; 
      break; 
    }
  }

  if (lambdamin < 0) {
    if (lambdaminratio < 0) {
      // Set lambdamin to the last knot
      lambdamin = std::min(gs.crbegin()->first, lambdamax);
    } else {
      lambdamin = lambdamax * lambdaminratio;
    }
  }

  if (maxnvar <= 0) { 
    return;
  }

  // Find the first knot at which the maximum block size exceeds maxnvar
  auto i = std::lower_bound(gs.cbegin(), gs.cend(), 2 * maxnvar, 
    [](const GraphSeq::value_type& a, const int& b) {
      uword maxsize = 0;
      for (const auto& j : a.second) {
        maxsize = std::max(maxsize, j.second.n_elem);
      }
      return maxsize < (uword) b;
    }
  );
  if (i != gs.cbegin()) { --i; }
  lambdamin = std::min(i->first, lambdamax);
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
//' @param maxnvar        Suggested maximum number of variables to include 
//'                       (ignored if \code{maxnvar <= 0})
//' @param lambdaminratio Minimum value of lambda as a fraction of 
//'                       the automatically determined maximum value of 
//'                       lambda (ignored if \code{lambdaminratio < 0})
//' @param lambdamin      Minimum value of lambda (set automatically if 
//'                       \code{lambdamin < 0})
//' @param lambda         Vector of regularization parameter values; overrides //'                       nsol, maxnvar, and lambdamin if nonempty
//' @param maxiter        Maximum number of iterations for each solution
//' @param tolerance      Convergence threshold
//' @param verbose        Level of verbosity; Silent if \code{verbose = 0}, otherwise 
//'                       display more messages and progress indicators as \code{verbose} 
//'                       increases
//'
//' @return An S3 object of class \code{fps} which is a list with the 
//'         following components:
//'   \item{ndim}{trace (dimension) of the estimate}
//'   \item{lambda}{a vector containing the regularization parameters of each 
//'                 estimate}
//'   \item{projection}{a list containing the the projection matrix estimates}
//'   \item{leverage}{a matrix whose columns are the diagonal entries of the 
//'                   projection matrix estimates}
//'   \item{L1}{a vector of the sum of absolute values of each estimate}
//'   \item{var.explained}{variance explained by each estimate (trace inner 
//'                        product of the projection and \code{S})}
//'   \item{var.total}{total variance (trace of \code{S})}
//'   \item{niter}{a vector containing the number of ADMM iterations for each 
//'                estimate}
//'
//' @details
//' For large input matrices (1000-by-1000 or larger) it is recommended 
//' that the \code{maxnvar} argument be set to a reasonably small number.
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
//' out <- fps(cor(x), ndim = 2, maxnvar = 50, verbose = 1)
//'
//' \dontrun{
//' # Choose lambda by cross-validation (this may take a few minutes)
//' cvout <- cv(out, x, FUN = cor, verbose = 1)
//' plot(cvout)
//' v <- coef(out, lambda = cvout$lambda.cv)
//' print(v)
//' }
//'
// [[Rcpp::export]]
List fps(NumericMatrix S, double ndim, int nsol = 50, 
         int maxnvar = -1, double lambdaminratio = -1, double lambdamin = -1, 
         NumericVector lambda = NumericVector::create(), 
         int maxiter = 100, double tolerance = 1e-3, int verbose = 0) {

  // Sanity checks
  if (S.nrow() < 2) { stop("Expected S to be a matrix"); }
  if (ndim <= 0.0 || ndim >= S.nrow()) { stop("Expected 0 < ndim < nrow(S)"); }
  if (nsol < 1) { stop("Expected nsol > 0");
  if (maxiter < 1) { stop("Expected maxiter > 0");
  if (tolerance <= 0.0) { stop("Expected tolerance > 0");

  // Wrap the input matrix with an arma::mat
  const mat _S(S.begin(), S.nrow(), S.ncol(), false);

  // Compute the sequence of solution graphs
  GraphSeq gs(_S, std::max(0.0, lambdamin));

  // Generate lambda sequence if necessary
  vec _lambda;
  if (lambda.size() > 0) {
    _lambda = arma::sort(vec(lambda.begin(), lambda.size()), "descend");
    nsol = _lambda.n_elem;
    lambdamin = _lambda(nsol - 1);
  } else {
    double lambdamax;
    compute_lambdarange(gs, lambdamin, lambdamax, lambdaminratio, 
                        maxnvar, ndim);
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
  if (dimnames.size() > 0) {
    leverage.attr("dimnames") = List::create(_["variable"] = dimnames[0],
                                             _["lambda"] = R_NilValue);
  }

  // ADMM variables
  mat _z(_S.n_rows, _S.n_cols, fill::zeros),
      _u(_S.n_rows, _S.n_cols, fill::zeros);

  // ADMM parameters
  double tolerance_abs = std::sqrt(ndim) * tolerance,
         admm_penalty = arma::norm(vectorise(_S), "inf"),
         admm_adjust = 2.0;

  // Outer loop to compute solution path
  for (int i = 0; i < nsol; i++) {
    if (verbose > 0) { Rcout << "."; }

    // Find active vertex partition and construct block matrix
    const GraphSeq::partition_t& active = gs.get_active(_lambda[i]);

    SymBlockMap<GraphSeq::vertex_t> block_S(_S, active), 
                                    block_z(_z, active), 
                                    block_u(_u, active);

    // ADMM
    niter[i] = admm(FantopeProjection(ndim), 
                    EntrywiseSoftThreshold(_lambda[i]), 
                    block_S, block_z, block_u, 
                    admm_penalty, admm_adjust,
                    maxiter, tolerance_abs);

    // Restore dense matrices
    block_z.copy_to(_z);
    block_u.copy_to(_u);

    // Store solution
    NumericMatrix p(_S.n_rows, _S.n_cols);
    p.attr("dimnames") = S.attr("dimnames");
    mat _p(p.begin(), p.nrow(), p.ncol(), false);
    block_z.copy_to(_p);
    projection(i) = p;

    L1(i) = sumabs(block_z);
    varexplained(i) = dot(block_S, block_z);
    _leverage.col(i) = _p.diag();

    if (verbose > 1) { Rcout << niter[i]; }
    if (verbose > 2) { Rcout << "(" << admm_penalty << ")"; }
  }

  if (verbose > 0) { Rcout << std::endl; }

  // Return
  List out = List::create(
    Named("ndim") = ndim,
    Named("lambda") = wrap(_lambda.memptr(), 
                           _lambda.memptr() + _lambda.n_elem),
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
