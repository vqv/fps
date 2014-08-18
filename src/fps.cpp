//
// fps.cpp
// 
// Created by Vincent Q. Vu on 2014-03-14
// Copyright 2014 Vincent Q. Vu. All rights reserved
// 

#include <RcppArmadillo.h>
#include <cmath>

#include "admm.h"
#include "block/block.h"
#include "graphseq/graphseq.h"
#include "lambda.h"
#include "projection.h"
#include "softthreshold.h"
#include "distance.h"

using namespace Rcpp;
using namespace arma;

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
//' @param maxblocksize   Suggested maximum block size; ignored if \code{== 0}
//' @param lambdaminratio Minimum value of lambda as a fraction of 
//'                       the automatically determined maximum value of 
//'                       lambda; ignored if \code{< 0}
//' @param lambdamin      Minimum value of lambda; determined automatically if 
//'                       \code{< 0}
//' @param lambda         Vector of regularization parameter values; overrides //'                       nsol, maxblocksize, and lambdamin if nonempty
//' @param maxiter        Maximum number of iterations for each solution
//' @param tolerance      Convergence threshold
//' @param verbose        Level of verbosity; silent if \code{= 0}; otherwise 
//'                       display more messages and progress indicators as 
//'                       \code{verbose} increases
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
//' that the \code{maxblocksize} argument be set to a reasonably small number.
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
//' out <- fps(cor(x), ndim = 2, maxblocksize = 50, verbose = 1)
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
List fps(NumericMatrix S, double ndim, unsigned int nsol = 50, 
         unsigned int maxblocksize = 0, 
         double lambdaminratio = -1, double lambdamin = -1, 
         NumericVector lambda = NumericVector::create(), 
         int maxiter = 100, double tolerance = 1e-3, int verbose = 0) {

  // Sanity checks
  if (S.nrow() < 2) { stop("Expected S to be a matrix"); }
  if (ndim <= 0.0 || ndim >= S.nrow()) { stop("Expected 0 < ndim < nrow(S)"); }
  if (nsol < 1) { stop("Expected nsol > 0"); }
  if (maxiter < 1) { stop("Expected maxiter > 0"); }
  if (tolerance <= 0.0) { stop("Expected tolerance > 0"); }

  // Wrap the input matrix with an arma::mat
  const mat _S(S.begin(), S.nrow(), S.ncol(), false);

  vec _lambda;

  if(lambda.size() > 0) {
    _lambda = arma::vec(lambda.begin(), lambda.size());
    std::sort(_lambda.begin(), _lambda.end(), std::greater<double>());
    lambdamin = _lambda[_lambda.n_elem - 1];
    nsol = _lambda.n_elem;
  }

  // Compute the sequence of solution graphs
  GraphSeq gs(_S, std::fmax(0.0, lambdamin), 
              maxblocksize > 0 ? (uword) maxblocksize : _S.n_cols);

  // Generate lambda sequence if necessary
  if (lambda.size() == 0) {
    _lambda = compute_lambda(gs, lambdamin, lambdaminratio, nsol);
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
         admm_penalty = gs.cbegin()->first,
         admm_adjust = 2.0;

  // Outer loop to compute solution path
  for (unsigned int i = 0; i < nsol; i++) {
    if (verbose > 0) { Rcout << "."; }

#ifdef FPS_DONT_USE_GRAPH_OPTIMIZATION
    // ADMM
    niter[i] = admm(FantopeProjection(ndim), 
                    EntrywiseSoftThreshold(_lambda[i]), 
                    FrobeniusDistance(), 
                    _S, _z, _u, admm_penalty, admm_adjust, maxiter, 
                    tolerance_abs);

    // Store solution
    NumericMatrix p(_S.n_rows, _S.n_cols);
    p.attr("dimnames") = S.attr("dimnames");
    mat _p(p.begin(), p.nrow(), p.ncol(), false);
    _p = _z;
    projection(i) = p;

    L1(i) = norm(vectorise(_z), 1);
    varexplained(i) = dot(_S, _z);
    _leverage.col(i) = _p.diag();
#else
    // Find active vertex partition and construct block matrix
    const GraphSeq::partition_t& active = gs.get_active(_lambda[i]);

    block::symmap<GraphSeq::partition_t> b_S(_S, active), 
                                         b_z(_z, active), 
                                         b_u(_u, active);

    // ADMM
    niter[i] = admm(FantopeProjection(ndim), 
                    EntrywiseSoftThreshold(_lambda[i]), 
                    FrobeniusDistance(), 
                    static_cast<block::mat&>(b_S), 
                    static_cast<block::mat&>(b_z), 
                    static_cast<block::mat&>(b_u), 
                    admm_penalty, admm_adjust, maxiter, tolerance_abs);

    // Restore dense matrices
    b_z.copy_to(_z);
    b_u.copy_to(_u);

    // Store solution
    NumericMatrix p(_S.n_rows, _S.n_cols);
    p.attr("dimnames") = S.attr("dimnames");
    mat _p(p.begin(), p.nrow(), p.ncol(), false);
    b_z.copy_to(_p);
    projection(i) = p;

    L1(i) = block::sumabs(b_z);
    varexplained(i) = block::dot(b_S, b_z);
    _leverage.col(i) = _p.diag();
#endif

    if (verbose > 1) { Rcout << niter[i]; }
    if (verbose > 2) { Rcout << "(" << admm_penalty << ")"; }
  }

  if (verbose > 0) { Rcout << std::endl; }

  // Store ordering of rows/columns
  arma::uvec order(_S.n_cols);
  auto i = order.begin();
  for (const auto& b : gs.crbegin()->second) {
    for (auto bi : b.second) { *i++ = bi + 1; }
  }

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
    Named("order") = wrap(order.memptr(), 
                          order.memptr() + order.n_elem),
    Named("niter") = niter
  );
  out.attr("class") = "fps";

  return out;
}
