//
// admm.h
// 
// Created by Vincent Q. Vu on 2014-07-08
// Copyright 2014 Vincent Q. Vu. All rights reserved
//

#ifndef __ADMM_H
#define __ADMM_H

#include <RcppArmadillo.h>
#include "block/mat.h"

/**
 * @brief               Projection and selection ADMM algorithm
 * @details 
 * Implements an ADMM algorithm for solving the optimization problem:
 * \f[
 *   \max_{x \in C} \langle input, x \rangle - R(x)
 * \f]
 * This can be interpreted as a regularized support function where the 
 * regularizer is the function R(x). The working memory for this algorithm 
 * is passed in by reference to the function.
 * 
 * 
 * @param projection    A functor operator()(M&) that implements 
 *                      Euclidean projection onto a convex set
 * @param selection     A functor operator()(M&, const double&) that 
 *                      implements the proximal operator of scaled regularizer 
 * @param input         Input matrix
 * @param z             Solution matrix. Reference to an M of the same  
 *                      dimension as input
 * @param u             Dual variable matrix. Reference to an M of the same 
 *                      dimension as input.
 * @param admm_penalty  Reference to the ADMM penalty parameter; its value will 
 *                      be modified by the function.
 * @param admm_adjust   Factor by which the ADMM penalty can increase/decrease
 * @param maxiter       Maximum number of iterations
 * @param tolerance     Convergence tolerance level for the primal and dual 
 *                      residual norms
 * @return The number of iterations
 */
template <typename F, typename G>
int admm(F projection, G selection, 
         const arma::mat& input, arma::mat& z, arma::mat& u, 
         double& admm_penalty, const double admm_adjust, 
         const int maxiter, const double tolerance)
{
  double rr, ss;
  arma::mat x, z_old;

  int niter = 1;
  do {
    // Store previous value of z
    z_old = z;

    // Projection
    x = z - u + (input / admm_penalty);
    projection(x);

    // Selection
    z = x + u;
    selection(z, 1.0 / admm_penalty);

    // Dual variable update
    u = u + x - z;

    // Compute primal and dual residual norms
    rr = arma::norm(arma::vectorise(x - z));
    ss = admm_penalty * arma::norm(arma::vectorise(z - z_old));

    // Check convergence criterion and return if converged
    if (rr < tolerance && ss < tolerance) { return niter; }

    // Penalty adjustment (Boyd, et al. 2010)
    if (rr > 10.0 * ss) {
        admm_penalty = admm_penalty * admm_adjust;
        u *= (1.0 / admm_adjust);
    } else if (ss > 10.0 * rr) {
        admm_penalty = admm_penalty / admm_adjust;
        u *= admm_adjust;
    }
  } while (niter++ < maxiter);

  // Only if we reached the maximum number of iterations
  return -1;
}

template <typename F, typename G> 
int admm(F projection, G selection, 
         const BlockMat& input, BlockMat& z, BlockMat& u,
         double& admm_penalty, const double admm_adjust, 
         const int maxiter, const double tolerance)
{
  BlockMat x(input.size()), z_old(input.size());

  struct block_data {
    arma::mat const *input;
    arma::mat *z, *u, *x;
  };
  std::vector<block_data> data;

  data.reserve(input.size());

  auto ii = input.cbegin();
  for(auto zi = z.begin(), ui = u.begin(), xi = x.begin();
      ii != input.cend();
      ++ii, ++zi, ++ui, ++xi) {
    block_data d = {&*ii, &*zi, &*ui, &*xi};
    data.push_back(std::move(d));
  }

  double rr, ss;
  int niter = 1;
  do {
    // Store previous value of z
    z_old = z;

    // Projection
    for ( auto& d : data ) { *d.x = *d.z - *d.u + (*d.input / admm_penalty); }
    projection(x);

    // Selection
    for ( auto& d : data ) { *d.z = *d.x + *d.u; }
    selection(z, 1.0 / admm_penalty);

    // Dual variable update
    for ( auto& d : data ) { *d.u += *d.x - *d.z; }

    // Compute primal and dual residual norms
    rr = dist(x, z), 
    ss = admm_penalty * dist(z, z_old);

    // Check convergence criterion and return if converged
    if(rr < tolerance && ss < tolerance) { return niter; }

    // Penalty adjustment (Boyd, et al. 2010)
    if(rr > 10.0 * ss) {
        admm_penalty = admm_penalty * admm_adjust;
        u *= (1.0 / admm_adjust);
    } else if(ss > 10.0 * rr) {
        admm_penalty = admm_penalty / admm_adjust;
        u *= admm_adjust;
    }
  } while(niter++ < maxiter);

  // Only if we reached the maximum number of iterations
  return -1;
}

#endif
