//
// admm.h
// 
// Created by Vincent Q. Vu on 2014-07-08
// Copyright 2014 Vincent Q. Vu. All rights reserved
//

#ifndef __ADMM_H
#define __ADMM_H

#include <RcppArmadillo.h>

template <typename F, typename G> 
int admm(F projection, G selection, 
         const arma::mat& input, const double& lambda,
         arma::mat& x, arma::mat& z, arma::mat& u, arma::mat& z_old, 
         double& admm_penalty, const double& admm_adjust,
         int maxiter, const double& tolerance)
{
  int niter;
  double rr, ss;

  for(niter = 0; niter < maxiter; niter++) {
    // Store previous value of z
    z_old = z;

    // Projection
    x = z - u + (input / admm_penalty);
    projection(x);

    // Selection
    z = x + u;
    selection(z, lambda / admm_penalty);

    // Dual variable update
    u = u + x - z;

    // Compute primal and dual residual norms
    rr = arma::norm(x - z, "fro");
    ss = admm_penalty * arma::norm(z - z_old, "fro");

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

#endif
