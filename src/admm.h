//
// admm.h
// 
// Created by Vincent Q. Vu on 2014-07-08
// Copyright 2014 Vincent Q. Vu. All rights reserved
//

#ifndef __ADMM_H
#define __ADMM_H

template <typename P, typename S, typename D, typename M> 
int admm(P Projection, S Selection, D Distance, 
         const M& input, M& z, M& u,
         double& admm_penalty, const double admm_adjust, 
         const int maxiter, const double tolerance)
{
  M x, z_old;

  double rr, ss;
  int niter = 1;
  do {
    // Store previous value of z
    z_old = z;

    // Projection
    x = z - u;
    Projection(x);

    // Selection
    z = x + u + input / admm_penalty;
    Selection(z, 1.0 / admm_penalty);

    // Dual variable update
    u += x - z;

    // Compute primal and dual residual distances
    rr = Distance(x, z), 
    ss = admm_penalty * Distance(z, z_old);

    // Check convergence criterion and return if converged
    if(rr < tolerance && ss < tolerance) { return niter; }

    // Penalty adjustment (Boyd, et al. 2010)
    if(rr > 10.0 * ss) {
        admm_penalty = admm_penalty * admm_adjust;
        u /= admm_adjust;
    } else if(ss > 10.0 * rr) {
        admm_penalty = admm_penalty / admm_adjust;
        u *= admm_adjust;
    }
  } while(niter++ < maxiter);

  // Only if we reached the maximum number of iterations
  return -1;
}

#endif
