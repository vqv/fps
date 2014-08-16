#ifndef __LAMBDA_H
#define __LAMBDA_H

#include <RcppArmadillo.h>

template <typename G>
arma::vec compute_lambda(const G& gs, double lambdamin, 
                         double lambdaminratio, arma::uword n) {

  double lambdamax = gs.cbegin()->first;

  if (lambdamin < 0) {
    if (lambdaminratio < 0) {
      // Set lambdamin to the last knot with 2 or more blocks
      auto i = gs.crbegin();
      while (i->second.size() == 1) { ++i; }
      if (i == gs.crend()) { 
        lambdamin = lambdamax;
      } else {
        lambdamin = i->first;
      }
    } else if (lambdaminratio <= 1.0) {
      lambdamin = lambdamax * lambdaminratio;
    }
  }

  return arma::linspace(lambdamax, lambdamin, n);
}

#endif
