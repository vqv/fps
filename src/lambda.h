#ifndef __LAMBDA_H
#define __LAMBDA_H

#include <RcppArmadillo.h>

template <typename G>
arma::vec compute_lambda(const G& gs, double lambdamin, 
                         double lambdaminratio, arma::uword n) {

  // The first knot has infinite weight
  double lambdamax = gs.size() > 1 ? (++gs.cbegin())->first : 0;

  if (lambdamin < 0) {
    if (lambdaminratio < 0) {
      lambdamin = std::min(gs.crbegin()->first, lambdamax);
    } else if (lambdaminratio <= 1.0) {
      lambdamin = lambdamax * lambdaminratio;
    }
  }

  return arma::linspace(lambdamax, lambdamin, n);
}

#endif
