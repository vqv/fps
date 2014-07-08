#ifndef __PROJECTION_H
#define __PROJECTION_H

#include <RcppArmadillo.h>

void fantope_projection(arma::mat& x, double d);
void singularvalue_projection(arma::mat& x, double d);

#endif
