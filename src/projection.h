//
// projection.h
// 
// Copyright 2014 Vincent Q. Vu. All rights reserved
//

#ifndef __PROJECTION_H
#define __PROJECTION_H

#include <RcppArmadillo.h>

void fantope_projection(arma::mat& x, double d);
void singularvalue_projection(arma::mat& x, double d);

#endif
