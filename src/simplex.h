//
// simplex.h
// 
// Copyright 2014 Vincent Q. Vu. All rights reserved
//

#ifndef __SIMPLEX_H
#define __SIMPLEX_H

#include <RcppArmadillo.h>
#include "block/block"

arma::uword simplex(arma::vec& x, double d, bool interior = false);
arma::uvec simplex(block::vec& x, double d, bool interior = false);

#endif
