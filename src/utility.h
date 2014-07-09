//
// utility.h
// 
// Copyright 2014 Vincent Q. Vu. All rights reserved
//

#ifndef __UTILITY_H
#define __UTILITY_H

#include <RcppArmadillo.h>

void compute_maxoffdiag(arma::vec& maxoffdiag, const arma::mat& x);
void loglinearseq(arma::vec& out, double min, double max, arma::uword n);

#endif
