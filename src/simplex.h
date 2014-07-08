#ifndef __SIMPLEX_H
#define __SIMPLEX_H

#include <RcppArmadillo.h>

int simplex(arma::vec& x, double d, bool interior = false);

#endif
