//
// mat.cpp
// 
// Created by Vincent Q. Vu on 2014-08-07
// Copyright 2014 Vincent Q. Vu. All rights reserved
// 
#include "mat.h"
#include <cmath>

// Trace inner product
double dot(const BlockMat& a, const BlockMat& b) {
  double x = 0;
  for (auto ai = a.cbegin(), bi = b.cbegin(); ai != a.cend(); ++ai, ++bi) {
    x += dot(*ai, *bi);
  }
  return x;
}

// dot(a'b, a'b) = trace(a'b b' a)
double dotsquare(const BlockMat& a, const BlockMat& b) {
  double x = 0;
  for (auto ai = a.cbegin(), bi = b.cbegin(); ai != a.cend(); ++ai, ++bi) {
    x += arma::accu(arma::square((*ai).t() * (*bi)));
  }
  return x;
}

// dot(a b', a b') = trace(a'a b' b)
double tdotsquare(const BlockMat& a, const BlockMat& b) {
  double x = 0;
  for (auto ai = a.cbegin(), bi = b.cbegin(); ai != a.cend(); ++ai, ++bi) {
    x += arma::accu(arma::square((*ai) * (*bi).t()));
  }
  return x;
}

// Frobenius distance
double dist(const BlockMat& a, const BlockMat& b) {
  double x = 0;
  for (auto ai = a.cbegin(), bi = b.cbegin(); ai != a.cend(); ++ai, ++bi) {
    x += arma::accu(arma::square(*ai - *bi));
  }
  return std::sqrt(x);
}

double sumabs(const BlockMat& x) {
  double s = 0;
  for(auto& b : x) {
    s += norm(vectorise(b), 1);
  }
  return s;
}

void svd(BlockMat& u, BlockVec& s, BlockMat& v, const BlockMat& x) {
  u.blocks.clear(); s.blocks.clear(); v.blocks.clear();
  for (auto& b : x) {
    u.blocks.push_back(arma::mat());
    s.blocks.push_back(arma::vec());
    v.blocks.push_back(arma::mat());
    arma::svd(u.blocks.back(), s.blocks.back(), v.blocks.back(), b);
  }
}

void eig_sym(BlockVec& eigval, BlockMat& eigvec, const BlockMat& x) {
  eigval.blocks.clear(); eigvec.blocks.clear();
  for (auto& b : x) {
    eigval.blocks.push_back(arma::vec());
    eigvec.blocks.push_back(arma::mat());
    arma::eig_sym(eigval.blocks.back(), eigvec.blocks.back(), b);
  }
}
