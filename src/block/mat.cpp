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
  u.resize(x.size()); s.resize(x.size()); v.resize(x.size());
  auto ui = u.begin(); auto si = s.begin(); auto vi = v.begin();
  for (auto& b : x) {
    arma::svd(*ui, *si, *vi, b);
    ++ui; ++si; ++vi;
  }
}

void eig_sym(BlockVec& eigval, BlockMat& eigvec, const BlockMat& x) {
  eigval.resize(x.size()); eigvec.resize(x.size());
  auto di = eigval.begin(); auto vi = eigvec.begin();
  for (auto& b : x) {
    arma::eig_sym(*di, *vi, b);
    ++di; ++vi;
  }
}
