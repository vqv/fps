//
// blockmat.h
// 
// Copyright 2014 Vincent Q. Vu. All rights reserved
//

#ifndef __BLOCKMAT_H
#define __BLOCKMAT_H

#include <RcppArmadillo.h>
#include <map>

class BlockMat {
public:
  typedef std::pair<arma::uvec, arma::uvec> index_t;

  typedef std::list<arma::mat>::iterator iterator;
  typedef std::list<arma::mat>::const_iterator const_iterator;
  typedef std::list<arma::mat>::size_type size_type;

  arma::uword n_rows, n_cols;

  size_type size() const { return blocks.size(); }

  const_iterator cbegin() const { return blocks.cbegin(); }
  const_iterator cend() const { return blocks.cend(); }
  const_iterator begin() const { return blocks.begin(); }
  const_iterator end() const { return blocks.end(); }
  iterator begin() { return blocks.begin(); }
  iterator end() { return blocks.end(); }

  void copy_to(arma::mat& X) const {
    auto i = indices.cbegin();
    for (auto b = blocks.cbegin(); b != blocks.cend(); ++b, ++i) {
      X.submat(i->first, i->second) = *b;
    }
  }

  BlockMat() : n_rows(0), n_cols(0), blocks(), indices() {}
  BlockMat(const BlockMat& b) {
    n_rows = b.n_rows;
    n_cols = b.n_cols;
    blocks = b.blocks;
    indices = b.indices;
  }

  double sumabs() const;
  void svd(std::list<arma::mat>& u, std::list<arma::vec>& s, 
           std::list<arma::mat>& v) const;
  void eig_sym(std::list<arma::vec>& eigval, 
               std::list<arma::mat>& eigvec) const;

  BlockMat& operator*=(const double& rhs) {
    for (auto& x : blocks) {
      x *= rhs;
    }
    return *this;
  }

  friend double dot(const BlockMat& a, const BlockMat& b);
  friend double dotsquare(const BlockMat& a, const BlockMat& b);
  friend double tdotsquare(const BlockMat& a, const BlockMat& b);
  friend double dist(const BlockMat& a, const BlockMat& b);

protected:
  std::list<arma::mat> blocks;
  std::list<index_t> indices;

};

template <typename Key>
class BlockMap : public BlockMat {

public:
  typedef typename std::map<Key, BlockMat::index_t> indexmap_t;

  BlockMap(const arma::mat& X, const indexmap_t& indexmap) {
    n_rows = 0;
    n_cols = 0;
    for (const auto& i : indexmap) {
      if(i.second.first.n_elem == 0 || i.second.second.n_elem == 0) {
        continue;
      }
      arma::mat b = X.submat(i.second.first, i.second.second);
      blocks.push_back(std::move(b));
      indices.push_back(i.second);
      n_rows += b.n_rows;
      n_cols += b.n_cols;
    }
  }

};

template <typename Key>
class SymBlockMap : public BlockMat {

public:
  typedef arma::uvec index_t;
  typedef typename std::map<Key, index_t> indexmap_t;

  SymBlockMap(const arma::mat& X, const indexmap_t& indexmap) {
    n_rows = 0;
    n_cols = 0;
    for (const auto& i : indexmap) {
      arma::mat b = X.submat(i.second, i.second);
      blocks.push_back(std::move(b));
      indices.push_back(std::make_pair(i.second, i.second));
      n_rows += b.n_rows;
      n_cols += b.n_cols;
    }
  }

};

#endif
