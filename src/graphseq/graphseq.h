//
// graphseq.h
// 
// Copyright 2014 Vincent Q. Vu. All rights reserved
//

#ifndef __GRAPHSEQ_H
#define __GRAPHSEQ_H

#include "base.h"
#include <boost/pending/disjoint_sets.hpp>

struct GraphBlock : public arma::uvec {
public:
  GraphBlock(const GraphBlock& a, const GraphBlock& b) 
    : arma::uvec(a.n_elem + b.n_elem) 
  {
    auto i = begin();
    for (auto ai : a) { *i++ = ai; }
    for (auto bi : b) { *i++ = bi; }
  }

  GraphBlock(arma::uword i) : arma::uvec(1) { at(0) = i; }

  arma::uword size() const { return n_elem; }
  void sort() { std::sort(begin(), end()); }
};

struct GraphSeq : GraphSeqBase<arma::uword, GraphBlock> {

  GraphSeq(const arma::mat& x, double minweight, 
           arma::uword maxblocksize, arma::uword minblocknum = 1) {

    // Construct edges
    std::priority_queue<edge_t> edges;
    for (arma::uword j = 0; j < x.n_cols; ++j) {
      for (arma::uword i = j + 1; i < x.n_rows; ++i) {
        double weight = std::fabs(x(i,j));
        if (weight > minweight) {
          edges.push( std::make_pair(weight, std::make_pair(i, j)) );
        }
      }
    }

    init(x.n_cols, edges, minweight, maxblocksize, minblocknum);
  }

  // Sparse matrix constructor - expects x to contain at least a lower triangle
  GraphSeq(const arma::sp_mat& x, double minweight, 
           arma::uword maxblocksize, arma::uword minblocknum = 1) {

    // Construct edges
    std::priority_queue<edge_t> edges;
    for (auto i = x.begin(); i != x.end(); ++i) {
      double weight = std::fabs(*i);
      if (weight > minweight && i.row() > i.col()) {
        edges.push( std::make_pair(weight, std::make_pair(i.row(), i.col())) );
      }
    }

    init(x.n_cols, edges, minweight, maxblocksize, minblocknum);
  }

protected:
  void init(arma::uword n, queue_t& edges, double minweight, 
            arma::uword maxblocksize, arma::uword minblocknum = 1) 
  {
    // Initialize disjoint sets structure
    std::vector<arma::uword> rank(n);
    std::vector<arma::uword> parent(n);
    boost::disjoint_sets<arma::uword*, arma::uword*> ds(&rank[0], &parent[0]);

    // Initialize singleton partition
    for (arma::uword v = 0; v < n; ++v) {
      ds.make_set(v);
      current.insert(current.cend(), 
                     std::make_pair(ds.find_set(v), 
                                    GraphBlock(v)));
    }

    // Initialize graph sequence
    GraphSeqBase<arma::uword, GraphBlock>::init(
      ds, edges, minweight, maxblocksize, minblocknum);
  }
};


#endif
