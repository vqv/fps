//
// bigraphseq.h
// 
// Copyright 2014 Vincent Q. Vu. All rights reserved
//

#ifndef __BIGRAPHSEQ_H
#define __BIGRAPHSEQ_H

#include "base.h"
#include <boost/pending/disjoint_sets.hpp>

// (row?, index)
typedef std::pair<bool, arma::uword> BiGraphVertex;

struct BiGraphBlock {

  arma::uvec first, second;

  BiGraphBlock(const BiGraphVertex& v) {
    arma::uvec *i = v.first ? &first : &second;
    i->set_size(1);
    i->at(0) = v.second;
  }

  BiGraphBlock(const BiGraphBlock& a, const BiGraphBlock& b) 
    : first(a.first.n_elem + b.first.n_elem),
      second(a.second.n_elem + b.second.n_elem) 
  {
    arma::uvec::iterator i;

    i = first.begin();
    for (auto ai : a.first) { *i++ = ai;}
    for (auto bi : b.first) { *i++ = bi;}

    i = second.begin();
    for (auto ai : a.second) { *i++ = ai;}
    for (auto bi : b.second) { *i++ = bi;}
  }

  arma::uword size() const { return first.n_elem + second.n_elem; }
};

struct BiGraphSeq : GraphSeqBase<BiGraphVertex, BiGraphBlock> {
public:
  BiGraphSeq(const arma::mat& x, double minweight, 
             arma::uword maxblocksize, arma::uword minblocknum) {

    // Construct edges
    std::priority_queue<edge_t> edges;
    for (arma::uword j = 0; j < x.n_cols; ++j) {
      for (arma::uword i = 0 ; i < x.n_rows; ++i) {
        double weight = std::fabs(x(i,j));
        if (weight > minweight) {
          edges.push(std::make_pair(weight, 
                      std::make_pair(BiGraphVertex(true, i), 
                                     BiGraphVertex(false, j)) ));
        }
      }
    }

    init(x.n_rows, x.n_cols, edges, minweight, maxblocksize, minblocknum);
  }

  // Sparse matrix constructor
  BiGraphSeq(const arma::sp_mat& x, double minweight, 
             arma::uword maxblocksize, arma::uword minblocknum) {

    // Construct edges
    std::priority_queue<edge_t> edges;
    for (auto i = x.begin(); i != x.end(); ++i) {
      double weight = std::fabs(*i);
      if (weight > minweight) {
        edges.push(std::make_pair(weight, 
                    std::make_pair(BiGraphVertex(true, i.row()), 
                                   BiGraphVertex(false, i.col())) ));        
      }
    }

    init(x.n_rows, x.n_cols, edges, minweight, maxblocksize, minblocknum);
  }

protected:
  void init(arma::uword n_rows, arma::uword n_cols, queue_t& edges, 
            double minweight, arma::uword maxblocksize, 
            arma::uword minblocknum) 
  {
    // Initialize disjoint sets structure
    typedef std::map<BiGraphVertex, std::size_t> rank_t;
    typedef std::map<BiGraphVertex, BiGraphVertex> parent_t;
    rank_t rank_map;
    parent_t parent_map;
    typedef boost::associative_property_map<rank_t> Rank;
    typedef boost::associative_property_map<parent_t> Parent;
    Rank rank_pmap(rank_map);
    Parent parent_pmap(parent_map);
    boost::disjoint_sets<Rank, Parent> ds(rank_pmap, parent_pmap);

    // Initialize singleton partition
    partition_t &singletons = sequence.begin()->second;
    for (arma::uword i = 0; i < n_rows; ++i) {
      BiGraphVertex v(true, i);
      ds.make_set(v);
      singletons.insert(singletons.cend(), 
                     std::make_pair(ds.find_set(v), BiGraphBlock(v)));
    }
    for (arma::uword j = 0; j < n_cols; ++j) {
      BiGraphVertex v(false, j);
      ds.make_set(v);
      singletons.insert(singletons.cend(), 
                     std::make_pair(ds.find_set(v), BiGraphBlock(v)));
    }

    // Initialize graph sequence
    GraphSeqBase<BiGraphVertex, BiGraphBlock>::init(
      ds, edges, maxblocksize, minblocknum);

    // Remove singletons
    for (auto& g : sequence) {
      for (auto i = g.second.begin(); i != g.second.end(); ) {
        if (i->second.first.n_elem == 0 || i->second.second.n_elem == 0) {
          i = g.second.erase(i);
        } else {
          ++i;
        }
      }
    }

  }
};

#endif
