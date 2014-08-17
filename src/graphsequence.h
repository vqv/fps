//
// graphsequence.h
// 
// Copyright 2014 Vincent Q. Vu. All rights reserved
//

#ifndef __GRAPHSEQUENCE_H
#define __GRAPHSEQUENCE_H

#include <RcppArmadillo.h>
#include <queue>
#include <map>
#include <utility>
#include <algorithm>
#include <boost/pending/disjoint_sets.hpp>

// A map representing a sequence of vertex partitions induced by 
// the connected components that form as edges are added in decreasing 
// order of weight.
// A component is represented by std::set<V>
// A partition is represented by std::map<V, std::set<V>>
template <typename Vertex, typename Block>
struct GraphSeqBase {
  typedef Vertex vertex_t;
  typedef Block block_t;
  typedef std::map<Vertex, Block> partition_t;
  typedef std::map<double, partition_t, std::greater<double>> sequence_t;

  // Expose a subset of the map interface to sequence
  typedef typename sequence_t::const_iterator const_iterator;
  typedef typename sequence_t::const_reverse_iterator const_reverse_iterator;
  typedef typename sequence_t::value_type value_type;
  typedef typename sequence_t::key_type key_type;
  typedef typename sequence_t::mapped_type mapped_type;
  typedef typename sequence_t::size_type size_type;

  GraphSeqBase() : 
    current_max(0),
    last_weight(std::numeric_limits<double>::infinity())
    {}

  size_type size() const { return sequence.size(); };

  const_iterator begin() const { return sequence.begin(); }
  const_iterator end() const { return sequence.end(); }
  const_iterator cbegin() const { return sequence.cbegin(); }
  const_iterator cend() const { return sequence.cend(); }
  const_reverse_iterator rbegin() const { return sequence.rbegin(); }
  const_reverse_iterator rend() const { return sequence.rend(); }
  const_reverse_iterator crbegin() const { return sequence.crbegin(); }
  const_reverse_iterator crend() const { return sequence.crend(); }

  // The first partition whose weight is not > is active
  const partition_t& get_active(const double weight) const {
    return sequence.lower_bound(weight)->second;
  }

protected:
  typedef std::pair<double, std::pair<Vertex, Vertex>> edge_t;
  typedef std::priority_queue<edge_t> queue_t;

  sequence_t sequence;

  partition_t current;
  arma::uword current_max;
  std::set<Vertex> newblocks;
  double last_weight;

  void flush_blocks(const double weight) {
    for (const auto& m : newblocks) {
      auto p = current.find(m);
      if (p != current.end()) { 
        p->second.sort();
      }
    }
    newblocks.clear();
    sequence.insert(sequence.cend(), std::make_pair(weight, current));
  }

  void merge_blocks(const double weight, const Vertex& a, const Vertex& b, 
                    const Vertex& c) {

    // New knot so store the current partition
    if (weight < last_weight) {
      flush_blocks(weight);
      last_weight = weight;
    }

    auto pa = current.find(a);
    auto pb = current.find(b);

    Block newblock(pa->second, pb->second);
    if (newblock.size() > current_max) { current_max = newblock.size(); }

    current.erase(pa);
    current.erase(pb);
    current.insert(std::make_pair(c, std::move(newblock)));
  }

  // Returns a set of identifiers for blocks that merged. 
  // Note that these blocks could possibly merge into more than 
  // one new block.
  template <typename DisjointSets>
  void merge(DisjointSets& ds, queue_t& edges) {

    Vertex a, b, c;
    double weight;

    // Pop edges until we find an edge between two components
    while (!edges.empty()) {
      a = ds.find_set( edges.top().second.first );
      b = ds.find_set( edges.top().second.second );
      if (a == b) {
        edges.pop();
        continue;
      }

      weight = edges.top().first;

      // Merge components and pop the edge
      ds.link(a, b);
      c = ds.find_set(edges.top().second.first);
      merge_blocks(weight, a, b, c);
      edges.pop();

      // Add all edges with the same weight in case of ties
      while (!edges.empty() && edges.top().first == weight) {
        a = ds.find_set(edges.top().second.first);
        b = ds.find_set(edges.top().second.second);
        if (a != b) {
          ds.link(a, b);
          c = ds.find_set(edges.top().second.first);
          merge_blocks(weight, a, b, c);
        }
        edges.pop();
      }
      break;
    }
  }

  template <typename DisjointSets>
  void init(DisjointSets& ds, queue_t& edges, double minweight, 
            arma::uword maxblocksize, arma::uword minblocknum = 1) {

    // Merge components until there is only one
    while (!edges.empty() && current.size() > minblocknum && 
           current_max < maxblocksize) {
      merge(ds, edges);
    }

    // Store the final partition
    if (edges.empty() || current.size() == 1) {
      flush_blocks(minweight);
    } else {
      flush_blocks(std::nexttoward(last_weight, 0.0));
    }
  }

};

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
  void sort() { 
    std::sort(first.begin(), first.end());
    std::sort(second.begin(), second.end());
  }
};

struct BiGraphSeq : GraphSeqBase<BiGraphVertex, BiGraphBlock> {
public:
  BiGraphSeq(const arma::mat& x, double minweight, 
             arma::uword maxblocksize, arma::uword minblocknum = 1) {

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
             arma::uword maxblocksize, arma::uword minblocknum = 1) {

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
            arma::uword minblocknum = 1) 
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
    for (arma::uword i = 0; i < n_rows; ++i) {
      BiGraphVertex v(true, i);
      ds.make_set(v);
      current.insert(current.cend(), 
                     std::make_pair(ds.find_set(v), BiGraphBlock(v)));
    }
    for (arma::uword j = 0; j < n_cols; ++j) {
      BiGraphVertex v(false, j);
      ds.make_set(v);
      current.insert(current.cend(), 
                     std::make_pair(ds.find_set(v), BiGraphBlock(v)));
    }

    // Initialize graph sequence
    GraphSeqBase<BiGraphVertex, BiGraphBlock>::init(
      ds, edges, minweight, maxblocksize, minblocknum);

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
