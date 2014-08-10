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
#include <boost/pending/disjoint_sets.hpp>
#include <sstream>

template <typename Derived>
struct GraphSeqBase_traits;

// A map representing a sequence of vertex partitions induced by 
// the connected components that form as edges are added in decreasing 
// order of weight.
// A component is represented by std::set<V>
// A partition is represented by std::map<V, std::set<V>>
template <typename Derived>
struct GraphSeqBase {
  typedef typename GraphSeqBase_traits<Derived>::vertex_t vertex_t;
  typedef std::pair<double, std::pair<vertex_t, vertex_t>> edge_t;

  typedef typename GraphSeqBase_traits<Derived>::block_t block_t;
  typedef std::map<vertex_t, block_t> partition_t;
  typedef std::map<double, partition_t, std::greater<double>> sequence_t;

  // The first partition whose weight is not > is active
  const partition_t& get_active(const double weight) const {
    return sequence.lower_bound(weight)->second;
  }

  // Expose a subset of the map interface to sequence
  typedef typename sequence_t::const_iterator const_iterator;
  typedef typename sequence_t::const_reverse_iterator const_reverse_iterator;
  typedef typename sequence_t::value_type value_type;
  typedef typename sequence_t::key_type key_type;
  typedef typename sequence_t::mapped_type mapped_type;
  typedef typename sequence_t::size_type size_type;

  size_type size() const { return sequence.size(); };

  mapped_type& operator[](const key_type& key) { return sequence[key]; }
  mapped_type& operator[](key_type&& key) { return sequence[key]; }
  const_iterator begin() const { return sequence.begin(); }
  const_iterator end() const { return sequence.end(); }
  const_iterator cbegin() const { return sequence.cbegin(); }
  const_iterator cend() const { return sequence.cend(); }
  const_reverse_iterator rbegin() const { return sequence.rbegin(); }
  const_reverse_iterator rend() const { return sequence.rend(); }
  const_reverse_iterator crbegin() const { return sequence.crbegin(); }
  const_reverse_iterator crend() const { return sequence.crend(); }

protected:
  sequence_t sequence;

  typedef std::priority_queue<edge_t> queue_t;

  // Returns a set of identifiers for blocks that merged. 
  // Note that these blocks could possibly merge into more than 
  // one new block.
  template <typename DisjointSets>
  std::set<vertex_t> merge(DisjointSets& ds, queue_t& edges) {

    vertex_t a, b, key;
    double weight;
    std::set<vertex_t> merged;

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
      key = ds.find_set(edges.top().second.first);
      static_cast<Derived*>(this)->merge_blocks(weight, a, b, key);
      merged.insert(a);
      merged.insert(b);
      merged.insert(key);
      edges.pop();

      // Add all edges with the same weight in case of ties
      while (!edges.empty() && edges.top().first == weight) {
        a = ds.find_set(edges.top().second.first);
        b = ds.find_set(edges.top().second.second);
        if (a != b) {
          ds.link(a, b);
          key = ds.find_set(edges.top().second.first);
          static_cast<Derived*>(this)->merge_blocks(weight, a, b, key);
          merged.insert(a);
          merged.insert(b);
          merged.insert(key);
        }
        edges.pop();
      }
      break;
    }

    return merged;
  }
};

struct GraphSeq;

template <>
struct GraphSeqBase_traits<GraphSeq> {
public:
  typedef arma::uword vertex_t;
  typedef arma::uvec block_t;
};

struct GraphSeq : GraphSeqBase<GraphSeq> {
public:
  friend struct GraphSeqBase<GraphSeq>;
  typedef typename GraphSeqBase_traits<GraphSeq>::vertex_t vertex_t;
  typedef typename GraphSeqBase_traits<GraphSeq>::block_t block_t;

  GraphSeq(const arma::mat& x, double lambdamin = 0.0) {

    // Construct edges
    std::priority_queue<edge_t> edges;
    for (arma::uword i = 0; i < x.n_rows; ++i) {
      for (arma::uword j = i + 1; j < x.n_cols; ++j) {
        double weight = std::fabs(x(i,j));
        if (weight > lambdamin) {
          edges.emplace(weight, std::make_pair(i, j));
        }
      }
    }

    // Initialize disjoint sets structure
    std::vector<arma::uword> rank(x.n_cols);
    std::vector<arma::uword> parent(x.n_cols);
    boost::disjoint_sets<arma::uword*, arma::uword*> ds(&rank[0], &parent[0]);

    // Initialize singleton partition
    for (arma::uword v = 0; v < x.n_cols; ++v) {
      arma::uvec b(1);
      b[0] = v;
      ds.make_set(v);
      current.emplace_hint(current.cend(), ds.find_set(v), std::move(b));
    }
    last_weight = std::numeric_limits<double>::infinity();

    // Merge components until there is only one
    while (!edges.empty() && sequence.crbegin()->second.size() > 1) {
      std::set<vertex_t> merged = merge(ds, edges);
      for(const auto& m : merged) {
        partition_t::iterator p = sequence.rbegin()->second.find(m);
        if(p == sequence.rbegin()->second.end()) { continue; }
        p->second = arma::sort(p->second);
      }
    }

    // Store the final partition
    sequence.emplace_hint(sequence.cend(), lambdamin, current);
  }

protected:
  partition_t current;
  double last_weight;

  void merge_blocks(const double& weight, 
                    const vertex_t& a, const vertex_t& b, 
                    const vertex_t& key) {

    // New knot so store the current partition
    if (weight < last_weight) {
      sequence.emplace_hint(sequence.cend(), weight, current);
      last_weight = weight;
    }

    partition_t::iterator pa, pb;

    pa = current.find(a), 
    pb = current.find(b);
    block_t newblock = arma::join_vert(pa->second, pb->second);

    current.erase(pa);
    current.erase(pb);
    current.emplace(key, std::move(newblock));
  }
};

struct BiGraphSeq;

template <>
struct GraphSeqBase_traits<BiGraphSeq> {
public:
  typedef std::pair<arma::uword, arma::uword> vertex_t;
  typedef std::pair<arma::uvec, arma::uvec> block_t;
};

struct BiGraphSeq : GraphSeqBase<BiGraphSeq> {
public:
  friend struct GraphSeqBase<BiGraphSeq>;
  typedef typename GraphSeqBase_traits<BiGraphSeq>::vertex_t vertex_t;
  typedef typename GraphSeqBase_traits<BiGraphSeq>::block_t block_t;

  BiGraphSeq(const arma::mat& x, double lambdamin = 0.0) {

    // Construct edges
    std::priority_queue<edge_t> edges;
    for (arma::uword i = 0; i < x.n_rows; ++i) {
      for (arma::uword j = 0 ; j < x.n_cols; ++j) {
        double weight = std::fabs(x(i,j));
        if (weight > lambdamin) {
          edges.emplace(weight, 
            std::make_pair(std::make_pair(0, i), std::make_pair(1, j))
          );
        }
      }
    }

    // Initialize disjoint sets structure
    typedef std::map<vertex_t, std::size_t> rank_t;
    typedef std::map<vertex_t, vertex_t> parent_t;
    rank_t rank_map;
    parent_t parent_map;
    typedef boost::associative_property_map<rank_t> Rank;
    typedef boost::associative_property_map<parent_t> Parent;
    Rank rank_pmap(rank_map);
    Parent parent_pmap(parent_map);
    boost::disjoint_sets<Rank, Parent> ds(rank_pmap, parent_pmap);

    // Initialize singleton partition
    for (arma::uword i = 0; i < x.n_rows; ++i) {
      vertex_t v = std::make_pair(0, i);
      block_t b;
      b.first.set_size(1);
      b.first[0] = i;
      ds.make_set(v);
      current.emplace_hint(current.cend(), ds.find_set(v), std::move(b));
    }
    for (arma::uword j = 0; j < x.n_cols; ++j) {
      vertex_t v = std::make_pair(1, j);
      block_t b;
      b.second.set_size(1);
      b.second[0] = j;
      ds.make_set(v);
      current.emplace_hint(current.cend(), ds.find_set(v), std::move(b));
    }
    last_weight = std::numeric_limits<double>::infinity();

    // Merge components until there is only one
    while (!edges.empty() && sequence.crbegin()->second.size() > 1) {
      std::set<vertex_t> merged = merge(ds, edges);
      for(const auto& m : merged) {
        partition_t::iterator p = sequence.rbegin()->second.find(m);
        if(p == sequence.rbegin()->second.end()) { continue; }
        p->second.first = arma::sort(p->second.first);
        p->second.second = arma::sort(p->second.second);
      }
    }

    // Store the final partition
    sequence.emplace_hint(sequence.cend(), lambdamin, current);
  }

protected:


  partition_t current;
  double last_weight;

  void merge_blocks(const double& weight, 
                    const vertex_t& a, const vertex_t& b, 
                    const vertex_t& key) {

    // New knot so store the current partition
    if (weight < last_weight) {
      sequence.emplace_hint(sequence.cend(), weight, current);
      last_weight = weight;
    }

    block_t newblock;
    partition_t::iterator pa, pb;

    pa = current.find(a);
    pb = current.find(b);
    newblock.first = arma::join_vert(pa->second.first, pb->second.first);
    newblock.second = arma::join_vert(pa->second.second, pb->second.second);

    current.erase(pa);
    current.erase(pb);
    current.emplace(key, std::move(newblock));
  }
};

#endif
