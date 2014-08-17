//
// base.h
// 
// Copyright 2014 Vincent Q. Vu. All rights reserved
//

#ifndef __GRAPHSEQ_BASE_H
#define __GRAPHSEQ_BASE_H

#include <RcppArmadillo.h>
#include <map>
#include <queue>
#include <utility>
#include <limits>
#include <cmath>

// A map representing a sequence of vertex partitions induced by 
// the connected components that form as edges are added in decreasing 
// order of weight.
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


#endif
