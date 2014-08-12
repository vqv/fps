//
// map.h
// 
// Copyright 2014 Vincent Q. Vu. All rights reserved
//

#include <map>

template <typename Key, typename BM = BlockMat<arma::mat> >
class map : public BM {

public:
  typedef std::pair<arma::uvec, arma::uvec> index_t;
  typedef typename std::map<Key, index_t> indexmap_t;

  map(const typename BM::block_type& X, const indexmap_t& indexmap) {
    for (const auto& i : indexmap) {
      index_t index = i.second;
      if(index.first.n_elem == 0 || index.second.n_elem == 0) { continue; }

      typename BM::block_type b = X.submat(index.first, index.second);
      BM::blocks.push_back(std::move(b));
      indices.push_back(std::move(index));
    }
  }

  void copy_to(typename BM::block_type& X) const {
    auto b = BM::blocks.cbegin();
    for (auto& i : indices) { X.submat(i.first, i.second) = *b++; }
  }

protected:
  std::deque<index_t> indices;
};

template <typename Key, typename BM = BlockMat<arma::mat> >
class symmap : public BM {

public:
  typedef arma::uvec index_t;
  typedef typename std::map<Key, index_t> indexmap_t;

  symmap(const typename BM::block_type& X, const indexmap_t& indexmap) {
    for (const auto& i : indexmap) {
      index_t index = i.second;
      if(index.n_elem == 0) { continue; }

      typename BM::block_type b = X.submat(index, index);
      BM::blocks.push_back(std::move(b));
      indices.push_back(std::move(index));
    }
  }

  void copy_to(typename BM::block_type& X) const {
    auto b = BM::blocks.cbegin();
    for (auto& i : indices) { X.submat(i, i) = *b++; }
  }

protected:
  std::deque<index_t> indices;
};
