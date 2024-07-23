#pragma once

#include "config.h"
#include "utils.h"
#include <vector>

using namespace std;

class VoxelPair {
  public:
    const vector<index_t> birth;
    const vector<index_t> death;

    VoxelPair(const vector<index_t> &birth, const vector<index_t> &death);
    template <typename... Types>
    VoxelPair(const std::tuple<Types...> &_birth,
              const std::tuple<Types...> &_death)
        : birth(tupleToVector(_birth)), death(tupleToVector(_death)) {}

    void print() const;
};

class VoxelMatch {
  public:
    const VoxelPair pair0;
    const VoxelPair pair1;

    VoxelMatch(const VoxelPair &pair0, const VoxelPair &pair1);
    void print() const;
};
