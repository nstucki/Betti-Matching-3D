#pragma once

#include "../config.h"

#include <cstddef> // For std::ptrdiff_t
#include <cstdint>
#include <iostream>
#include <iterator> // For std::forward_iterator_tag
#include <optional>
#include <tuple>
#include <vector>

#include <sstream>

using namespace std;

namespace dim3 {
typedef tuple<index_t, index_t, index_t> Coordinate;
typedef vector<Coordinate> RepresentativeCycle;

class Cube {
  public:
    Cube();
    Cube(value_t birth, index_t x, index_t y, index_t z, uint8_t type);
    Cube(value_t birth, vector<index_t> coordinates, uint8_t type);
    Cube(const Cube &cube);
    bool operator==(const Cube &rhs) const;
    index_t x() const;
    index_t y() const;
    index_t z() const;
    uint8_t type() const;
    void print() const;
    value_t birth;
    uint64_t index;
};

struct CubeComparator {
    bool operator()(const Cube &Cube1, const Cube &Cube2) const;
};

class Pair {
  public:
    Pair();
    Pair(const Cube &birth, const Cube &death);
    Pair(const Pair &pair);
    bool operator==(const Pair &rhs) const;
    void print() const;
    Cube birth;
    Cube death;
};

class Match {
  public:
    Match(Pair pair0, Pair pair1);
    void print() const;
    Pair pair0;
    Pair pair1;
};

class CubicalGridComplex {
  public:
    CubicalGridComplex(const vector<value_t> &image,
                       const vector<index_t> &shape);
    CubicalGridComplex(CubicalGridComplex &&other);
    ~CubicalGridComplex();
    size_t getNumberOfCubes(const uint8_t &dim) const;
    value_t getBirth(const index_t &x, const index_t &y,
                     const index_t &z) const;
    value_t getBirth(const index_t &x, const index_t &y, const index_t &z,
                     const uint8_t &type, const uint8_t &dim) const;
    Coordinate getParentVoxel(const Cube &c, const uint8_t &dim) const;
    void printImage() const;
    void printRepresentativeCycle(const RepresentativeCycle &reprCycle) const;
    const vector<index_t> shape;
    const index_t m_x;
    const index_t m_y;
    const index_t m_z;
    const index_t m_yz;
    const index_t m_xyz;
    const index_t n_yz;
    const index_t n_xyz;

  private:
    value_t ***allocateMemory() const;
    void getGridFromVector(const vector<value_t> &vector);
    value_t ***grid;
};

class UnionFind {
  public:
    UnionFind(const CubicalGridComplex &cgc);
    index_t find(index_t x);
    index_t link(index_t x, index_t y);
    value_t getBirth(const index_t &idx) const;
    Coordinate getCoordinates(index_t idx) const;
    vector<index_t> getBoundaryIndices(const Cube &edge) const;
    void reset();

  private:
    vector<index_t> parent;
    vector<value_t> birthtime;
    const CubicalGridComplex &cgc;
};

class UnionFindDual {
  public:
    UnionFindDual(const CubicalGridComplex &cgc);
    index_t find(index_t x);
    index_t link(index_t x, index_t y);
    value_t getBirth(const index_t &idx) const;
    Coordinate getCoordinates(index_t idx) const;
    vector<index_t> getBoundaryIndices(const Cube &edge) const;
    void reset();

  private:
    vector<index_t> parent;
    vector<value_t> birthtime;
    const CubicalGridComplex &cgc;
};

template <int _dim, class _Tp> class CubeMap {
    /// Datastructure to map cube indices to values efficiently via a
    /// 1-dimensional array by mapping cube indices to (x,y,z,type) coordinates
    /// in a 4-dimensional space and representing this as a 1-dimensional array.

  public:
    CubeMap(vector<index_t> shape);
    void emplace(uint64_t cube_index, _Tp element);
    const std::optional<_Tp> &find(uint64_t cube_index) const;
    void clear();
    optional<_Tp> &operator[](uint64_t cube_index);

  private:
    vector<std::optional<_Tp>> elements;
    uint64_t computeCoordinateIndex(uint64_t cube_index) const;
    vector<index_t> shape;
    std::optional<_Tp> none = {};
    static const int NUM_TYPES = (_dim == 1 || _dim == 2) ? 3 : 1;
    const int stride_x_direction;
    const int stride_y_direction;
    const int stride_z_direction;
};

template <int _dim, class _Tp>
CubeMap<_dim, _Tp>::CubeMap(vector<index_t> shape)
    : shape(shape), stride_x_direction(shape[1] * shape[2] * NUM_TYPES),
      stride_y_direction(shape[2] * NUM_TYPES), stride_z_direction(NUM_TYPES),
      elements(shape[0] * shape[1] * shape[2] * NUM_TYPES) {}

template <int _dim, class _Tp>
void CubeMap<_dim, _Tp>::emplace(uint64_t cube_index, _Tp element) {
    if (cube_index != NONE_INDEX) {
        elements[computeCoordinateIndex(cube_index)] = element;
    } else {
        throw runtime_error(
            "CubeMap::emplace may not be called with NONE_INDEX magic number");
    }
}

template <int _dim, class _Tp>
const std::optional<_Tp> &CubeMap<_dim, _Tp>::find(uint64_t cube_index) const {
    if (cube_index != NONE_INDEX) {
        return elements[computeCoordinateIndex(cube_index)];
    }
    return none;
}

template <int _dim, class _Tp>
optional<_Tp> &CubeMap<_dim, _Tp>::operator[](uint64_t cube_index) {
    if (cube_index != NONE_INDEX) {
        return elements[computeCoordinateIndex(cube_index)];
    }
    throw runtime_error("CubeMap subscript operator may not be called with "
                        "NONE_INDEX magic number");
}

template <int _dim, class _Tp>
uint64_t CubeMap<_dim, _Tp>::computeCoordinateIndex(uint64_t cube_index) const {
    int x = (cube_index >> 44) & 0xfffff;
    int y = (cube_index >> 24) & 0xfffff;
    int z = (cube_index >> 4) & 0xfffff;
    int type = (_dim == 1 || _dim == 2) ? (cube_index & 0xf) : 0;

    return x * stride_x_direction + y * stride_y_direction +
           z * stride_z_direction + type;
}

template <int _dim, class _Tp> void CubeMap<_dim, _Tp>::clear() {
    elements.clear();
    elements.resize(shape[0] * shape[1] * shape[2] * NUM_TYPES);
}
} // namespace dim3
