#pragma once

#include "../config.h"

#include <cstdint>
#include <vector>

using namespace std;

namespace dim2 {
typedef tuple<index_t, index_t> Coordinate;
typedef vector<Coordinate> RepresentativeCycle;

class Cube {
  public:
    value_t birth;
    uint64_t index;

    Cube();
    Cube(value_t birth, index_t x, index_t y, uint8_t type);
    Cube(const Cube &cube);
    index_t x() const;
    index_t y() const;
    uint8_t type() const;
    bool operator==(const Cube &rhs) const;
    void print() const;
};

struct CubeComparator {
    bool operator()(const Cube &Cube1, const Cube &Cube2) const;
};

class Pair {
  public:
    const Cube birth;
    const Cube death;

    Pair();
    Pair(const Cube &birth, const Cube &death);
    Pair(const Pair &pair);
    bool operator==(const Pair &rhs) const;
    void print() const;
};

class Match {
  public:
    Pair pair0;
    Pair pair1;

    Match(Pair pair0, Pair pair1);
    void print() const;
};

class CubicalGridComplex {
  public:
    CubicalGridComplex(const vector<value_t> &image,
                       const vector<index_t> &shape);
    CubicalGridComplex(CubicalGridComplex &&other);
    ~CubicalGridComplex();
    size_t getNumberOfCubes(const uint8_t &dim) const;
    value_t getBirth(const index_t &x, const index_t &y) const;
    value_t getBirth(const index_t &x, const index_t &y, const uint8_t &type,
                     const uint8_t &dim) const;
    dim2::Coordinate getParentVoxel(const Cube &c, const uint8_t &dim) const;
    void printImage() const;
    void
    printRepresentativeCycle(const dim2::RepresentativeCycle &reprCycle) const;
    const vector<index_t> shape;
    const index_t m_x;
    const index_t m_y;
    const index_t m_xy;
    const index_t n_xy;

  private:
    value_t **allocateMemory() const;
    void getGridFromVector(const vector<value_t> &vector);
    value_t **grid;
};

class UnionFind {
  public:
    UnionFind(const CubicalGridComplex &cgc);
    index_t find(index_t x);
    index_t link(index_t x, index_t y);
    value_t getBirth(const index_t &idx) const;
    dim2::Coordinate getCoordinates(index_t idx) const;
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
    dim2::Coordinate getCoordinates(index_t idx) const;
    vector<index_t> getBoundaryIndices(const Cube &edge) const;
    void reset();

  private:
    vector<index_t> parent;
    vector<value_t> birthtime;
    const CubicalGridComplex &cgc;
};
} // namespace dim2
