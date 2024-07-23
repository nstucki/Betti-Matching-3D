#pragma once

#include "../config.h"

#include <cstdint>
#include <vector>

using namespace std;

namespace dim1 {
typedef tuple<index_t> Coordinate;
typedef vector<Coordinate> RepresentativeCycle;

class Cube {
  public:
    value_t birth;
    index_t index;

    Cube();
    Cube(const Cube &cube);
    Cube(value_t birth, index_t x);
    index_t x() const;
    bool operator==(const Cube &rhs) const;
    void print() const;
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
    const Cube birth;
    const Cube death;
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
    value_t getBirth(const index_t &x) const;
    value_t getBirth(const index_t &x, const uint8_t &dim) const;
    index_t getParentVoxel(const Cube &c, const uint8_t &dim) const;
    void printImage() const;
    void
    printRepresentativeCycle(const dim1::RepresentativeCycle &reprCycle) const;
    const vector<index_t> shape;
    const index_t m_x;

  private:
    value_t *allocateMemory() const;
    void getGridFromVector(const vector<value_t> &vector);
    value_t *grid;
};

class UnionFind {
  public:
    UnionFind(const CubicalGridComplex &cgc);
    index_t find(index_t x);
    index_t link(index_t x, index_t y);
    value_t getBirth(const index_t &idx) const;
    index_t getCoordinate(index_t idx) const;
    vector<index_t> getBoundaryIndices(const Cube &edge) const;
    void reset();

  private:
    vector<index_t> parent;
    vector<value_t> birthtime;
    const CubicalGridComplex &cgc;
};
} // namespace dim1
