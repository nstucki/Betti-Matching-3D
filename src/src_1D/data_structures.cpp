#include "data_structures.h"
#include <algorithm>
#include <iostream>

using namespace dim1;
using namespace std;

Cube::Cube() : birth(0), index(NONE) {}

Cube::Cube(const Cube &cube) : birth(cube.birth), index(cube.index) {}

Cube::Cube(value_t _birth, index_t x) : birth(_birth), index(x) {}

bool Cube::operator==(const Cube &rhs) const { return (index == rhs.index); }

index_t Cube::x() const { return index; }

void Cube::print() const { cout << "(" << birth << "," << x() << ")"; }

bool CubeComparator::operator()(const Cube &cube1, const Cube &cube2) const {
    if (cube1.birth == cube2.birth) {
        return (cube1.index < cube2.index);
    } else {
        return (cube1.birth < cube2.birth);
    }
}

Pair::Pair() {}

Pair::Pair(const Cube &_birth, const Cube &_death)
    : birth(_birth), death(_death) {}

Pair::Pair(const Pair &pair) : birth(pair.birth), death(pair.death) {}

bool Pair::operator==(const Pair &rhs) const {
    return (birth == rhs.birth && death == rhs.death);
}

void Pair::print() const {
    cout << "(";
    birth.print();
    cout << ";";
    death.print();
    cout << ")";
}

Match::Match(Pair _pair0, Pair _pair1) : pair0(_pair0), pair1(_pair1) {}

void Match::print() const {
    pair0.print();
    cout << " <-> ";
    pair1.print();
    cout << endl;
}

CubicalGridComplex::CubicalGridComplex(const vector<value_t> &image,
                                       const vector<index_t> &_shape)
    : shape(_shape), m_x(shape[0] - 1) {
    getGridFromVector(image);
}

CubicalGridComplex::CubicalGridComplex(CubicalGridComplex &&other)
    : m_x(other.m_x), shape(std::move(other.shape)) {
    grid = other.grid;
    other.grid = nullptr;
}

CubicalGridComplex::~CubicalGridComplex() {
    if (grid != nullptr) {
        delete[] grid;
    }
}

size_t CubicalGridComplex::getNumberOfCubes(const uint8_t &dim) const {
    switch (dim) {
    case 0:
        return shape[0];

    case 1:
        return m_x;
    }
    throw runtime_error("No cubes in dim " + std::to_string(unsigned(dim)));
}

value_t CubicalGridComplex::getBirth(const index_t &x) const {
    return grid[x + 1];
}

value_t CubicalGridComplex::getBirth(const index_t &x,
                                     const uint8_t &dim) const {
    switch (dim) {
    case 0:
        return getBirth(x);

    case 1:
        return max(getBirth(x), getBirth(x + 1));
    }
    throw runtime_error("Birth not found!");
}

index_t CubicalGridComplex::getParentVoxel(const Cube &cube,
                                           const uint8_t &dim) const {
    index_t x = cube.x();
    switch (dim) {
    case 0:
        return x;

    case 1:
        if (cube.birth == getBirth(x + 1)) {
            return x + 1;
        } else {
            return x;
        }
    }
    throw runtime_error("Parent voxel not found!");
}

void CubicalGridComplex::printImage() const {
    value_t birth;
    for (index_t x = 0; x < shape[0]; ++x) {
        birth = getBirth(x);
        if (birth < 10) {
            cout << ' ' << birth << ' ';
        } else {
            cout << birth << ' ';
        }
    }
    cout << endl << endl;
}

void CubicalGridComplex::printRepresentativeCycle(
    const dim1::RepresentativeCycle &reprCycle) const {
    for (index_t x = 0; x < shape[0]; ++x) {
        auto it = find(reprCycle.begin(), reprCycle.end(), dim1::Coordinate{x});
        if (it == reprCycle.begin()) {
            cout << "2  ";
        } else if (it == reprCycle.end() - 1) {
            cout << "-1 ";
        } else if (it != reprCycle.end()) {
            cout << "1  ";
        } else {
            cout << "0  ";
        }
    }
}

value_t *CubicalGridComplex::allocateMemory() const {
    value_t *g = new value_t[shape[0] + 2];
    if (g == NULL) {
        throw runtime_error("Out of memory!");
    }
    return g;
}

void CubicalGridComplex::getGridFromVector(const vector<value_t> &vec) {
    size_t counter = 0;
    grid = allocateMemory();
    for (index_t x = 0; x < shape[0] + 2; ++x) {
        if (x == 0 || x == shape[0] + 1) {
            grid[x] = INFTY;
        } else {
            grid[x] = vec[counter++];
        }
    }
}

UnionFind::UnionFind(const CubicalGridComplex &_cgc) : cgc(_cgc) {
    size_t n = cgc.getNumberOfCubes(0);
    parent.reserve(n);
    birthtime.reserve(n);
    for (index_t x = 0; x < _cgc.shape[0]; ++x) {
        parent.push_back(x);
        birthtime.push_back(cgc.getBirth(x));
    }
}

index_t UnionFind::find(index_t x) {
    index_t y = x, z = parent[y];
    while (z != y) {
        y = z;
        z = parent[y];
    }
    y = parent[x];
    while (z != y) {
        parent[x] = z;
        x = y;
        y = parent[x];
    }
    return z;
}

index_t UnionFind::link(index_t x, index_t y) {
    if (birthtime[x] > birthtime[y]) {
        parent[x] = y;
        return x;
    } else if (birthtime[x] < birthtime[y]) {
        parent[y] = x;
        return y;
    } else {
        if (x > y) {
            parent[x] = y;
            return x;
        } else {
            parent[y] = x;
            return y;
        }
    }
}

value_t UnionFind::getBirth(const index_t &idx) const { return birthtime[idx]; }

index_t UnionFind::getCoordinate(index_t idx) const { return idx; }

vector<index_t> UnionFind::getBoundaryIndices(const Cube &edge) const {
    return {edge.x(), edge.x() + 1};
}

void UnionFind::reset() {
    for (size_t i = 0; i < parent.size(); ++i) {
        parent[i] = i;
    }
}
