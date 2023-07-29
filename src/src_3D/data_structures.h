#pragma once

#include "config.h"

#include <cstdint>
#include <vector>

using namespace std;


class Cube {
    public:
	value_t birth;
	uint64_t index;

	Cube();
    Cube(value_t birth, index_t x, index_t y, index_t z, uint8_t type);
	Cube(const Cube& cube);
	index_t x() const;
	index_t y() const;
	index_t z() const;
	uint8_t type() const;
    bool operator==(const Cube& rhs) const;
	void print() const;
};


struct CubeComparator{ bool operator()(const Cube& Cube1, const Cube& Cube2) const; };


class Pair {
    public:
	const Cube birth;
	const Cube death;

	Pair();
    Pair(const Cube& birth, const Cube& death);
	Pair(const Pair& pair);
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
	const vector<index_t> shape;
	const index_t n_x;
	const index_t n_y;
	const index_t n_z;

	const index_t n_yz;
	const index_t n_xyz;
	const index_t m_yz;
	const index_t m_xyz;

	// vector
	CubicalGridComplex(const vector<value_t>& image, const vector<index_t>& shape);
	~CubicalGridComplex();
	index_t getNumberOfCubes(uint8_t dim) const;
	value_t getBirth(index_t x, index_t y, index_t z) const;
	value_t getBirth(index_t x, index_t y, index_t z, uint8_t type, uint8_t dim) const;
	vector<index_t> getParentVoxel(const Cube& c, uint8_t dim) const;
	void printImage() const;

	private:
	value_t*** grid;

	value_t*** allocateMemory() const;
	void getGridFromVector(const vector<value_t> vector);
};


class UnionFind {
    public:
	UnionFind(const CubicalGridComplex* const cgc);
	index_t find(index_t x);
	index_t link(index_t x, index_t y);
	value_t getBirth(index_t x) const;
	vector<index_t> getCoordinates(index_t x) const;
	vector<index_t> getBoundaryIndices(const Cube& edge) const;
	void reset();

    private:
	vector<index_t> parent;
	vector<value_t> birthtime;
	const CubicalGridComplex* const cgc;
};

class UnionFindDual {
    public:
	UnionFindDual(const CubicalGridComplex* const cgc);
	index_t find(index_t x);
	index_t link(index_t x, index_t y);
	value_t getBirth(index_t x) const;
	vector<index_t> getCoordinates(index_t x) const;
	vector<index_t> getBoundaryIndices(const Cube& edge) const;
	void reset();

    private:
	vector<index_t> parent;
	vector<value_t> birthtime;
	const CubicalGridComplex* const cgc;
};