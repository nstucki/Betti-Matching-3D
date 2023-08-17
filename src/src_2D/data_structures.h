#pragma once

#include "../config.h"

#include <cstdint>
#include <vector>

using namespace std;


class Cube_2D {
    public:
	value_t birth;
	uint64_t index;

	Cube_2D();
    Cube_2D(value_t birth, index_t x, index_t y, index_t z, uint8_t type);
	Cube_2D(const Cube& cube);
	index_t x() const;
	index_t y() const;
	index_t z() const;
	uint8_t type() const;
    bool operator==(const Cube& rhs) const;
	void print() const;
};


struct CubeComparator_2D { bool operator()(const Cube& Cube1, const Cube& Cube2) const; };


class Pair_2D {
    public:
	const Cube birth;
	const Cube death;

	Pair_2D();
    Pair_2D(const Cube& birth, const Cube& death);
	Pair_2D(const Pair& pair);
	bool operator==(const Pair &rhs) const;
	void print() const;
};


class Match_2D {
	public:
	Pair pair0;
	Pair pair1;

	Match_2D(Pair pair0, Pair pair1);
	void print() const;
};


class VoxelPair_2D {
	public:
	const vector<index_t> birth;
	const vector<index_t> death;

	VoxelPair_2D(const vector<index_t>& birth, const vector<index_t>& death);
	void print() const;
};


class VoxelMatch_2D {
	public:
	const VoxelPair pair0;
	const VoxelPair pair1;

	VoxelMatch_2D(const VoxelPair& pair0, const VoxelPair& pair1);
	void print() const;
};


class CubicalGridComplex_2D {
    public:
	const vector<index_t> shape;
	const index_t n_x;
	const index_t n_y;
	const index_t n_z;
	const index_t n_yz;
	const index_t n_xyz;
	const index_t m_yz;
	const index_t m_xyz;

	CubicalGridComplex_2D(const vector<value_t>& image, const vector<index_t>& shape);
	~CubicalGridComplex_2D();
	index_t getNumberOfCubes(const uint8_t& dim) const;
	value_t getBirth(const index_t& x, const index_t& y, const index_t& z) const;
	value_t getBirth(const index_t& x, const index_t& y, const index_t& z, const uint8_t& type, const uint8_t& dim) const;
	vector<index_t> getParentVoxel(const Cube& c, const uint8_t& dim) const;
	void printImage() const;

	private:
	value_t*** grid;

	value_t*** allocateMemory() const;
	void getGridFromVector(const vector<value_t> vector);
};


class UnionFind_2D {
    public:
	UnionFind_2D(const CubicalGridComplex_2D* const cgc);
	index_t find(index_t x);
	index_t link(index_t x, index_t y);
	value_t getBirth(index_t x) const;
	vector<index_t> getCoordinates(index_t x) const;
	vector<index_t> getBoundaryIndices(const Cube& edge) const;
	void reset();

    private:
	vector<index_t> parent;
	vector<value_t> birthtime;
	const CubicalGridComplex_2D* const cgc;
};

class UnionFindDual_2D {
    public:
	UnionFindDual_2D(const CubicalGridComplex_2D* const cgc);
	index_t find(index_t x);
	index_t link(index_t x, index_t y);
	value_t getBirth(index_t x) const;
	vector<index_t> getCoordinates(index_t x) const;
	vector<index_t> getBoundaryIndices(const Cube& edge) const;
	void reset();

    private:
	vector<index_t> parent;
	vector<value_t> birthtime;
	const CubicalGridComplex_2D* const cgc;
};
