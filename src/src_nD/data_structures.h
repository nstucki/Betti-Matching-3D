#pragma once

#include "../config.h"

#include <vector>

using namespace std;


namespace dimN {
class Cube {
	public:
	value_t birth;
	vector<index_t> coordinates;

	Cube();
	Cube(value_t birth, vector<index_t> coordinates);
	Cube(const Cube& cube);
	bool operator==(const Cube& rhs) const;
	bool operator!=(const Cube& rhs) const;
	void print() const;
};


struct CubeComparator {
	bool operator()(const Cube& Cube1, const Cube& Cube2) const;
};


class Pair {
	public:
	Cube birth;
	Cube death;

	Pair();
	Pair(const Cube& birth, const Cube& death);
	Pair(const Pair& pair);
	bool operator==(const Pair& rhs) const;
	void print() const;
};


class Match{
	public:
	const Pair pair0;
	const Pair pair1;

	Match(const Pair& pair_0, const Pair& pair_1);
	void print() const;
};


class CubicalGridComplex {
	public:
	const vector<index_t> shape;
	const index_t dim;

	CubicalGridComplex(const vector<value_t> _image, const vector<index_t>& _shape);
	index_t getNumberOfCubes(const index_t dim) const;
	index_t getCubeIndex(const Cube& cube) const;
	index_t getCubeIndex(const vector<index_t>& coordinates) const;
	vector<index_t> getCubeCoordinates(index_t idx) const;
	value_t getBirth (const vector<index_t>& coordinates) const;
	vector<index_t> getParentVoxel(const Cube& c) const;
	void printImage() const;
	void printCubes() const;

	private:
	const vector<value_t> image;
	vector<index_t> pixelCoordFactor;
	vector<index_t> cubeCoordFactor;

	index_t getIndex(const vector<index_t>& pixelCoordinates) const;
	value_t getValue(const vector<index_t>& pixelCoordinates) const;
	bool getParentVoxelRecursion(const value_t& birth, const size_t& axis, 
									vector<index_t>& parentVoxel, vector<size_t>& nonDegen) const;
};


class UnionFind{
	public:
	UnionFind(const CubicalGridComplex& cgc);
	index_t find(index_t x);
	index_t link(const index_t& x, const index_t& y);
	value_t getBirth(const index_t& idx) const;
	index_t getIndex(const vector<index_t>& coordinates) const;
	vector<index_t> getCoordinates(index_t idx) const;
	vector<index_t> getBoundaryIndices(const Cube& edge) const;
	void reset();

	private:
	vector<index_t> parent;
	vector<value_t> birthtime;
	const CubicalGridComplex& cgc;
};


class UnionFindDual {
	public:
	UnionFindDual(const CubicalGridComplex& cgc);
	index_t find(index_t x);
	index_t link(const index_t& x, const index_t& y);
	value_t getBirth(const index_t& idx) const;
	index_t getIndex(const vector<index_t>& coordinates) const;
	vector<index_t> getCoordinates(index_t idx) const;
	vector<index_t> getBoundaryIndices(const Cube& edge) const;
	void reset();

	private:
	vector<index_t> parent;
	vector<value_t> birthtime;
	const CubicalGridComplex& cgc;
	index_t star;
};
}