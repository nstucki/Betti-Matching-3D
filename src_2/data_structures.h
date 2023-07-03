#pragma once

#include <vector>

using namespace std;

#define NONE numeric_limits<uint64_t>::max()


class Cube {
	public:
	float birth;
	vector<uint64_t> coordinates;

	Cube();
    Cube(float birth, vector<uint64_t> coordinates);
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
	private:
	const vector<float> image;

	uint64_t getIndex(const vector<uint64_t>& coordinates) const;
	float getValue(const vector<uint64_t>& coordinates) const;

	public:
	const vector<uint64_t> shape;
	const uint64_t dim;
	
	CubicalGridComplex(const vector<float> _image, const vector<uint64_t>& _shape);
	uint64_t getNumberOfCubes(const uint64_t dim) const;
	uint64_t getCubeIndex(const Cube& cube) const;
	uint64_t getCubeIndex(const vector<uint64_t>& coordinates) const;
	float getBirth (const vector<uint64_t>& coordinates) const;
	void printImage() const;
	void printCubes() const;
};


class UnionFind{
	private:
	vector<uint64_t> parent;
	vector<float> birthtime;
	const CubicalGridComplex& cgc;

	public:
	UnionFind(const CubicalGridComplex& cgc);
	uint64_t find(uint64_t x);
	uint64_t link(const uint64_t& x, const uint64_t& y);
	float getBirth(const uint64_t& idx) const;
	vector<uint64_t> getCoordinates(uint64_t idx) const;
};


class UnionFindDual {
	private:
	const CubicalGridComplex& cgc;
	vector<uint64_t> parent;
	vector<float> birthtime;

	public:
	uint64_t star;
	
	UnionFindDual(const CubicalGridComplex& cgc);
	uint64_t find(uint64_t x);
	uint64_t link(const uint64_t& x, const uint64_t& y);
	float getBirth(const uint64_t& idx) const;
	vector<uint64_t> getCoordinates(uint64_t idx) const;
	uint64_t getIndex(const vector<uint64_t>& coordinates) const;
	
};