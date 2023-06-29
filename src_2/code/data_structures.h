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
    void copyCube(const Cube& cube);
    bool operator==(const Cube& rhs) const;
    void print() const;
};


struct CubeComparator {
    bool operator()(const Cube& Cube1, const Cube& Cube2) const;
};


class Pair {
	public:
	const Cube birth;
	const Cube death;

    Pair(const Cube& birth, const Cube& death);
	bool operator==(const Pair& rhs) const;
	void print() const;
};


class Match{
	public:
	const Pair& pair0;
	const Pair& pair1;

    Match(const Pair& pair_0, const Pair& pair_1);
	void print() const;
};


class CubicalGridComplex {
private:
	const vector<float> image;

public:
	const vector<uint64_t> shape;
	const uint64_t dim;
	
	CubicalGridComplex(const vector<float>& _image, const vector<uint64_t>& _shape);
	uint64_t getIndex(const vector<uint64_t>& coordinates) const;
	uint64_t getCubeIndex(const Cube& c) const;
	uint64_t getCubeIndex(const vector<uint64_t> coordinates) const;
	float getValue(const vector<uint64_t>& coordinates) const;
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
	float getBirth(uint64_t x) const;
	vector<uint64_t> getCoordinates(uint64_t x) const;
	uint64_t find(uint64_t x);
	uint64_t link(uint64_t x, uint64_t y);
};


class UnionFindDual {
private:
	vector<uint64_t> parent;
	vector<float> birthtime;
	const CubicalGridComplex& cgc;

public:
	uint64_t n;
	
	UnionFindDual(const CubicalGridComplex& cgc);
	float getBirth(uint64_t x) const;
	vector<uint64_t> getCoordinates(uint64_t x) const;
	uint64_t getIndex(vector<uint64_t> coordinates) const;
	uint64_t find(uint64_t x);
	uint64_t link(uint64_t x, uint64_t y);
};