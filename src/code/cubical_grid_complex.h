#pragma once

#include "cube.h"
#include <vector>

using namespace std;

class CubicalGridComplex{
public:
	uint8_t dim;
	vector<uint32_t> shape;
	uint32_t n_x;
	uint32_t n_y;
	uint32_t n_z;
	//testen ob sinnvoll
	uint32_t n_yz;
	uint32_t n_xyz;
	uint32_t m_yz;
	uint32_t m_xyz;
	// vector
	double*** grid;

	CubicalGridComplex(const vector<double>& image, const vector<uint32_t> shape);
	~CubicalGridComplex();
	double*** allocate_memory();
	void gridFromVector(const vector<double> vec);
	double getBirth(uint32_t x, uint32_t y, uint32_t z);
	double getBirth(uint32_t x, uint32_t y, uint32_t z, uint8_t t, uint8_t dim);
	vector<uint32_t> ParentVoxel(Cube &c);
};