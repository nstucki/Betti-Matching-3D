#pragma once

#include "cubical_grid_complex.h"

class UnionFind{
private:
	vector<uint32_t> parent;
	vector<double> time_max;
public:
	vector<double> birthtime;

	UnionFind(CubicalGridComplex* const _cgc);
	uint32_t find(uint32_t x);
	uint32_t link(uint32_t x, uint32_t y);
};

class UnionFindDual{
private:
	vector<uint32_t> parent;
	vector<double> time_max;
public:
	vector<double> birthtime;

	UnionFindDual(CubicalGridComplex* const _cgc);
	uint32_t find(uint32_t x);
	uint32_t link(uint32_t x, uint32_t y);
};
