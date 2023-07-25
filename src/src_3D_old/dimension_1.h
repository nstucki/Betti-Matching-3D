#pragma once

#include "cubical_grid_complex.h"
#include "write_pair.h"
#include "config.h"
#include <queue>
#include <unordered_map>

typedef priority_queue<Cube, vector<Cube>, CubeComparator> CubeQue;

class Dimension1{
private:
	CubicalGridComplex* cgc;
	vector<WritePair>* pairs;
	Config config;

public:
	Dimension1(CubicalGridComplex* _cgc, vector<WritePair> &_pairs, Config &_config);
	void compute_pairs(vector<Cube> &ctr, bool image=false);
	Cube pop_pivot(CubeQue& column);
	Cube get_pivot(CubeQue& column);
	void add_cache(uint32_t i, CubeQue &working_boundary, unordered_map<uint32_t, CubeQue> &cache);
};