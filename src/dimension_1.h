#pragma once

#include "cubical_grid_complex.h"
#include "write_pair.h"
#include <queue>

typedef priority_queue<Cube, vector<Cube>, CubeComparator> CubeQue;

class Dimension1{
private:
	CubicalGridComplex* cgc;
	vector<WritePair>* pairs;

public:
	Dimension1(CubicalGridComplex* _cgc, vector<WritePair>& _pairs);
	void compute_pairs(vector<Cube>& ctr, bool image=false);
	Cube pop_pivot(CubeQue& column);
	Cube get_pivot(CubeQue& column);
};