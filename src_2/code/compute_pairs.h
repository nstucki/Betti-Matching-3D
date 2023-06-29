#pragma once

#include "data_structures.h"
#include "config.h"
#include <queue>
#include <unordered_map>

typedef priority_queue<Cube, vector<Cube>, CubeComparator> CubeQue;

class ComputePairs{
private:
	Config& config;

public:
	vector<Pair> pairs0;
	vector<Pair> pairs1;
	vector<Pair> pairsComp;
	vector<tuple<Pair*, Pair*>> match;
	vector<tuple<Cube,Pair*>> match0;
	vector<tuple<Cube,Pair*>> match1;

	ComputePairs(Config& config);
	void computePairs(const CubicalGridComplex& cgc, vector<Cube>& ctr, bool image);
	Cube popPivot(CubeQue& column) const;
	Cube getPivot(CubeQue& column) const;
	void addCache(uint64_t i, CubeQue &working_boundary, unordered_map<uint64_t, CubeQue> &cache) const;
};