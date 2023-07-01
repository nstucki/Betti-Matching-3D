#pragma once

#include "data_structures.h"
#include "config.h"
#include <queue>
#include <unordered_map>

typedef priority_queue<Cube, vector<Cube>, CubeComparator> CubeQue;

class InterDimensions{
private:
	const CubicalGridComplex& cgc0;
	const CubicalGridComplex& cgc1;
	const CubicalGridComplex& cgcComp;
	const Config& config;
	uint64_t computeDim;
	unordered_map<uint64_t, uint64_t> pivot_column_index;
	unordered_map<uint64_t, CubeQue > cache;
	unordered_map<uint64_t,Pair> matchMap0;
	unordered_map<uint64_t,Pair> matchMap1;
	unordered_map<uint64_t,bool> matchMapComp;
	unordered_map<uint64_t,Cube> matchMapIm0;
	unordered_map<uint64_t,Cube> matchMapIm1;
	
	Cube popPivot(CubeQue& column) const;
	Cube getPivot(CubeQue& column) const;
	void addCache(uint64_t i, CubeQue &working_boundary, unordered_map<uint64_t, CubeQue> &cache) const;
	void computePairsComp(vector<Cube>& ctr);
	void computePairs(uint8_t k, vector<Cube>& ctr);
	void computePairsImage(uint8_t k, vector<Cube>& ctr);
	void assembleNewColumns(vector<Cube>& ctr);
	void computeMatching();

public:
	vector<vector<Pair>> pairs0;
	vector<vector<Pair>> pairs1;
	vector<vector<Pair>> pairsComp;
	vector<vector<Match>> matches;

	InterDimensions(const CubicalGridComplex& cgc0, const CubicalGridComplex& cgc1, const CubicalGridComplex& cgcComp, const Config& config);
	void computePairsAndMatch(vector<Cube>& ctr0, vector<Cube>& ctr1, vector<Cube>& ctrComp);
};