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
	vector<vector<Pair>>& pairs0;
	vector<vector<Pair>>& pairs1;
	vector<vector<Pair>>& pairsComp;
	vector<vector<Match>>& matches;
	uint64_t computeDim;
	unordered_map<uint64_t, uint64_t> pivot_column_index;
	unordered_map<uint64_t, CubeQue > cache;
	unordered_map<uint64_t, Pair> matchMap0;
	unordered_map<uint64_t, Pair> matchMap1;
	unordered_map<uint64_t, bool> matchMapComp;
	unordered_map<uint64_t, Cube> matchMapIm0;
	unordered_map<uint64_t, Cube> matchMapIm1;
	
	Cube popPivot(CubeQue& column) const;
	Cube getPivot(CubeQue& column) const;
	void addCache(uint64_t i, CubeQue &working_boundary);
	void computePairsComp(const vector<Cube>& ctr);
	void computePairs(const vector<Cube>& ctr, uint8_t k);
	void computePairsImage(const vector<Cube>& ctr, uint8_t k);
	void assembleColumnsToReduce(const CubicalGridComplex& cgc, vector<Cube>& ctr) const;
	void computeMatching();

public:

	InterDimensions(const CubicalGridComplex& cgc0, const CubicalGridComplex& cgc1, const CubicalGridComplex& cgcComp, 
					const Config& config, vector<vector<Pair>>& pairs0, vector<vector<Pair>>& pairs1, vector<vector<Pair>>& pairsComp,
					vector<vector<Match>>& matches);
	void computePairsAndMatch(vector<Cube>& ctr0, vector<Cube>& ctr1, vector<Cube>& ctrComp);
};