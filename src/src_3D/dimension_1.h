#pragma once

#include "data_structures.h"
#include "config.h"

#include <queue>
#include <unordered_map>

typedef priority_queue<Cube, vector<Cube>, CubeComparator> CubeQue;

class Dimension1 {
	public:
	Dimension1(const CubicalGridComplex* const cgc0, const CubicalGridComplex* const cgc1, const CubicalGridComplex* const cgcComp, 
				const Config& config, vector<Pair>& pairs0, vector<Pair>& pairs1, vector<Pair>& pairsComp,
				vector<Match>& matches, unordered_map<index_t, bool>& isMatched0, unordered_map<index_t, bool>& isMatched1);

	private:
	const CubicalGridComplex* const cgc0;
	const CubicalGridComplex* const cgc1;
	const CubicalGridComplex* const cgcComp;
	const Config& config;
	vector<Pair>& pairs0;
	vector<Pair>& pairs1;
	vector<Pair>& pairsComp;
	vector<Match>& matches;
	unordered_map<uint64_t, bool>& isMatched0;
	unordered_map<uint64_t, bool>& isMatched1;
	unordered_map<uint64_t, bool> matchMapComp;
	unordered_map<uint64_t, Pair*> matchMap0;
	unordered_map<uint64_t, Pair*> matchMap1;
	unordered_map<uint64_t, Cube*> matchMapIm0;
	unordered_map<uint64_t, Cube*> matchMapIm1;
	unordered_map<uint64_t, uint64_t> pivotColumnIndex;
	unordered_map<uint64_t, CubeQue> cache;

	void computePairsComp(vector<Cube>& ctr);
	void computePairs(const vector<Cube>& ctr, uint8_t k);
	void computeImagePairs(const vector<Cube>& ctr, uint8_t k);
	void computeMatching();
	void assembleColumnsToReduce(const CubicalGridComplex& cgc, vector<Cube>& ctr) const;
	Cube popPivot(CubeQue& column) const;
	Cube getPivot(CubeQue& column) const;
	void addCache(index_t i, CubeQue &working_boundary);
};