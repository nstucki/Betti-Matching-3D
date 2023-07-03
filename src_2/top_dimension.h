#pragma once

#include "data_structures.h"
#include "config.h"

#include <tuple>
#include <unordered_map>


class TopDimension {
	private:
	const CubicalGridComplex& cgc0;
	const CubicalGridComplex& cgc1;
	const CubicalGridComplex& cgcComp;
	vector<Pair>& pairs0;
	vector<Pair>& pairs1;
	vector<Pair>& pairsComp;
	vector<Match>& matches;
	unordered_map<uint64_t, bool>& isMatched0;
	unordered_map<uint64_t, bool>& isMatched1;
	const Config& config;
	unordered_map<uint64_t, Pair> matchMap0;
	unordered_map<uint64_t, Pair> matchMap1;

	void enumerateDualEdges(const CubicalGridComplex& cgc, vector<Cube>& edges) const;
	void computePairsComp(vector<Cube>& ctr);
	void computePairsImage(vector<Cube>& ctr, uint8_t k);
	void computeMatching();

	public:
	TopDimension(const CubicalGridComplex& cgc0, const CubicalGridComplex& cgc1, const CubicalGridComplex& cgcComp, 
					vector<Pair>& pairs0, vector<Pair>& pairs1, vector<Pair>& pairsComp, vector<Match>& matches,
					unordered_map<uint64_t, bool>& isMatched0, unordered_map<uint64_t, bool>& isMatched1, const Config& config);
	void computePairsAndMatch(vector<Cube>& ctr0, vector<Cube>& ctr1, vector<Cube>& ctrComp);
	
};