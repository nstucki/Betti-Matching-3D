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
	const Config& config;
	unordered_map<uint64_t,Pair> matchMap0;
	unordered_map<uint64_t,Pair> matchMap1;

	void enumerateDualEdges(const CubicalGridComplex& cgc, vector<Cube>& edges) const;
	void computePairsComp(vector<Cube>& ctr);
	void computePairsImage(uint8_t k, vector<Cube>& ctr);
	void computeMatching();

	public:
	vector<Pair> pairs0;
	vector<Pair> pairs1;
	vector<Pair> pairsComp;
	vector<Match> matches;
	
	TopDimension(const CubicalGridComplex& cgc0, const CubicalGridComplex& cgc1, const CubicalGridComplex& cgcComp, 
					const Config& config);
	void computePairsAndMatch(vector<Cube>& ctr0, vector<Cube>& ctr1, vector<Cube>& ctrComp);
	
};