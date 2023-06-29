#pragma once

#include "data_structures.h"
#include "config.h"

#include <tuple>
#include <unordered_map>


class TopDimension {
	private:
	const Config& config;
	const CubicalGridComplex& cgc0;
	const CubicalGridComplex& cgc1;
	const CubicalGridComplex& cgcComp;

	void enumerateDualEdges(const CubicalGridComplex& cgc, vector<Cube>& edges) const;

	public:
	vector<Pair> pairs0;
	vector<Pair> pairs1;
	vector<Pair> pairsComp;
	vector<Match> matches;
	unordered_map<uint64_t,Pair&> matchMap0;
	unordered_map<uint64_t,Pair&> matchMap1;
	

	TopDimension(const CubicalGridComplex& cgc0, const CubicalGridComplex& cgc1, const CubicalGridComplex& cgcComp, 
					const Config& config);
	void computePairsComp(vector<Cube>& ctr);
	//void computePairsImage(const CubicalGridComplex& cgc, const CubicalGridComplex& cgc_comp, vector<Pair>& pairs, 
	//						unordered_map<uint64_t,Pair&>& match, vector<Cube>& ctr);
	void computePairsImage(uint8_t k, vector<Cube>& ctr);
	void computeMatching();
	
};