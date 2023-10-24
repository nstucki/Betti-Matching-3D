#pragma once

#include "data_structures.h"

#include <unordered_map>
#include <cstdint>


namespace dimN {
class TopDimension {
	public:
	TopDimension(const CubicalGridComplex& cgc0, const CubicalGridComplex& cgc1, const CubicalGridComplex& cgcComp, 
					const Config& config, vector<Pair>& pairs0, vector<Pair>& pairs1, vector<Pair>& pairsComp,
					vector<Match>& matches, unordered_map<index_t, bool>& isMatched0, unordered_map<index_t, bool>& isMatched1);
	void computePairsAndMatch(vector<Cube>& ctr0, vector<Cube>& ctr1, vector<Cube>& ctrComp);

	private:
	const CubicalGridComplex& cgc0;
	const CubicalGridComplex& cgc1;
	const CubicalGridComplex& cgcComp;
	const Config& config;
	vector<Pair>& pairs0;
	vector<Pair>& pairs1;
	vector<Pair>& pairsComp;
	vector<Match>& matches;
	unordered_map<index_t, bool>& isMatched0;
	unordered_map<index_t, bool>& isMatched1;
	unordered_map<index_t, Pair> matchMap0;
	unordered_map<index_t, Pair> matchMap1;
	UnionFindDual uf0;
	UnionFindDual uf1;
	UnionFindDual ufComp;

	void enumerateDualEdges(const CubicalGridComplex& cgc, vector<Cube>& dualEdges) const;
	void computeImagePairs(vector<Cube>& dualEdges, uint8_t k);
	void computePairsCompAndMatch(vector<Cube>& dualEdges);
};
}