#pragma once

#include "data_structures.h"

#include <unordered_map>


namespace dim3 {
class Dimension0 {
	public:
	Dimension0(const CubicalGridComplex& cgc0, const CubicalGridComplex& cgc1, const CubicalGridComplex& cgcComp, 
				const Config& config, vector<Pair>& pairs0, vector<Pair>& pairs1, vector<Pair>& pairsComp, 
				vector<Match>& matches, unordered_map<uint64_t, bool>& isMatched0, unordered_map<uint64_t, bool>& isMatched1);
	void computePairsAndMatch(vector<Cube>& ctr0, vector<Cube>& ctr1, vector<Cube>& ctrComp);
	void computeInput0Pairs(vector<Cube>& ctr0);
	vector<vector<index_t>> getRepresentativeCycle(const Pair& pair, const CubicalGridComplex& cgc) const;

	private:
	void computePairs(vector<Cube>& ctr, uint8_t k);
	void computeImagePairsAndMatch(vector<Cube>& ctr);
	void enumerateEdges(vector<Cube>& edges, const CubicalGridComplex& cgc) const;
	const CubicalGridComplex& cgc0;
	const CubicalGridComplex& cgc1;
	const CubicalGridComplex& cgcComp;
	const Config& config;
	vector<Pair>& pairs0;
	vector<Pair>& pairs1;
	vector<Pair>& pairsComp;
	vector<Match>& matches;
	unordered_map<uint64_t, bool>& isMatched0;
	unordered_map<uint64_t, bool>& isMatched1;
	unordered_map<index_t, Pair> matchMap0;
	unordered_map<index_t, Pair> matchMap1;
	UnionFind uf0;
	UnionFind uf1;
	UnionFind ufComp;
};
}
