#pragma once

#include "data_structures.h"
#include "enumerators.h"

#include <queue>
#include <unordered_map>


namespace dim3 {
	typedef priority_queue<Cube, vector<Cube>, CubeComparator> CubeQueue;
	

	class Dimension1 {
		public:
		Dimension1(const CubicalGridComplex& cgc0, const CubicalGridComplex& cgc1, const CubicalGridComplex& cgcComp, 
					const Config& config, vector<Pair>& pairs0, vector<Pair>& pairs1, vector<Pair>& pairsComp,
					vector<Match>& matches, unordered_map<uint64_t, bool>& isMatched0, unordered_map<uint64_t, bool>& isMatched1);
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
		unordered_map<uint64_t, bool>& isMatched0;
		unordered_map<uint64_t, bool>& isMatched1;
		unordered_map<uint64_t, bool> isMatchedComp;
		unordered_map<uint64_t, Pair> matchMap0;
		unordered_map<uint64_t, Pair> matchMap1;
		unordered_map<uint64_t, uint64_t> matchMapIm0;
		unordered_map<uint64_t, uint64_t> matchMapIm1;
		unordered_map<uint64_t, uint64_t> pivotColumnIndex;
		unordered_map<uint64_t, CubeQueue> cache;

		void computePairs(const vector<Cube>& ctr, uint8_t k);
		void computePairsComp(vector<Cube>& ctr);
		void computeImagePairs(vector<Cube>& ctr, uint8_t k);
		void computeMatching();
		void enumerateEdges(const CubicalGridComplex& cgc, vector<Cube>& edges) const;
		Cube popPivot(CubeQueue& column) const;
		Cube getPivot(CubeQueue& column) const;
		void addCache(index_t i, CubeQueue &working_boundary);
	};
}
