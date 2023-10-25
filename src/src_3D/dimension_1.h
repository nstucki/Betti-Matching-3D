#pragma once

#include "data_structures.h"
#include "enumerators.h"

#include <cstdint>
#include <optional>
#include <queue>
#include <unordered_map>
#include <vector>


namespace dim3 {
	typedef priority_queue<Cube, vector<Cube>, CubeComparator> CubeQueue;


	class Dimension1 {
		public:
		Dimension1(const CubicalGridComplex& cgc0, const CubicalGridComplex& cgc1, const CubicalGridComplex& cgcComp, 
					const Config& config, vector<Pair>& pairs0, vector<Pair>& pairs1, vector<Pair>& pairsComp,
					vector<Match>& matches,
					Cube1Map<bool>& _isMatched0, Cube1Map<bool>& _isMatched1
					);
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
		Cube1Map<bool> &isMatched0;
		Cube1Map<bool> &isMatched1;
		Cube1Map<bool> isMatchedComp;

		Cube1Map<Pair> matchMap0;
		Cube1Map<Pair> matchMap1;
		Cube1Map<uint64_t> matchMapIm0;
		Cube1Map<uint64_t> matchMapIm1;
		Cube1Map<uint64_t> pivotColumnIndex0;
		Cube1Map<uint64_t> pivotColumnIndex1;
		Cube1Map<uint64_t> pivotColumnIndexComp;

		void computePairs(const vector<Cube>& ctr, uint8_t k);
		void computePairsComp(vector<Cube>& ctr);
		void computeImagePairs(vector<Cube>& ctr, uint8_t k);
		void computeMatching();
		void enumerateEdges(const CubicalGridComplex& cgc, vector<Cube>& edges, Cube1Map<uint64_t> &pivotColumnIndex) const;
		Cube popPivot(CubeQueue& column) const;
		Cube getPivot(CubeQueue& column) const;
		void addCache(const index_t& i, CubeQueue& working_boundary, vector<std::optional<vector<Cube>>> &cache);
	};
}
