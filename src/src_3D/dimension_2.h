#pragma once

#include "../config.h"
#include "data_structures.h"
#include "enumerators.h"

#include <unordered_map>


namespace dim3 {
	class Dimension2 {
		public:
		Dimension2(const CubicalGridComplex& cgc0, const CubicalGridComplex& cgc1, const CubicalGridComplex& cgcComp, 
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
		unordered_map<index_t, Pair> matchMap0;
		unordered_map<index_t, Pair> matchMap1;
		UnionFindDual uf0;
		UnionFindDual uf1;
		UnionFindDual ufComp;

		void enumerateDualEdges(const CubicalGridComplex& cgc, vector<Cube>& dualEdges) const;
		void computeImagePairs(vector<Cube>& dualEdges, const uint8_t& k);
		void computeCompPairsAndMatch(vector<Cube>& dualEdges);
#ifdef USE_APPARENT_PAIRS
		bool isApparentPair(BoundaryEnumerator& enumerator, CoboundaryEnumerator& coEnumerator, const Cube& dualEdge) const;
#endif
	};
}