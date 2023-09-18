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
		unordered_map<uint64_t, bool> isPairedComp;
		unordered_map<uint64_t, Pair> matchMap0;
		unordered_map<uint64_t, Pair> matchMap1;
		unordered_map<uint64_t, uint64_t> matchMapIm0;
		unordered_map<uint64_t, uint64_t> matchMapIm1;
		unordered_map<uint64_t, uint64_t> pivotColumnIndex;
#ifdef USE_CACHE
		unordered_map<uint64_t, CubeQueue> cache;
#endif

		void computePairs(const vector<Cube>& ctr, uint8_t k);
		void computePairsComp(vector<Cube>& ctr);
		void computePairsImage(vector<Cube>& ctr, uint8_t k);
		void computeMatching();
		void enumerateEdges(const CubicalGridComplex& cgc, vector<Cube>& edges) const;
		Cube popPivot(CubeQueue& column) const;
		Cube getPivot(CubeQueue& column) const;
#ifdef USE_CACHE
		bool tryCache(const size_t& j, CubeQueue& workingBoundary) const;
		void addCache(const index_t& i, CubeQueue& working_boundary, queue<index_t>& cachedColumnIdx);
#endif
#ifdef USE_APPARENT_PAIRS
		bool isApparentPair(const Cube& face, vector<Cube>& faces, 
								BoundaryEnumerator& enumerator, CoboundaryEnumerator& coEnumerator) const;
		bool isApparentPairImage(const Cube& face, vector<Cube>& faces, 
								BoundaryEnumerator& enumerator, CoboundaryEnumerator& coEnumerator) const;
		bool pivotIsApparentPair(const Cube& pivot, const Cube& column, CoboundaryEnumerator& coEnumerator) const;
		bool pivotIsApparentPairImage(const Cube& pivot, const Cube& column, CoboundaryEnumerator& coEnumerator) const;
#endif
#ifdef USE_EMERGENT_PAIRS
		bool isEmergentPair(const Cube&column, Cube& pivot, vector<Cube>& faces, bool& checkEmergentPair,
								BoundaryEnumerator& enumerator, CoboundaryEnumerator& coEnumerator) const;
		bool isEmergentPairImage(const Cube&column, Cube& pivot, vector<Cube>& faces, bool& checkEmergentPair,
									BoundaryEnumerator& enumerator, CoboundaryEnumerator& coEnumerator, 
									const CubicalGridComplex& cgc) const;
#endif
	};
}
