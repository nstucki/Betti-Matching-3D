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
	void computePairsAndMatch(vector<Cube>& ctr0, vector<Cube>& ctr1, vector<Cube>& ctrComp, vector<Cube>& ctrImage);
	void computeInput0Pairs(vector<Cube>& ctr0);
	vector<vector<index_t>> getRepresentativeCycle(const Pair& pair, const CubicalGridComplex& cgc);

	private:
	enum ComputePairsMode { INPUT_PAIRS, COMPARISON_PAIRS, IMAGE_PAIRS };

	void computePairs(vector<Cube>& ctr, uint8_t k);
	void computePairsComp(vector<Cube>& ctr);
	void computePairsImage(vector<Cube>& ctr, uint8_t k);

	template <ComputePairsMode computePairsMode>
	void computePairsUnified(vector<Cube>& ctr, uint8_t k);
	void computeMatching();
	void enumerateEdges(vector<Cube>& edges, const CubicalGridComplex& cgc, CubeMap<1, size_t>& pivotColumnIndex) const;
	void enumerateColumnsToReduce(vector<Cube>& ctr, const CubicalGridComplex& cgc) const;
	Cube popPivot(CubeQueue& column) const;
	Cube getPivot(CubeQueue& column) const;
#ifdef USE_REDUCTION_MATRIX
	void useReductionMatrix(const Cube& column, CubeQueue& workingBoundary, BoundaryEnumerator& enumerator
#ifdef USE_CACHE
							, CubeMap<2, vector<Cube>>& cache
#endif
							) const;
#endif
#ifdef USE_CACHE
	bool columnIsCached(const Cube& column, CubeQueue& workingBoundary, CubeMap<2, vector<Cube>>& cache) const;
	void addCache(const Cube& column, CubeQueue& workingBoundary, queue<uint64_t>& cachedColumnIdx, CubeMap<2, vector<Cube>>& cache);
#endif
#if defined(USE_APPARENT_PAIRS) or defined(USE_APPARENT_PAIRS_COMP)
	bool pivotIsApparentPair(const Cube& pivot, vector<Cube>& faces, 
								BoundaryEnumerator& enumerator, CoboundaryEnumerator& coEnumerator) const;
	bool pivotOfColumnIsApparentPair(const Cube& pivot, const Cube& column, vector<Cube>& faces, 
										BoundaryEnumerator& enumerator, CoboundaryEnumerator& coEnumerator) const;
#endif
#ifdef USE_EMERGENT_PAIRS
	template <ComputePairsMode computePairsMode>
	bool isEmergentPair(const Cube&column, Cube& pivot, size_t& j, vector<Cube>& faces, bool& checkEmergentPair,
								const CubicalGridComplex& cgc, BoundaryEnumerator& enumerator, BoundaryEnumerator& enumeratorAP, 
								CoboundaryEnumerator& coEnumeratorAP, CubeMap<1, size_t>& pivotColumnIndex) const;
#endif
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
#ifdef USE_ISPAIRED
	unordered_map<uint64_t, bool> isPairedComp;
#endif
	CubeMap<1, Pair> matchMap0;
	CubeMap<1, Pair> matchMap1;
	CubeMap<1, uint64_t> matchMapIm0;
	CubeMap<1, uint64_t> matchMapIm1;

	CubeMap<1, size_t> pivotColumnIndexInput0; // to be used for input 0 pairs
	CubeMap<1, size_t> pivotColumnIndexInput1; // to be used for input 1 pairs
	CubeMap<1, size_t> pivotColumnIndexComp; // to be used for comparison pairs
	CubeMap<1, size_t> pivotColumnIndexImage0; // to be used for image pairs 0
	CubeMap<1, size_t> pivotColumnIndexImage1; // to be used for image pairs 1
#ifdef USE_REDUCTION_MATRIX
	CubeMap<2, vector<Cube>> reductionMatrix;
#endif
};
}
