#include "top_dimension.h"
#include "enumerators.h"

#include <iostream>
#include <chrono>
#include <limits>
#include <algorithm>
#include <cstdint>

using namespace dimN;
using namespace std;
using namespace std::chrono;


TopDimension::TopDimension(const CubicalGridComplex& _cgc0, const CubicalGridComplex& _cgc1, const CubicalGridComplex& _cgcComp,
					const Config& _config, vector<Pair>& _pairs0, vector<Pair>& _pairs1, vector<Pair>& _pairsComp,
					vector<Match>& _matches, unordered_map<index_t, bool>& _isMatched0, unordered_map<index_t, bool>& _isMatched1) : 
					cgc0(_cgc0), cgc1(_cgc1), cgcComp(_cgcComp), config(_config), 
					pairs0(_pairs0), pairs1(_pairs1), pairsComp(_pairsComp), 
					matches(_matches), isMatched0(_isMatched0), isMatched1(_isMatched1),
					uf0(UnionFindDual(cgc0)), uf1(UnionFindDual(cgc1)), ufComp(UnionFindDual(cgcComp)) {}

void TopDimension::computePairsAndMatch(vector<Cube>& ctr0, vector<Cube>& ctr1, vector<Cube>& ctrComp) {
	size_t actualDim = cgc0.dim;
	for (index_t s : cgc0.shape) { if (s == 1) { --actualDim; } }
	bool needToCompute = (cgc0.dim == actualDim);
	if (!needToCompute) {
#ifdef RUNTIME
		cout << endl << "input & image 0: ";
#endif
		enumerateDualEdges(cgc0, ctr0);
#ifdef RUNTIME
		cout << "barcodes not computed";
#endif

#ifdef RUNTIME
	cout << endl << "input & image 1: ";
#endif
		enumerateDualEdges(cgc1, ctr1);
#ifdef RUNTIME
		cout << "barcodes not computed";
#endif

#ifdef RUNTIME
	cout << endl << "comparison & matching: ";
#endif
		enumerateDualEdges(cgcComp, ctrComp);
#ifdef RUNTIME
		cout << "barcode & matching not computed";
#endif
	} else {
		#ifdef RUNTIME
	cout << endl << "input & image 0: ";
#endif
	enumerateDualEdges(cgc0, ctr0);
	computeImagePairs(ctr0, 0);

#ifdef RUNTIME
	cout << endl << "input & image 1: ";
#endif
	enumerateDualEdges(cgc1, ctr1);
	ufComp.reset();
	computeImagePairs(ctr1, 1);

#ifdef RUNTIME
	cout << endl << "comparison & matching: ";
#endif
	enumerateDualEdges(cgcComp, ctrComp);
	ufComp.reset();
	computePairsCompAndMatch(ctrComp);
	}
}

void TopDimension::enumerateDualEdges(const CubicalGridComplex& cgc, vector<Cube>& dualEdges) const {
#ifdef RUNTIME
	cout << "enumeration ";
	auto start = high_resolution_clock::now();
#endif
	dualEdges.reserve(cgc.getNumberOfCubes(cgc.dim-1));
	CubeEnumerator cubeEnum(cgc, cgc.dim-1);
	Cube cube = cubeEnum.getNextCube();
	if (cube.birth < config.threshold) { dualEdges.push_back(cube); }
	while (cubeEnum.hasNextCube()) {
		cube = cubeEnum.getNextCube();
		if (cube.birth < config.threshold) { dualEdges.push_back(cube); }
	}
	sort(dualEdges.begin(), dualEdges.end(), CubeComparator());
#ifdef RUNTIME
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<milliseconds>(stop - start);
	cout << duration.count() << " ms, ";
#endif
}

void TopDimension::computeImagePairs(vector<Cube>& dualEdges, uint8_t k) {
#ifdef RUNTIME
cout << "barcodes ";
auto start = high_resolution_clock::now();
#endif
	const CubicalGridComplex& cgc = (k == 0) ? cgc0 : cgc1;
	UnionFindDual& uf = (k == 0) ? uf0 : uf1; 
	vector<Pair>& pairs = (k == 0) ? pairs0 : pairs1;
	unordered_map<index_t, Pair>& matchMap = (k==0) ? matchMap0 : matchMap1;
	
	vector<index_t> boundaryIndices;
	index_t parentIdx0;
	index_t parentIdx1;
	index_t birthIdx;
	index_t birthIdxComp;
	value_t birth;
	for (auto edge = dualEdges.rbegin(), last = dualEdges.rend(); edge != last; ++edge) {
		boundaryIndices = uf.getBoundaryIndices(*edge);
		parentIdx0 = uf.find(boundaryIndices[0]);
		parentIdx1 = uf.find(boundaryIndices[1]);
		if (parentIdx0 != parentIdx1) {
			birthIdx = uf.link(parentIdx0, parentIdx1);
			birth = uf.getBirth(birthIdx);
			parentIdx0 = ufComp.find(boundaryIndices[0]);
			parentIdx1 = ufComp.find(boundaryIndices[1]);
			birthIdxComp = ufComp.link(parentIdx0, parentIdx1);
			if (edge->birth != birth) {
				pairs.push_back(Pair(*edge, Cube(birth, uf.getCoordinates(birthIdx))));
				matchMap.emplace(birthIdxComp, pairs.back());
			} 
			edge->coordinates[0] = NONE;
		}
	}
	auto new_end = remove_if(dualEdges.begin(), dualEdges.end(), [](const Cube& edge){ return edge.coordinates[0] == NONE; });
	dualEdges.erase(new_end, dualEdges.end());
#ifdef RUNTIME
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<milliseconds>(stop - start);
	cout << duration.count() << " ms";
#endif
}

void TopDimension::computePairsCompAndMatch(vector<Cube>& dualEdges) {
#ifdef RUNTIME
	cout << "barcode and matching ";
	auto start = high_resolution_clock::now();
#endif 
	vector<index_t> boundaryIndices;
	index_t parentIdx0;
	index_t parentIdx1;
	index_t birthIdx;
	value_t birth;
	for (auto edge = dualEdges.rbegin(), last = dualEdges.rend(); edge != last; ++edge) {
		boundaryIndices = ufComp.getBoundaryIndices(*edge);
		parentIdx0 = ufComp.find(boundaryIndices[0]);
		parentIdx1 = ufComp.find(boundaryIndices[1]);
		if (parentIdx0 != parentIdx1) {
			birthIdx = ufComp.link(parentIdx0, parentIdx1);
			birth = ufComp.getBirth(birthIdx);
			if (edge->birth != birth) { 
#ifdef COMPUTE_COMPARISON
				pairsComp.push_back(Pair(*edge, Cube(birth, ufComp.getCoordinates(birthIdx))));
#endif
				auto find0 = matchMap0.find(birthIdx);
				auto find1 = matchMap1.find(birthIdx);
				if (find0 != matchMap0.end() && find1 != matchMap1.end()) {
					matches.push_back(Match(find0->second, find1->second));
					isMatched0.emplace(cgc0.getCubeIndex(find0->second.birth), true);
					isMatched1.emplace(cgc1.getCubeIndex(find1->second.birth), true);
				}
			}
			edge->coordinates[0] = NONE;
		}
	}
	auto new_end = remove_if(dualEdges.begin(), dualEdges.end(), [](const Cube& edge){ return edge.coordinates[0] == NONE; });
	dualEdges.erase(new_end, dualEdges.end());
#ifdef RUNTIME
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<milliseconds>(stop - start);
	cout << duration.count() << " ms";
#endif
}