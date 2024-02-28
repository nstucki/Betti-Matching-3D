#include "dimension_1.h"

#include <iostream>
#include <chrono>
#include <algorithm>
#include <set>

using namespace dim2;
using namespace std;
using namespace std::chrono;



Dimension1::Dimension1(const CubicalGridComplex& _cgc0, const CubicalGridComplex& _cgc1, 
						const CubicalGridComplex& _cgcComp, const Config& _config, 
						vector<Pair>& _pairs0, vector<Pair>& _pairs1, vector<Pair>& _pairsComp, vector<Match>& _matches, 
						unordered_map<uint64_t, bool>& _isMatched0, unordered_map<uint64_t, bool>& _isMatched1) :
						cgc0(_cgc0), cgc1(_cgc1), cgcComp(_cgcComp), config(_config), 
						pairs0(_pairs0), pairs1(_pairs1), pairsComp(_pairsComp),
						matches(_matches), isMatched0(_isMatched0), isMatched1(_isMatched1), 
						uf0(UnionFindDual(cgc0)), uf1(UnionFindDual(cgc1)), ufComp(UnionFindDual(cgcComp)) {}


void Dimension1::computePairsAndMatch(vector<Cube>& ctr0, vector<Cube>& ctr1, vector<Cube>& ctrComp) {
#ifdef RUNTIME
	cout << endl << "input & image 0: ";
#endif
	enumerateDualEdges(ctr0, cgc0);
    computeImagePairs(ctr0, 0);

#ifdef RUNTIME
	cout << endl << "input & image 1: ";
#endif
	enumerateDualEdges(ctr1, cgc1);
	ufComp.reset();
    computeImagePairs(ctr1, 1);
    
#ifdef RUNTIME
	cout << endl << "comparison & matching: ";
#endif
	enumerateDualEdges(ctrComp, cgcComp);
	ufComp.reset();
	computeCompPairsAndMatch(ctrComp);
}

void Dimension1::computeInput0Pairs(vector<Cube>& ctr0) {
	enumerateDualEdges(ctr0, cgc0);
    computeImagePairs(ctr0, 0);
}

vector<vector<index_t>> Dimension1::getRepresentativeCycle(const Pair& pair, const CubicalGridComplex& cgc) const {
	vector<Cube> dualEdges;
	enumerateDualEdges(dualEdges, cgc);

	UnionFindDual uf(cgc);
	vector<index_t> boundaryIndices(2);
	index_t parentIdx0;
	index_t parentIdx1;
	index_t birthIdx;
	for (auto edge = dualEdges.rbegin(), last = dualEdges.rend(); edge != last; ++edge) {
		if (*edge == pair.birth) { break; }
		boundaryIndices = uf.getBoundaryIndices(*edge);
		parentIdx0 = uf.find(boundaryIndices[0]);
		parentIdx1 = uf.find(boundaryIndices[1]);
		if (parentIdx0 != parentIdx1) { birthIdx = uf.link(parentIdx0, parentIdx1); }
	}

	set<vector<index_t>> cubeCoordinates;
	parentIdx0 = uf.find(pair.death.x()*cgc.m_y + pair.death.y());
	for (size_t i = 0; i < cgc.getNumberOfCubes(2); ++i) {
		parentIdx1 = uf.find(i);
		if (parentIdx0 == parentIdx1) { cubeCoordinates.insert(uf.getCoordinates(i)); }
	}

	multiset<vector<index_t>> boundaryVertices;
	vector<index_t> vertex;
	for (const vector<index_t>& c : cubeCoordinates) {
		for (uint8_t x = 0; x < 2; ++x) {
			for (uint8_t y = 0; y < 2; ++y) {
				vertex = c;
				vertex[0] += x;
				vertex[1] += y;
				boundaryVertices.insert(vertex);
			}
		}
	}

	vector<vector<index_t>> reprCycle;
	reprCycle.push_back(cgc.getParentVoxel(pair.birth, 1));
	for (const vector<index_t>& vertex : boundaryVertices) {
		auto lower = boundaryVertices.lower_bound(vertex);
		auto upper = boundaryVertices.upper_bound(vertex);
		int multiplicity = distance(lower, upper);
		if (multiplicity < 4) { if(find(reprCycle.begin(), reprCycle.end(), vertex) == reprCycle.end()) { reprCycle.push_back(vertex); } }
	}
	reprCycle.push_back(cgc.getParentVoxel(pair.death, 2));

	return reprCycle;
}


void Dimension1::enumerateDualEdges(vector<Cube>& dualEdges, const CubicalGridComplex& cgc) const {
#ifdef RUNTIME
	cout << "enumeration ";
	auto start = high_resolution_clock::now();
#endif

	dualEdges.reserve(cgc.getNumberOfCubes(1));
	value_t birth;
	for (index_t x = 0; x < cgc.shape[0]; ++x) {
		for (index_t y = 0; y < cgc.shape[1]; ++y) {
			for (uint8_t type = 0; type < 2; ++type) {
				birth = cgc.getBirth(x, y, type, 1);
				if (birth < config.threshold) { dualEdges.push_back(Cube(birth, x, y, type)); }
			}
		}
	}

	sort(dualEdges.begin(), dualEdges.end(), CubeComparator());

#ifdef RUNTIME
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<milliseconds>(stop - start);
	cout << duration.count() << " ms, ";
#endif
}


void Dimension1::computeImagePairs(vector<Cube>& dualEdges, const uint8_t& k) {
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
	vector<index_t> birthCoordinates(2);
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
				birthCoordinates = uf.getCoordinates(birthIdx);
				pairs.push_back(Pair(*edge, Cube(birth, birthCoordinates[0], birthCoordinates[1], 0)));
				matchMap.emplace(birthIdxComp, pairs.back());
			}
#ifdef USE_CLEARING_DIM0
			edge->index = NONE_INDEX;
#endif
		}
	}

#ifdef USE_CLEARING_DIM0
	auto new_end = remove_if(dualEdges.begin(), dualEdges.end(), [](const Cube& cube){ return cube.index == NONE_INDEX; });
	dualEdges.erase(new_end, dualEdges.end());
#endif

#ifdef RUNTIME
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<milliseconds>(stop - start);
	cout << duration.count() << " ms";
#endif
}


void Dimension1::computeCompPairsAndMatch(vector<Cube>& dualEdges) {
#ifdef RUNTIME
	cout << "barcode and matching ";
	auto start = high_resolution_clock::now();
#endif 

	vector<index_t> boundaryIndices;
	index_t parentIdx0;
	index_t parentIdx1;
	index_t birthIdx;
	value_t birth;
	vector<index_t> birthCoordinates(2);
	for (auto edge = dualEdges.rbegin(), last = dualEdges.rend(); edge != last; ++edge){
		boundaryIndices = ufComp.getBoundaryIndices(*edge);
		parentIdx0 = ufComp.find(boundaryIndices[0]);
		parentIdx1 = ufComp.find(boundaryIndices[1]);
		if (parentIdx0 != parentIdx1) {
			birthIdx = ufComp.link(parentIdx0, parentIdx1);
			birth = ufComp.getBirth(birthIdx);
			if (edge->birth != birth) {
#ifdef COMPUTE_COMPARISON
				birthCoordinates = ufComp.getCoordinates(birthIdx);
				pairsComp.push_back(Pair(*edge, Cube(birth, birthCoordinates[0], birthCoordinates[1], 0)));
#endif
				auto find0 = matchMap0.find(birthIdx);
				auto find1 = matchMap1.find(birthIdx);
				if (find0 != matchMap0.end() && find1 != matchMap1.end()) {
					matches.push_back(Match(find0->second, find1->second));
					isMatched0.emplace(find0->second.birth.index, true);
					isMatched1.emplace(find1->second.birth.index, true);
				}
			}
#ifdef USE_CLEARING_DIM0
			edge->index = NONE_INDEX;
#endif
		}
	}

#ifdef USE_CLEARING_DIM0
	auto new_end = std::remove_if(dualEdges.begin(), dualEdges.end(), [](const Cube& cube){ return cube.index == NONE_INDEX; });
	dualEdges.erase(new_end, dualEdges.end());
#endif

#ifdef RUNTIME
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<milliseconds>(stop - start);
	cout << duration.count() << " ms";
#endif
}