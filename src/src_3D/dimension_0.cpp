#include "dimension_0.h"

#include <iostream>
#include <chrono>
#include <algorithm>

using namespace std::chrono;


Dimension0::Dimension0(const CubicalGridComplex* const _cgc0, const CubicalGridComplex* const _cgc1, 
						const CubicalGridComplex* const _cgcComp, const Config& _config, 
						vector<Pair>& _pairs0, vector<Pair>& _pairs1, vector<Pair>& _pairsComp, vector<Match>& _matches, 
						unordered_map<uint64_t, bool>& _isMatched0, unordered_map<uint64_t, bool>& _isMatched1) : 
						cgc0(_cgc0), cgc1(_cgc1), cgcComp(_cgcComp), config(_config), 
						pairs0(_pairs0), pairs1(_pairs1), pairsComp(_pairsComp),
						matches(_matches), isMatched0(_isMatched0), isMatched1(_isMatched1),
						uf0(UnionFind(cgc0)), uf1(UnionFind(cgc1)), ufComp(UnionFind(cgcComp)) {}

void Dimension0::computePairsAndMatch(vector<Cube>& ctr0, vector<Cube>& ctr1, vector<Cube>& ctrComp) {
	#ifdef RUNTIME
	cout << "input 0:" << endl;
	#endif
	#ifndef USE_CLEARING_DIM_0
	ctr0.clear();
	enumerateEdges(cgc0, ctr0);
	#endif
	computePairs(ctr0, 0);

    #ifdef RUNTIME
	cout << "input 1:" << endl;
	#endif
	#ifndef USE_CLEARING_DIM_0
	ctr1.clear();
	enumerateEdges(cgc1, ctr1);
	#endif
    computePairs(ctr1, 1);

	#ifdef RUNTIME
	cout << "comparison image, image 0, image 1 & matching:" << endl;
	#endif
	#ifndef USE_CLEARING_DIM_0
	ctrComp.clear();
	enumerateEdges(cgcComp, ctrComp);
	#endif
	uf0.reset();
	uf1.reset();
	computeImagePairsAndMatch(ctrComp);
}

void Dimension0::computePairs(vector<Cube>& edges, uint8_t k) {
	#ifdef RUNTIME
	cout <<"barcode: ";
	auto start = high_resolution_clock::now();
	#endif
	const CubicalGridComplex* const cgc = (k == 0) ? cgc0 : cgc1;
	UnionFind& uf = (k == 0) ? uf0 : uf1; 
	vector<Pair>& pairs = (k == 0) ? pairs0 : pairs1;
	unordered_map<index_t, Pair>& matchMap = (k == 0) ? matchMap0 : matchMap1;

	vector<index_t> boundaryIndices(2);
	index_t parentIdx0;
	index_t parentIdx1;
	index_t birthIdx;
	value_t birth;
	vector<index_t> birthCoordinates(3);
	for (auto edge = edges.begin(), last = edges.end(); edge != last; ++edge) {
		boundaryIndices = uf.getBoundaryIndices(*edge);
		parentIdx0 = uf.find(boundaryIndices[0]);
		parentIdx1 = uf.find(boundaryIndices[1]);
		if (parentIdx0 != parentIdx1) {
			birthIdx = uf.link(parentIdx0, parentIdx1);
			birth = uf.getBirth(birthIdx);
			if (birth != edge->birth) {
				birthCoordinates = uf.getCoordinates(birthIdx);
				pairs.push_back(Pair(Cube(birth, birthCoordinates[0], birthCoordinates[1], birthCoordinates[2], 0), *edge));
				matchMap.emplace(birthIdx, pairs.back());
			}
		}
	}
	#ifdef RUNTIME
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<milliseconds>(stop - start);
	cout << duration.count() << " ms" << endl;
	#endif
}

void Dimension0::computeImagePairsAndMatch(vector<Cube>& edges) {
	#ifdef RUNTIME
	cout << "barcodes & matching: ";
	auto start = high_resolution_clock::now();
	#endif
	vector<index_t> boundaryIndices(2);
	index_t parentIdx0;
	index_t parentIdx1;
	index_t birthIdx0;
	index_t birthIdx1;
	index_t birthIdxComp;
	value_t birth;
	vector<index_t> birthCoordinates(3);
	for (auto edge = edges.begin(), last = edges.end(); edge != last; ++edge) {
		boundaryIndices = ufComp.getBoundaryIndices(*edge);
		parentIdx0 = ufComp.find(boundaryIndices[0]);
		parentIdx1 = ufComp.find(boundaryIndices[1]);
		if (parentIdx0 != parentIdx1) {
			birthIdxComp = ufComp.link(parentIdx0, parentIdx1);
			birth = ufComp.getBirth(birthIdxComp);
			parentIdx0 = uf0.find(boundaryIndices[0]);
			parentIdx1 = uf0.find(boundaryIndices[1]);
			birthIdx0 = uf0.link(parentIdx0, parentIdx1);
			parentIdx0 = uf1.find(boundaryIndices[0]);
			parentIdx1 = uf1.find(boundaryIndices[1]);
			birthIdx1 = uf1.link(parentIdx0, parentIdx1);
			if (birth != edge->birth) {
				birthCoordinates = ufComp.getCoordinates(birthIdxComp);
				#ifdef COMPUTE_COMPARISON
				pairsComp.push_back(Pair(Cube(birth, birthCoordinates[0], birthCoordinates[1], birthCoordinates[2], 0), *edge));
				#endif
				auto find0 = matchMap0.find(birthIdx0);
				auto find1 = matchMap1.find(birthIdx1);
				if (find0 != matchMap0.end() && find1 != matchMap1.end()) {
					matches.push_back(Match(find0->second, find1->second));
					isMatched0.emplace(find0->second.birth.index, true);
					isMatched1.emplace(find1->second.birth.index, true);
				}
			} 
		}
	}
	#ifdef RUNTIME
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<milliseconds>(stop - start);
	cout << duration.count() << " ms" << endl;
	#endif
}

void Dimension0::enumerateEdges(const CubicalGridComplex* const cgc, vector<Cube>& edges) const {
	#ifdef RUNTIME
	cout << "enumeration: ";
	auto start = high_resolution_clock::now();
	#endif
	edges.reserve(cgc->getNumberOfCubes(1));
	value_t birth;
	for (index_t x = 0; x < cgc->shape[0]; x++) {
		for (index_t y = 0; y < cgc->shape[1]; y++) {
			for (index_t z = 0; z < cgc->shape[2]; z++) {
				for (uint8_t type = 0; type < 3; type++) {
					birth = cgc->getBirth(x, y, z, type, 1);
					if (birth < config.threshold) { edges.push_back(Cube(birth, x, y, z, type)); }	
				}				
			}
		}
	}
	sort(edges.begin(), edges.end(), CubeComparator());
	#ifdef RUNTIME
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<milliseconds>(stop - start);
	cout << duration.count() << " ms" << endl;
	#endif
}