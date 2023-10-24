#include "dimension_0.h"

#include <iostream>
#include <chrono>
#include <algorithm>

using namespace dim1;
using namespace std::chrono;



Dimension0::Dimension0(const CubicalGridComplex& _cgc0, const CubicalGridComplex& _cgc1, 
						const CubicalGridComplex& _cgcComp, const Config& _config, 
						vector<Pair>& _pairs0, vector<Pair>& _pairs1, vector<Pair>& _pairsComp, vector<Match>& _matches, 
						unordered_map<uint64_t, bool>& _isMatched0, unordered_map<uint64_t, bool>& _isMatched1) : 
						cgc0(_cgc0), cgc1(_cgc1), cgcComp(_cgcComp), config(_config), 
						pairs0(_pairs0), pairs1(_pairs1), pairsComp(_pairsComp),
						matches(_matches), isMatched0(_isMatched0), isMatched1(_isMatched1),
						uf0(UnionFind(cgc0)), uf1(UnionFind(cgc1)), ufComp(UnionFind(cgcComp)) {}


void Dimension0::computePairsAndMatch(vector<Cube>& ctr0, vector<Cube>& ctr1, vector<Cube>& ctrComp) {
#ifdef RUNTIME
	cout << endl << "input 0: ";
#endif

	enumerateEdges(ctr0, cgc0);
	computePairs(ctr0, 0);

#ifdef RUNTIME
	cout << endl << "input 1: ";
#endif

	enumerateEdges(ctr1, cgc1);
    computePairs(ctr1, 1);

#ifdef RUNTIME
	cout << endl << "comparison & image 0 & image 1 & matching: ";
#endif

	enumerateEdges(ctrComp, cgcComp);
	computeImagePairsAndMatch(ctrComp);
}


vector<vector<index_t>> Dimension0::getRepresentativeCycle(const Pair& pair, const CubicalGridComplex& cgc) const {
	vector<Cube> edges;
	enumerateEdges(edges, cgc);
	UnionFind uf(cgc);
	vector<index_t> boundaryIndices(2);
	index_t parentIdx0;
	index_t parentIdx1;
	index_t birthIdx;
	
	for (Cube& edge : edges) {
		if (edge == pair.death) { break; }
		boundaryIndices = uf.getBoundaryIndices(edge);
		parentIdx0 = uf.find(boundaryIndices[0]);
		parentIdx1 = uf.find(boundaryIndices[1]);
		if (parentIdx0 != parentIdx1) { birthIdx = uf.link(parentIdx0, parentIdx1); }
	}

	vector<vector<index_t>> reprCycle;
	reprCycle.push_back({cgc.getParentVoxel(pair.birth, 0)});
	vector<index_t> vertex;
	parentIdx0 = uf.find(pair.birth.x());
	for (size_t i = 0; i < cgc.getNumberOfCubes(0); ++i) {
		parentIdx1 = uf.find(i);
		if (parentIdx0 == parentIdx1) { 
			vertex = {uf.getCoordinate(i)};
			if(find(reprCycle.begin(), reprCycle.end(), vertex) == reprCycle.end()) { reprCycle.push_back(vertex); }
		}
	}
	reprCycle.push_back({cgc.getParentVoxel(pair.death, 1)});

	return reprCycle;
}


void Dimension0::enumerateEdges(vector<Cube>& edges, const CubicalGridComplex& cgc) const {
#ifdef RUNTIME
	cout << "enumeration ";
	auto start = high_resolution_clock::now();
#endif

	edges.reserve(cgc.getNumberOfCubes(1));
	value_t birth;
	for (index_t x = 0; x < cgc.m_x; ++x) {
		birth = cgc.getBirth(x, 1);
		if (birth < config.threshold) { edges.push_back(Cube(birth, x)); }
	}

	sort(edges.begin(), edges.end(), CubeComparator());

#ifdef RUNTIME
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<milliseconds>(stop - start);
	cout << duration.count() << " ms, ";
#endif
}


void Dimension0::computePairs(vector<Cube>& edges, uint8_t k) {
#ifdef RUNTIME
	cout << "barcode ";
	auto start = high_resolution_clock::now();
#endif

	const CubicalGridComplex& cgc = (k == 0) ? cgc0 : cgc1;
	UnionFind& uf = (k == 0) ? uf0 : uf1; 
	vector<Pair>& pairs = (k == 0) ? pairs0 : pairs1;
	unordered_map<index_t, Pair>& matchMap = (k == 0) ? matchMap0 : matchMap1;
	vector<index_t> boundaryIndices(2);
	index_t parentIdx0;
	index_t parentIdx1;
	index_t birthIdx;
	value_t birth;
	index_t birthCoordinate;

	for (Cube& edge : edges) {
		boundaryIndices = uf.getBoundaryIndices(edge);
		parentIdx0 = uf.find(boundaryIndices[0]);
		parentIdx1 = uf.find(boundaryIndices[1]);
		if (parentIdx0 != parentIdx1) {
			birthIdx = uf.link(parentIdx0, parentIdx1);
			birth = uf.getBirth(birthIdx);
			if (birth != edge.birth) {
				birthCoordinate = uf.getCoordinate(birthIdx);
				pairs.push_back(Pair(Cube(birth, birthCoordinate), edge));
				matchMap.emplace(birthIdx, pairs.back());
			}
		}
	}

#ifdef RUNTIME
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<milliseconds>(stop - start);
	cout << duration.count() << " ms";
#endif
}


void Dimension0::computeImagePairsAndMatch(vector<Cube>& edges) {
#ifdef RUNTIME
	cout << "barcodes & matching ";
	auto start = high_resolution_clock::now();
#endif

	uf0.reset();
	uf1.reset();
	vector<index_t> boundaryIndices(2);
	index_t parentIdx0;
	index_t parentIdx1;
	index_t birthIdx0;
	index_t birthIdx1;
	index_t birthIdxComp;
	value_t birth;
	index_t birthCoordinate;

	for (Cube& edge : edges) {
		boundaryIndices = ufComp.getBoundaryIndices(edge);
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
			if (birth != edge.birth) {
				birthCoordinate = ufComp.getCoordinate(birthIdxComp);
#ifdef COMPUTE_COMPARISON
				pairsComp.push_back(Pair(Cube(birth, birthCoordinate), edge));
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
	cout << duration.count() << " ms";
#endif
}