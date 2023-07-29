#include "dimension_2.h"

#include <iostream>
#include <chrono>

using namespace std;
using namespace std::chrono;


Dimension2::Dimension2(const CubicalGridComplex* const _cgc0, const CubicalGridComplex* const _cgc1, 
						const CubicalGridComplex* const _cgcComp, const Config& _config, 
						vector<Pair>& _pairs0, vector<Pair>& _pairs1, vector<Pair>& _pairsComp, vector<Match>& _matches, 
						unordered_map<uint64_t, bool>& _isMatched0, unordered_map<uint64_t, bool>& _isMatched1) :
						cgc0(_cgc0), cgc1(_cgc1), cgcComp(_cgcComp), config(_config), 
						pairs0(_pairs0), pairs1(_pairs1), pairsComp(_pairsComp),
						matches(_matches), isMatched0(_isMatched0), isMatched1(_isMatched1), 
						uf0(UnionFindDual(cgc0)), uf1(UnionFindDual(cgc1)), ufComp(UnionFindDual(cgcComp)) {}

void Dimension2::computePairsAndMatch(vector<Cube>& ctr0, vector<Cube>& ctr1, vector<Cube>& ctrComp) {
	if (config.verbose) { cout << "computing dimension 2 ... "; }
    auto start = high_resolution_clock::now();
	
	enumerateDualEdges(cgc0, ctr0);
	enumerateDualEdges(cgc1, ctr1);
	enumerateDualEdges(cgcComp, ctrComp);

    computeImagePairs(ctr0, 0);
	ufComp.reset();
    computeImagePairs(ctr1, 1);
	ufComp.reset();
	computePairsCompAndMatch(ctrComp);

	auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);
    if (config.verbose) { cout << "took " << duration.count() << " ms" << endl; }
}

void Dimension2::enumerateDualEdges(const CubicalGridComplex* const cgc, vector<Cube>& dualEdges) const {
	dualEdges.clear();
	dualEdges.reserve(cgc->getNumberOfCubes(2));
	value_t birth;
	for (index_t x = 0; x < cgc->shape[0]; x++) {
		for (index_t y = 0; y < cgc->shape[1]; y++) {
			for (index_t z = 0; z < cgc->shape[2]; z++) {
				for (uint8_t type = 0; type < 3; type++) {
					birth = cgc->getBirth(x, y, z, type, 2);
					if (birth < config.threshold) { dualEdges.push_back(Cube(birth, x, y, z, type)); }
				}
			}
		}
	}
	sort(dualEdges.begin(), dualEdges.end(), CubeComparator());
}

void Dimension2::computeImagePairs(vector<Cube>& dualEdges, uint8_t k) {
	const CubicalGridComplex* const cgc = (k == 0) ? cgc0 : cgc1;
	UnionFindDual& uf = (k == 0) ? uf0 : uf1; 
	vector<Pair>& pairs = (k == 0) ? pairs0 : pairs1;
	unordered_map<index_t, Pair>& matchMap = (k==0) ? matchMap0 : matchMap1;
	
	vector<index_t> boundaryIndices;
	index_t parentIdx0;
	index_t parentIdx1;
	index_t birthIdx;
	index_t birthIdxComp;
	value_t birth;
	vector<index_t> birthCoordinates;
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
				pairs.push_back(Pair(*edge, Cube(birth, birthCoordinates[0], birthCoordinates[1], birthCoordinates[2], 0)));
				matchMap.emplace(birthIdxComp, pairs.back());
			} 
			edge->index = NONE;
		}
	}
	auto new_end = remove_if(dualEdges.begin(), dualEdges.end(), [](const Cube& cube){ return cube.index == NONE; });
	dualEdges.erase(new_end, dualEdges.end());
}

void Dimension2::computePairsCompAndMatch(vector<Cube>& dualEdges) {	
	vector<index_t> boundaryIndices;
	index_t parentIdx0;
	index_t parentIdx1;
	index_t birthIdx;
	value_t birth;
	vector<index_t> birthCoordinates;
	for (auto edge = dualEdges.rbegin(), last = dualEdges.rend(); edge != last; ++edge){
		boundaryIndices = ufComp.getBoundaryIndices(*edge);
		parentIdx0 = ufComp.find(boundaryIndices[0]);
		parentIdx1 = ufComp.find(boundaryIndices[1]);
		if (parentIdx0 != parentIdx1) {
			birthIdx = ufComp.link(parentIdx0, parentIdx1);
			birth = ufComp.getBirth(birthIdx);
			if (edge->birth != birth) {
				birthCoordinates = ufComp.getCoordinates(birthIdx);
				pairsComp.push_back(Pair(*edge, Cube(birth, birthCoordinates[0], birthCoordinates[1], birthCoordinates[2], 0)));
				auto find0 = matchMap0.find(birthIdx);
				auto find1 = matchMap1.find(birthIdx);
				if (find0 != matchMap0.end() && find1 != matchMap1.end()) {
					matches.push_back(Match(find0->second, find1->second));
					isMatched0.emplace(find0->second.birth.index, true);
					isMatched1.emplace(find1->second.birth.index, true);
				}
			}
			edge->index = NONE;
		}
	}
	auto new_end = std::remove_if(dualEdges.begin(), dualEdges.end(), [](const Cube& cube){ return cube.index == NONE; });
	dualEdges.erase(new_end, dualEdges.end());
}