#include "dimension_0.h"
#include "enumerators.h"

#include <iostream>
#include <chrono>

using namespace std;
using namespace std::chrono;

Dimension0::Dimension0(const CubicalGridComplex& _cgc0, const CubicalGridComplex& _cgc1, const CubicalGridComplex& _cgcComp,
						const Config& _config, vector<Pair>& _pairs0, vector<Pair>& _pairs1, vector<Pair>& _pairsComp,
						vector<Match>& _matches, unordered_map<index_t, bool>& _isMatched0, unordered_map<index_t, bool>& _isMatched1) : 
						cgc0(_cgc0), cgc1(_cgc1), cgcComp(_cgcComp), pairs0(_pairs0), pairs1(_pairs1), pairsComp(_pairsComp), 
						matches(_matches), isMatched0(_isMatched0), isMatched1(_isMatched1), config(_config) {}

void Dimension0::enumerateEdges(const CubicalGridComplex& cgc, vector<Cube>& edges) const {
	edges.clear();
	edges.reserve(cgc.getNumberOfCubes(1));
	CubeEnumerator cubeEnum(cgc, 1);
	Cube cube = cubeEnum.getNextCube();
	if (cube.birth <= config.threshold) {
		edges.push_back(cube);
	}
	while (cubeEnum.hasNextCube()) {
		cube = cubeEnum.getNextCube();
		if (cube.birth <= config.threshold) {
			edges.push_back(cube);
		}
	}
	sort(edges.begin(), edges.end(), CubeComparator());
}

void Dimension0::computePairs(vector<Cube>& edges, uint8_t k) {
	const CubicalGridComplex& cgc = (k == 0) ? cgc0 : cgc1;
	vector<Pair>& pairs = (k == 0) ? pairs0 : pairs1;
	unordered_map<index_t, Pair>& matchMap = (k == 0) ? matchMap0 : matchMap1;

	UnionFind uf(cgc);
	vector<index_t> boundaryCoordinates;
	vector<index_t> birthCoordinates;
	index_t nonDegenAxis;
	index_t boundaryIdx0;
	index_t boundaryIdx1;
	index_t birthIdx0;
	index_t birthIdx1;
	index_t birthIdx;
	float_t birth;
	for (auto& edge : edges) {
		boundaryCoordinates = edge.coordinates;
		for (index_t i = 0; i < cgcComp.dim; i++) {
			if (boundaryCoordinates[i]%2 == 1) {
				nonDegenAxis = i;
				break;
			}
		}
		boundaryCoordinates[nonDegenAxis] -= 1;
		boundaryIdx0 = uf.getIndex(boundaryCoordinates);
		boundaryCoordinates[nonDegenAxis] += 2;
		boundaryIdx1 = uf.getIndex(boundaryCoordinates);
		birthIdx0 = uf.find(boundaryIdx0);
		birthIdx1 = uf.find(boundaryIdx1);
		if (birthIdx0 != birthIdx1) {
			birthIdx = uf.link(birthIdx0, birthIdx1);
			birth = uf.getBirth(birthIdx);
			if (birth != edge.birth) {
				birthCoordinates = uf.getCoordinates(birthIdx);
				pairs.push_back(Pair(Cube(birth, birthCoordinates), edge));
				matchMap.emplace(cgc.getCubeIndex(birthCoordinates), pairs.back());
			}
		}
	}
}

void Dimension0::computeImagePairsAndMatch(vector<Cube>& edges) {
	UnionFind uf0(cgc0);
	UnionFind uf1(cgc1);
	UnionFind ufComp(cgcComp);
	vector<index_t> boundaryCoordinates;
	Pair pair0;
	Pair pair1;
	index_t nonDegenAxis;
	index_t boundaryIdx0;
	index_t boundaryIdx1;
	index_t parentIdx0;
	index_t parentIdx1;
	index_t birthIdx0;
	index_t birthIdx1;
	index_t birthIdxComp;
	float_t birth;
	for (auto& edge : edges) {
		boundaryCoordinates = edge.coordinates;
		for (index_t i = 0; i < cgcComp.dim; i++) {
			if (boundaryCoordinates[i]%2 == 1) {
				nonDegenAxis = i;
				break;
			}
		}
		boundaryCoordinates[nonDegenAxis] -= 1;
		boundaryIdx0 = ufComp.getIndex(boundaryCoordinates);
		boundaryCoordinates[nonDegenAxis] += 2;
		boundaryIdx1 = ufComp.getIndex(boundaryCoordinates);
		parentIdx0 = ufComp.find(boundaryIdx0);
		parentIdx1 = ufComp.find(boundaryIdx1);
		if (parentIdx0 != parentIdx1) {
			birthIdxComp = ufComp.link(parentIdx0, parentIdx1);
			birth = ufComp.getBirth(birthIdxComp);
			parentIdx0 = uf0.find(boundaryIdx0);
			parentIdx1 = uf0.find(boundaryIdx1);
			birthIdx0 = uf0.link(parentIdx0, parentIdx1);
			parentIdx0 = uf1.find(boundaryIdx0);
			parentIdx1 = uf1.find(boundaryIdx1);
			birthIdx1 = uf1.link(parentIdx0, parentIdx1);
			if (birth != edge.birth) {
				pairsComp.push_back(Pair(Cube(birth, ufComp.getCoordinates(birthIdxComp)), edge));
				birthIdx0 = cgc0.getCubeIndex(uf0.getCoordinates(birthIdx0));
				birthIdx1 = cgc1.getCubeIndex(uf1.getCoordinates(birthIdx1));
				auto find0 = matchMap0.find(birthIdx0);
				auto find1 = matchMap1.find(birthIdx1);
				if (find0 != matchMap0.end() && find1 != matchMap1.end()) {
					pair0 = find0->second;
					pair1 = find1->second;
					matches.push_back(Match(pair0, pair1));
					isMatched0.emplace(birthIdx0, true);
					isMatched1.emplace(birthIdx1, true);
				}
			} 
		}
	}
}

void Dimension0::computePairsAndMatch(vector<Cube>& ctr0, vector<Cube>& ctr1, vector<Cube>& ctrComp) {
	if (config.verbose) { cout << "computing dimension 0 ... "; }
    auto start = high_resolution_clock::now();

	enumerateEdges(cgc0, ctr0);
	computePairs(ctr0, 0);
	
	enumerateEdges(cgc1, ctr1);
    computePairs(ctr1, 1);

	enumerateEdges(cgcComp, ctrComp);
    computeImagePairsAndMatch(ctrComp);

	auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);
    if (config.verbose) { cout << "took " << duration.count() << " ms" << endl; }
}