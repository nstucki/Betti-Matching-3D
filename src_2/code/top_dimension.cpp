#include "template_functions.h"

#include "top_dimension.h"
#include "enumerators.h"

#include <iostream>

using namespace std;



TopDimension::TopDimension(const CubicalGridComplex& _cgc0, const CubicalGridComplex& _cgc1, const CubicalGridComplex& _cgcComp, 
							const Config& _config) : cgc0(_cgc0), cgc1(_cgc1), cgcComp(_cgcComp), config(_config) {}


void TopDimension::enumerateDualEdges(const CubicalGridComplex& cgc, vector<Cube>& edges) const {
	vector<uint64_t> shape = cgc.shape;
	uint64_t product = 1;
	for (auto& s : cgc.shape) {
		product *= s-1;
	}
	uint64_t numEdges = 0;
	for (auto& s : cgc.shape) {
		numEdges += s * (product/(s-1));
	}

	edges.clear();
	edges.reserve(numEdges);
	DualEdgeEnumerator enumerator = DualEdgeEnumerator(cgc);
	edges.push_back(enumerator.getCube());
	while (enumerator.hasNextCube()) {
		if (enumerator.getCube().birth <= config.threshold) {
			edges.push_back(enumerator.getCube());
		}
	}
	sort(edges.begin(), edges.end(), CubeComparator());
}


void TopDimension::computePairsComp(vector<Cube>& ctr) {
	enumerateDualEdges(cgcComp, ctr);
	UnionFindDual uf = UnionFindDual(cgcComp);
	uint64_t degenAxis;
	vector<uint64_t> boundaryCoordinates;
	uint64_t boundaryIdx0;
	uint64_t boundaryIdx1;
	uint64_t birthIdx0;
	uint64_t birthIdx1;
	uint64_t birthIdx;
	float birth;

	for (auto edge = ctr.rbegin(), last = ctr.rend(); edge != last; ++edge) {
		boundaryCoordinates = edge->coordinates;
		for (uint64_t i = 0; i < cgcComp.dim; i++) {
			if (boundaryCoordinates[i]%2 == 0) {
				degenAxis = i;
				break;
			}
		}
		if (boundaryCoordinates[degenAxis] == 0) {
			boundaryIdx0 = uf.n;
			boundaryCoordinates[degenAxis] += 1;
			boundaryIdx1 = uf.getIndex(boundaryCoordinates);
		} else if (boundaryCoordinates[degenAxis] == 2*cgcComp.dim-2) {
			boundaryCoordinates[degenAxis] -= 1;
			boundaryIdx0 = uf.getIndex(boundaryCoordinates);
			boundaryIdx1 = uf.n;
		} else {
			boundaryCoordinates[degenAxis] -= 1;
			boundaryIdx0 = uf.getIndex(boundaryCoordinates);
			boundaryCoordinates[degenAxis] += 2;
			boundaryIdx1 = uf.getIndex(boundaryCoordinates);
		}
		birthIdx0 = uf.find(boundaryIdx0);
		birthIdx1 = uf.find(boundaryIdx1);

		if (birthIdx0 != birthIdx1) {
			birthIdx = uf.link(birthIdx0, birthIdx1);
			birth = uf.getBirth(birthIdx);
			if (edge->birth != birth) {
				pairsComp.push_back(Pair(*edge, Cube(birth, uf.getCoordinates(birthIdx))));
			}
			edge->coordinates[0] = NONE;
		}
	}
	auto new_end = remove_if(ctr.begin(), ctr.end(), [](const Cube& edge){ return edge.coordinates[0] == NONE; });
	ctr.erase(new_end, ctr.end());
}


void TopDimension::computePairsImage(uint8_t k, vector<Cube>& ctr) {
	const CubicalGridComplex& cgc = (k == 0) ? cgc0 : cgc1;
	vector<Pair>& pairs = (k == 0) ? pairs0 : pairs1;
	unordered_map<uint64_t,Pair>& matchMap = (k==0) ? matchMap0 : matchMap1;
	
	enumerateDualEdges(cgc, ctr);
	UnionFindDual uf = UnionFindDual(cgc);
	UnionFindDual ufComp = UnionFindDual(cgcComp);
	uint64_t degenAxis;
	vector<uint64_t> boundaryCoordinates;
	uint64_t boundaryIdx0;
	uint64_t boundaryIdx1;
	uint64_t birthIdx0;
	uint64_t birthIdx1;
	uint64_t birthIdx;
	uint64_t birthIdxComp;
	float birth;

	for (auto edge = ctr.rbegin(), last = ctr.rend(); edge != last; ++edge) {
		boundaryCoordinates = edge->coordinates;
		for (uint64_t i = 0; i < cgc.dim; i++) {
			if (boundaryCoordinates[i]%2 == 0) {
				degenAxis = i;
				break;
			}
		}
		if (boundaryCoordinates[degenAxis] == 0) {
			boundaryIdx0 = uf.n;
			boundaryCoordinates[degenAxis] += 1;
			boundaryIdx1 = uf.getIndex(boundaryCoordinates);
		} else if (boundaryCoordinates[degenAxis] == 2*cgc.dim-2) {
			boundaryCoordinates[degenAxis] -= 1;
			boundaryIdx0 = uf.getIndex(boundaryCoordinates);
			boundaryIdx1 = uf.n;
		} else {
			boundaryCoordinates[degenAxis] -= 1;
			boundaryIdx0 = uf.getIndex(boundaryCoordinates);
			boundaryCoordinates[degenAxis] += 2;
			boundaryIdx1 = uf.getIndex(boundaryCoordinates);
		}
		birthIdx0 = uf.find(boundaryIdx0);
		birthIdx1 = uf.find(boundaryIdx1);

		if (birthIdx0 != birthIdx1) {
			birthIdx = uf.link(birthIdx0, birthIdx1);
			birth = uf.getBirth(birthIdx);

			birthIdx0 = ufComp.find(boundaryIdx0);
			birthIdx1 = ufComp.find(boundaryIdx1);
			birthIdxComp = ufComp.link(birthIdx0, birthIdx1);

			if (edge->birth != birth) {
				pairs.push_back(Pair(*edge, Cube(birth, uf.getCoordinates(birthIdx))));
				matchMap.emplace(cgcComp.getCubeIndex(ufComp.getCoordinates(birthIdxComp)), pairs.back());
			} 
			edge->coordinates[0] = NONE;
		}
	}
	auto new_end = remove_if(ctr.begin(), ctr.end(), [](const Cube& edge){ return edge.coordinates[0] == NONE; });
	ctr.erase(new_end, ctr.end());
}


void TopDimension::computeMatching() {
	for (auto& pair : pairsComp) {
		auto find0 = matchMap0.find(cgcComp.getCubeIndex(pair.death));
		auto find1 = matchMap1.find(cgcComp.getCubeIndex(pair.death));
		if (find0 != matchMap0.end() && find1 != matchMap1.end()) {
			matches.push_back(Match((find0 -> second), (find1 -> second)));
		}
	}
}


void TopDimension::computePairsAndMatch(vector<Cube>& ctr0, vector<Cube>& ctr1, vector<Cube>& ctrComp) {
	computePairsComp(ctrComp);
    computePairsImage(0, ctr0);
    computePairsImage(1, ctr1);
    computeMatching();
}