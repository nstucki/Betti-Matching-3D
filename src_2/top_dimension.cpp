#include "top_dimension.h"
#include "enumerators.h"


TopDimension::TopDimension(const CubicalGridComplex& _cgc0, const CubicalGridComplex& _cgc1, const CubicalGridComplex& _cgcComp,
					const Config& _config, vector<Pair>& _pairs0, vector<Pair>& _pairs1, vector<Pair>& _pairsComp,
					vector<Match>& _matches, unordered_map<uint64_t, bool>& _isMatched0, unordered_map<uint64_t, bool>& _isMatched1) : 
					cgc0(_cgc0), cgc1(_cgc1), cgcComp(_cgcComp), pairs0(_pairs0), pairs1(_pairs1), pairsComp(_pairsComp), 
					matches(_matches), isMatched0(_isMatched0), isMatched1(_isMatched1), config(_config) {}

void TopDimension::enumerateDualEdges(const CubicalGridComplex& cgc, vector<Cube>& dualEdges) const {
	dualEdges.clear();
	dualEdges.reserve(cgc.getNumberOfCubes(cgc.dim-1));
	CubeEnumerator cubeEnum(cgc, cgc.dim-1);
	Cube cube = cubeEnum.getNextCube();
	if (cube.birth <= config.threshold) {
		dualEdges.push_back(cube);
	}
	while (cubeEnum.hasNextCube()) {
		cube = cubeEnum.getNextCube();
		if (cube.birth <= config.threshold) {
			dualEdges.push_back(cube);
		}
	}
	sort(dualEdges.begin(), dualEdges.end(), CubeComparator());
}

void TopDimension::computePairsComp(vector<Cube>& dualEdges) {
	UnionFindDual uf(cgcComp);
	vector<uint64_t> boundaryCoordinates;
	uint64_t degenAxis;
	uint64_t boundaryIdx0;
	uint64_t boundaryIdx1;
	uint64_t parentIdx0;
	uint64_t parentIdx1;
	uint64_t birthIdx;
	float birth;
	for (auto edge = dualEdges.rbegin(), last = dualEdges.rend(); edge != last; ++edge) {
		boundaryCoordinates = edge->coordinates;
		for (uint64_t i = 0; i < cgcComp.dim; i++) {
			if (boundaryCoordinates[i]%2 == 0) {
				degenAxis = i;
				break;
			}
		}
		if (boundaryCoordinates[degenAxis] == 0) {
			boundaryIdx0 = uf.star;
			boundaryCoordinates[degenAxis] += 1;
			boundaryIdx1 = uf.getIndex(boundaryCoordinates);
		} else if (boundaryCoordinates[degenAxis] == 2*cgcComp.shape[degenAxis]-2) {
			boundaryCoordinates[degenAxis] -= 1;
			boundaryIdx0 = uf.getIndex(boundaryCoordinates);
			boundaryIdx1 = uf.star;
		} else {
			boundaryCoordinates[degenAxis] -= 1;
			boundaryIdx0 = uf.getIndex(boundaryCoordinates);
			boundaryCoordinates[degenAxis] += 2;
			boundaryIdx1 = uf.getIndex(boundaryCoordinates);
		}
		parentIdx0 = uf.find(boundaryIdx0);
		parentIdx1 = uf.find(boundaryIdx1);
		if (parentIdx0 != parentIdx1) {
			birthIdx = uf.link(parentIdx0, parentIdx1);
			birth = uf.getBirth(birthIdx);
			if (edge->birth != birth) {
				pairsComp.push_back(Pair(*edge, Cube(birth, uf.getCoordinates(birthIdx))));
			}
			edge->coordinates[0] = NONE;
		}
	}
	auto new_end = remove_if(dualEdges.begin(), dualEdges.end(), [](const Cube& edge){ return edge.coordinates[0] == NONE; });
	dualEdges.erase(new_end, dualEdges.end());
}

void TopDimension::computeImagePairs(vector<Cube>& dualEdges, uint8_t k) {
	const CubicalGridComplex& cgc = (k == 0) ? cgc0 : cgc1;
	vector<Pair>& pairs = (k == 0) ? pairs0 : pairs1;
	unordered_map<uint64_t,Pair>& matchMap = (k==0) ? matchMap0 : matchMap1;
	
	UnionFindDual uf(cgc);
	UnionFindDual ufComp(cgcComp);
	vector<uint64_t> boundaryCoordinates;
	uint64_t degenAxis;
	uint64_t boundaryIdx0;
	uint64_t boundaryIdx1;
	uint64_t parentIdx0;
	uint64_t parentIdx1;
	uint64_t birthIdx;
	uint64_t birthIdxComp;
	float birth;
	for (auto edge = dualEdges.rbegin(), last = dualEdges.rend(); edge != last; ++edge) {
		boundaryCoordinates = edge->coordinates;
		for (uint64_t i = 0; i < cgc.dim; i++) {
			if (boundaryCoordinates[i]%2 == 0) {
				degenAxis = i;
				break;
			}
		}
		if (boundaryCoordinates[degenAxis] == 0) {
			boundaryIdx0 = uf.star;
			boundaryCoordinates[degenAxis] += 1;
			boundaryIdx1 = uf.getIndex(boundaryCoordinates);
		} else if (boundaryCoordinates[degenAxis] == 2*cgc.shape[degenAxis]-2) {
			boundaryCoordinates[degenAxis] -= 1;
			boundaryIdx0 = uf.getIndex(boundaryCoordinates);
			boundaryIdx1 = uf.star;
		} else {
			boundaryCoordinates[degenAxis] -= 1;
			boundaryIdx0 = uf.getIndex(boundaryCoordinates);
			boundaryCoordinates[degenAxis] += 2;
			boundaryIdx1 = uf.getIndex(boundaryCoordinates);
		}

		parentIdx0 = uf.find(boundaryIdx0);
		parentIdx1 = uf.find(boundaryIdx1);
		if (parentIdx0 != parentIdx1) {
			birthIdx = uf.link(parentIdx0, parentIdx1);
			birth = uf.getBirth(birthIdx);
			parentIdx0 = ufComp.find(boundaryIdx0);
			parentIdx1 = ufComp.find(boundaryIdx1);
			birthIdxComp = ufComp.link(parentIdx0, parentIdx1);
			if (edge->birth != birth) {
				pairs.push_back(Pair(*edge, Cube(birth, uf.getCoordinates(birthIdx))));
				matchMap.emplace(cgcComp.getCubeIndex(ufComp.getCoordinates(birthIdxComp)), pairs.back());
			} 
			edge->coordinates[0] = NONE;
		}
	}
	auto new_end = remove_if(dualEdges.begin(), dualEdges.end(), [](const Cube& edge){ return edge.coordinates[0] == NONE; });
	dualEdges.erase(new_end, dualEdges.end());
}

void TopDimension::computeMatching() {
	Pair pair0;
	Pair pair1;
	for (auto& pair : pairsComp) {
		auto find0 = matchMap0.find(cgcComp.getCubeIndex(pair.death));
		auto find1 = matchMap1.find(cgcComp.getCubeIndex(pair.death));
		if (find0 != matchMap0.end() && find1 != matchMap1.end()) {
			pair0 = (find0->second);
			pair1 = (find1->second);
			matches.push_back(Match(pair0, pair1));
			isMatched0.emplace(cgc0.getCubeIndex(pair0.birth), true);
			isMatched1.emplace(cgc1.getCubeIndex(pair1.birth), true);
		}
	}
}

void TopDimension::computePairsAndMatch(vector<Cube>& ctr0, vector<Cube>& ctr1, vector<Cube>& ctrComp) {
	enumerateDualEdges(cgcComp, ctrComp);
	computePairsComp(ctrComp);
	
	enumerateDualEdges(cgc0, ctr0);
    computeImagePairs(ctr0, 0);

	enumerateDualEdges(cgc1, ctr1);
    computeImagePairs(ctr1, 1);

    computeMatching();
}