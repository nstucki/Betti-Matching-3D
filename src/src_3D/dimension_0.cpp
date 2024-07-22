#include "dimension_0.h"
#include "BettiMatching.h"
#include "../utils.h"

#include <algorithm>
#include <chrono>
#include <iostream>
#include <stdexcept>
#include <unordered_map>

using namespace dim3;
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

#ifndef USE_CLEARING_DIM0
	enumerateEdges(ctr0, cgc0);
#endif
	computePairs(ctr0, 0);

#ifdef RUNTIME
	cout << endl << "input 1: ";
#endif

#ifndef USE_CLEARING_DIM0
	enumerateEdges(ctr1, cgc1);
#endif
    computePairs(ctr1, 1);

#ifdef RUNTIME
	cout << endl << "comparison & image 0 & image 1 & match: ";
#endif

#ifndef USE_CLEARING_DIM0
	enumerateEdges(ctrComp, cgcComp);
#endif
	computeImagePairsAndMatch(ctrComp);
}

void Dimension0::computeInput0Pairs(vector<Cube>& ctr0) {
#ifndef USE_CLEARING_DIM0
	enumerateEdges(ctr0, cgc0);
#endif
	computePairs(ctr0, 0);
}

vector<vector<index_t>>
Dimension0::getRepresentativeCycle(const Pair &pair, const CubicalGridComplex &cgc) const {
    vector<Cube> edges;
    enumerateEdges(edges, cgc);
    UnionFind uf(cgc);
    vector<index_t> boundaryIndices(2);
    index_t parentIdx0;
    index_t parentIdx1;
    index_t birthIdx;

    for (Cube &edge : edges) {
    if (edge == pair.death) {
        break;
    }
    boundaryIndices = uf.getBoundaryIndices(edge);
    parentIdx0 = uf.find(boundaryIndices[0]);
    parentIdx1 = uf.find(boundaryIndices[1]);
        if (parentIdx0 != parentIdx1) {
            birthIdx = uf.link(parentIdx0, parentIdx1);
        }
    }

    vector<vector<index_t>> reprCycle;
    reprCycle.push_back(cgc.getParentVoxel(pair.birth, 0));
    vector<index_t> vertex;
    parentIdx0 = uf.find(pair.birth.x() * cgc.n_yz + pair.birth.y() * cgc.shape[2] + pair.birth.z());
    for (size_t i = 0; i < cgc.getNumberOfCubes(0); ++i) {
    parentIdx1 = uf.find(i);
        if (parentIdx0 == parentIdx1) {
            vertex = uf.getCoordinates(i);
            if (find(reprCycle.begin(), reprCycle.end(), vertex) == reprCycle.end()) {
                reprCycle.push_back(vertex);
            }
        }
    }
    reprCycle.push_back(cgc.getParentVoxel(pair.death, 1));

    return reprCycle;
}

tuple<vector<RepresentativeCycle>, vector<RepresentativeCycle>>
Dimension0::getAllRepresentativeCycles(uint8_t input, bool computeMatchedCycles, bool computeUnmatchedCycles) {
    const CubicalGridComplex &cgc = (input == 0) ? cgc0 : cgc1;
    UnionFind uf(cgc);
    vector<Pair> &pairs = (input == 0) ? pairs0 : pairs1;
    unordered_map<index_t, Pair> &matchMap = (input == 0) ? matchMap0 : matchMap1;
    auto &isMatched = (input == 0) ? isMatched0 : isMatched1;

    // Map from cube indices to union find indices (to match cycles with persistence pairs)
    CubeMap<1, index_t> unionFindIdxByDeathCube(cgc.shape);

    // Initialize singleton representative cycles
    vector<RepresentativeCycle> cycleByBirthIdx(cgc.shape[0] * cgc.shape[1] * cgc.shape[2]);
    for (int birthIdx = 0; birthIdx < cycleByBirthIdx.size(); birthIdx++) {
        cycleByBirthIdx[birthIdx].emplace_back(vectorToTuple<3>(uf.getCoordinates(birthIdx)));
    }

    vector<Cube> edges;
    enumerateEdges(edges, cgc);

    // Track the next pair to be found (we will find them in the same order as before) to know when to save a cycle
    auto currentPair = pairs.begin();

    // Follow the union-find algorithm:
    for (Cube &edge : edges) {
        vector<index_t> boundaryIndices = uf.getBoundaryIndices(edge);
        index_t parentIdx0 = uf.find(boundaryIndices[0]);
        index_t parentIdx1 = uf.find(boundaryIndices[1]);
        // When merging a younger component into an older one, extend the older component's cycle by the younger component's cycle
        if (parentIdx0 != parentIdx1) {
            auto youngerBirthIdx = uf.link(parentIdx0, parentIdx1);
            auto olderBirthIdx = (parentIdx1 == youngerBirthIdx) ? parentIdx0 : parentIdx1;
            cycleByBirthIdx[olderBirthIdx].insert(cycleByBirthIdx[olderBirthIdx].end(), cycleByBirthIdx[youngerBirthIdx].begin(), cycleByBirthIdx[youngerBirthIdx].end());
            unionFindIdxByDeathCube[edge.index] = youngerBirthIdx;

            // If the died component corresponds to a pair we'd like to save:
            if (currentPair != pairs.end() && edge.index == currentPair->death.index) {
                currentPair++;
            } else {
                // Else, delete the component's cycle to save memory
                cycleByBirthIdx[youngerBirthIdx].clear();
            }
        }
    }

    // Collect the representative cycles belonging to matched pairs
    vector<RepresentativeCycle> matchedCycles;
    if (computeMatchedCycles) {
        matchedCycles.reserve(matches.size());
        for (auto &match : matches) {
            auto &pair = (input == 0) ? match.pair0 : match.pair1;
            auto unionFindIdx = unionFindIdxByDeathCube.find(pair.death.index);
            if (!unionFindIdx.has_value()) {
                throw runtime_error("Union find index for matched pair cannot be found");
            }
            auto &cycle = cycleByBirthIdx[*unionFindIdx];
            matchedCycles.emplace_back(std::move(cycle));
        }
    }

    // Collect the representative cycles belonging to unmatched pairs
    vector<RepresentativeCycle> unmatchedCycles;
    if (computeUnmatchedCycles) {
        unmatchedCycles.reserve(pairs.size() - matches.size());
        for (auto &pair : pairs) {
            if (!isMatched[pair.birth.index]) {
                auto unionFindIdx = unionFindIdxByDeathCube.find(pair.death.index);
                if (!unionFindIdx.has_value()) {
                    throw runtime_error("Union find index for matched pair cannot be found");
                }
                auto &cycle = cycleByBirthIdx[*unionFindIdx];
                unmatchedCycles.emplace_back(std::move(cycle));
            }
        }
    }

    return {matchedCycles, unmatchedCycles};
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
	
	uf.reset();
	vector<index_t> boundaryIndices(2);
	index_t parentIdx0;
	index_t parentIdx1;
	index_t birthIdx;
	value_t birth;
	vector<index_t> birthCoordinates(3);

	for (Cube& edge : edges) {
		boundaryIndices = uf.getBoundaryIndices(edge);
		parentIdx0 = uf.find(boundaryIndices[0]);
		parentIdx1 = uf.find(boundaryIndices[1]);
		if (parentIdx0 != parentIdx1) {
			birthIdx = uf.link(parentIdx0, parentIdx1);
			birth = uf.getBirth(birthIdx);
			if (birth != edge.birth) {
				birthCoordinates = uf.getCoordinates(birthIdx);
				pairs.push_back(Pair(Cube(birth, birthCoordinates[0], birthCoordinates[1], birthCoordinates[2], 0), edge));
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
	cout << "barcodes & match ";
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
	vector<index_t> birthCoordinates(3);

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
				birthCoordinates = ufComp.getCoordinates(birthIdxComp);
#ifdef COMPUTE_COMPARISON
				pairsComp.push_back(Pair(Cube(birth, birthCoordinates[0], birthCoordinates[1], birthCoordinates[2], 0), edge));
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


void Dimension0::enumerateEdges(vector<Cube>& edges, const CubicalGridComplex& cgc) const {
#ifdef RUNTIME
	cout << "; enumeration ";
	auto start = high_resolution_clock::now();
#endif

	edges.clear();
	edges.reserve(cgc.getNumberOfCubes(1));
	value_t birth;
#ifdef USE_CLEARING_DIM0
	Cube cube;
#endif
#ifdef USE_STABLE_SORT_OR_STABLE_PARTITION
	bool binaryInputs = true;
#endif
	for (index_t x = 0; x < cgc.shape[0]; ++x) {
		for (index_t y = 0; y < cgc.shape[1]; ++y) {
			for (index_t z = 0; z < cgc.shape[2]; ++z) {
				for (uint8_t type = 0; type < 3; ++type) {
					birth = cgc.getBirth(x, y, z, type, 1);
					if (birth < config.threshold) {
						edges.push_back(Cube(birth, x, y, z, type));
#ifdef USE_STABLE_SORT_OR_STABLE_PARTITION
						if (binaryInputs && birth != 0 && birth != 1) binaryInputs = false;
#endif
					}	
				}				
			}
		}
	}

#ifdef USE_STABLE_SORT_OR_STABLE_PARTITION
	if (binaryInputs) {
		std::stable_partition(edges.begin(), edges.end(), [](Cube &cube) { return cube.birth == 0; });
	} else {
		std::stable_sort(edges.begin(), edges.end(), [](const Cube &cube1, const Cube &cube2) { return cube1.birth < cube2.birth; });
	}
#else
	std::sort(edges.begin(), edges.end(), CubeComparator());
#endif

#ifdef RUNTIME
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<milliseconds>(stop - start);
	cout << duration.count() << " ms";
#endif
}