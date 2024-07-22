#include "dimension_2.h"
#include "BettiMatching.h"
#include "data_structures.h"
#include "enumerators.h"

#include <algorithm>
#include <chrono>
#include <iostream>
#include <set>
#include <stdexcept>
#include <tuple>
#include <unordered_map>

using namespace dim3;
using namespace std;
using namespace std::chrono;



Dimension2::Dimension2(const CubicalGridComplex& _cgc0, const CubicalGridComplex& _cgc1, 
						const CubicalGridComplex& _cgcComp, const Config& _config, 
						vector<Pair>& _pairs0, vector<Pair>& _pairs1, vector<Pair>& _pairsComp, vector<Match>& _matches, 
						unordered_map<uint64_t, bool>& _isMatched0, unordered_map<uint64_t, bool>& _isMatched1) :
						cgc0(_cgc0), cgc1(_cgc1), cgcComp(_cgcComp), config(_config), 
						pairs0(_pairs0), pairs1(_pairs1), pairsComp(_pairsComp),
						matches(_matches), isMatched0(_isMatched0), isMatched1(_isMatched1), 
						uf0(UnionFindDual(cgc0)), uf1(UnionFindDual(cgc1)), ufComp(UnionFindDual(cgcComp)) {}


void Dimension2::computePairsAndMatch(vector<Cube>& ctr0, vector<Cube>& ctr1, vector<Cube>& ctrComp, vector<Cube>& ctrImage) {
	size_t actualDim = 3;
	for (index_t s : cgc0.shape) { if (s == 1) { --actualDim; } }
	bool needToCompute = ( actualDim == 3);
	if (!needToCompute) {
#ifdef RUNTIME
		cout << endl << "input 0 & image 0: ";
#endif
		enumerateDualEdges(ctr0, cgc0);
#ifdef RUNTIME
		cout << "barcodes not computed";
#endif

#ifdef RUNTIME
		cout << endl << "input 1 & image 1: ";
#endif
		enumerateDualEdges(ctr1, cgc1);
#ifdef RUNTIME
		cout << "barcodes not computed";
#endif

#ifdef RUNTIME
		cout << endl << "comparison & matching: ";
#endif
		enumerateDualEdgesComp(ctrComp);
#ifdef RUNTIME
		cout << "barcode & matching not computed";
#endif
	} else {
#ifdef RUNTIME
		cout << endl << "input 0 & image 0: ";
#endif
		enumerateDualEdges(ctr0, cgc0);
		computeInputAndImagePairs(ctr0, 0);

#ifdef RUNTIME
		cout << endl << "input 1 & image 1: ";
#endif
		enumerateDualEdges(ctr1, cgc1);
		computeInputAndImagePairs(ctr1, 1);
    
#ifdef RUNTIME
	cout << endl << "comparison & match: ";
#endif
		enumerateDualEdgesComp(ctrComp);
		computeCompPairsAndMatch(ctrComp, ctrImage);
	}
}

void Dimension2::computeInput0Pairs(vector<Cube>& ctr0)  {
	size_t actualDim = 3;
	for (index_t s : cgc0.shape) { if (s == 1) { --actualDim; } }
	bool needToCompute = (actualDim == 3);
	if (!needToCompute) {
		enumerateDualEdges(ctr0, cgc0);
	} else {
		enumerateDualEdges(ctr0, cgc0);
		computeInputAndImagePairs(ctr0, 0);
	}
}

// See https://stackoverflow.com/a/28411055
template <typename T, std::size_t... Indices>
auto vectorToTupleHelper(const std::vector<T> &v,
                         std::index_sequence<Indices...>) {
        return std::make_tuple(v[Indices]...);
}

template <std::size_t N, typename T>
auto vectorToTuple(const std::vector<T> &v) {
  return vectorToTupleHelper(v, std::make_index_sequence<N>());
}

template <typename Coordinate>
vector<Coordinate> extractRepresentativeCycle(const Pair &pair, const CubicalGridComplex &cgc, UnionFindDual &uf) {
    set<tuple<index_t, index_t, index_t>> cubeCoordinates;
    index_t parentIdx0 = uf.find(pair.death.x() * cgc.m_yz + pair.death.y() * cgc.m_z + pair.death.z());
    for (size_t i = 0; i < cgc.getNumberOfCubes(3); ++i) {
        index_t parentIdx1 = uf.find(i);
        if (parentIdx0 == parentIdx1) {
            cubeCoordinates.insert(
                vectorToTuple<3>(uf.getCoordinates(i)));
        }
    }

    multiset<tuple<index_t, index_t, index_t>> boundaryVertices;
    for (const tuple<index_t, index_t, index_t> &c : cubeCoordinates) {
        for (uint8_t x = 0; x < 2; ++x) {
            for (uint8_t y = 0; y < 2; ++y) {
                for (uint8_t z = 0; z < 2; ++z) {
                    boundaryVertices.insert({std::get<0>(c) + x, std::get<1>(c) + y, std::get<2>(c) + z});
                }
            }
        }
    }

    vector<Coordinate> reprCycle;
    for (const tuple<index_t, index_t, index_t> &vertex : boundaryVertices) {
        auto lower = boundaryVertices.lower_bound(vertex);
        auto upper = boundaryVertices.upper_bound(vertex);
        int multiplicity = distance(lower, upper);
        if (multiplicity < 8) {
            reprCycle.push_back({std::get<0>(vertex), std::get<1>(vertex), std::get<2>(vertex)});
        }
    }
    auto deathVoxel = cgc.getParentVoxel(pair.death, 3);
    reprCycle.push_back({deathVoxel[0], deathVoxel[1], deathVoxel[2]});

    return reprCycle;
}

tuple<vector<RepresentativeCycle>, vector<RepresentativeCycle>> Dimension2::getAllRepresentativeCycles(uint8_t input, bool computeMatchedCycles, bool computeUnmatchedCycles) const {
	const CubicalGridComplex &cgc = (input == 0) ? cgc0 : cgc1;
	UnionFindDual uf(cgc);
	vector<Pair> &pairs = (input == 0) ? pairs0 : pairs1;
	const unordered_map<index_t, Pair> &matchMap =
		(input == 0) ? matchMap0 : matchMap1;
	auto &isMatched = (input == 0) ? isMatched0 : isMatched1;

	vector<Cube> dualEdges;
	enumerateDualEdges(dualEdges, cgc);

	unordered_map<uint64_t, RepresentativeCycle> matchedCyclesByBirth;
	vector<RepresentativeCycle> unmatchedCycles;
	auto currentPair = pairs.begin();

	for (auto edge = dualEdges.rbegin(), last = dualEdges.rend(); edge != last; ++edge) {
		if (currentPair != pairs.end() && *edge == currentPair->birth) {
			RepresentativeCycle reprCycle = extractRepresentativeCycle<tuple<uint64_t, uint64_t, uint64_t>>(*currentPair, cgc, uf);
			if (isMatched[currentPair->birth.index]) {
				matchedCyclesByBirth[currentPair->birth.index] = reprCycle;
			} else {
                if (computeUnmatchedCycles) {
    				unmatchedCycles.emplace_back(reprCycle);
                }
			}
			currentPair++;
		}
		vector<index_t> boundaryIndices = uf.getBoundaryIndices(*edge);
		index_t parentIdx0 = uf.find(boundaryIndices[0]);
		index_t parentIdx1 = uf.find(boundaryIndices[1]);
		if (parentIdx0 != parentIdx1) {
			uf.link(parentIdx0, parentIdx1);
		}
	}

	vector<RepresentativeCycle> matchedCycles;

    if (computeMatchedCycles) {
        for (auto& match : matches) {
            auto& pair = (input == 0) ? match.pair0 : match.pair1;
            matchedCycles.emplace_back(std::move(matchedCyclesByBirth[pair.birth.index]));
        }
    }

	return {matchedCycles, unmatchedCycles};
}

vector<vector<index_t>> Dimension2::getRepresentativeCycle(const Pair& pair, const CubicalGridComplex& cgc) const {
    vector<Cube> dualEdges;
    enumerateDualEdges(dualEdges, cgc);
    UnionFindDual uf(cgc);
    vector<index_t> boundaryIndices(2);
    index_t parentIdx0;
    index_t parentIdx1;

    for (auto edge = dualEdges.rbegin(), last = dualEdges.rend();
            edge != last; ++edge) {
            if (*edge == pair.birth) {
                break;
            }
            boundaryIndices = uf.getBoundaryIndices(*edge);
            parentIdx0 = uf.find(boundaryIndices[0]);
            parentIdx1 = uf.find(boundaryIndices[1]);
            if (parentIdx0 != parentIdx1) {
                    uf.link(parentIdx0, parentIdx1);
            }
    }

    return extractRepresentativeCycle<vector<index_t>>(pair, cgc, uf);
}

void Dimension2::enumerateDualEdges(vector<Cube>& dualEdges, const CubicalGridComplex& cgc) const {
#ifdef RUNTIME
	cout << "enumeration ";
	auto start = high_resolution_clock::now();
#endif 

	dualEdges.reserve(cgc.getNumberOfCubes(2));
	value_t birth;
#ifdef USE_APPARENT_PAIRS
	Cube dualEdge;
	BoundaryEnumerator enumerator(cgc);
	CoboundaryEnumerator coEnumerator(cgc);
#ifdef RUNTIME
	size_t numApparentPairs = 0;
#endif
#endif
#ifdef USE_STABLE_SORT_OR_STABLE_PARTITION
	bool binaryInputs = true;
#endif

	for (index_t x = 0; x < cgc.shape[0]; ++x) {
		for (index_t y = 0; y < cgc.shape[1]; ++y) {
			for (index_t z = 0; z < cgc.shape[2]; ++z) {
				for (uint8_t type = 0; type < 3; ++type) {
					birth = cgc.getBirth(x, y, z, type, 2);
					if (birth < config.threshold) {
#ifdef USE_APPARENT_PAIRS
						dualEdge = Cube(birth, x, y, z, type);
						if (isApparentPair(dualEdge, enumerator, coEnumerator)) {
#ifdef RUNTIME
							++numApparentPairs;
#endif
						}
						else {
							dualEdges.push_back(dualEdge);
#ifdef USE_STABLE_SORT_OR_STABLE_PARTITION
							if (binaryInputs && birth != 0 && birth != 1) binaryInputs = false;
#endif
						}
#else
						dualEdges.push_back(Cube(birth, x, y, z, type));
#ifdef USE_STABLE_SORT_OR_STABLE_PARTITION
						if (binaryInputs && birth != 0 && birth != 1) binaryInputs = false;
#endif
#endif
					}
				}
			}
		}
	}
#ifdef USE_STABLE_SORT_OR_STABLE_PARTITION
	if (binaryInputs) {
		std::stable_partition(dualEdges.begin(), dualEdges.end(), [](Cube &cube) { return cube.birth == 0; });
	} else {
		std::stable_sort(dualEdges.begin(), dualEdges.end(), [](const Cube &cube1, const Cube &cube2) { return cube1.birth < cube2.birth; });
	}
#else
	std::sort(dualEdges.begin(), dualEdges.end(), CubeComparator());
#endif

#ifdef RUNTIME
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<milliseconds>(stop - start);
	cout << duration.count() << " ms, " << dualEdges.size() << " columns to reduce";
#ifdef USE_APPARENT_PAIRS
	cout << ", " << numApparentPairs << " apparent pairs";
#endif
	cout << "; ";
#endif
}


void Dimension2::enumerateDualEdgesComp(vector<Cube>& dualEdges) const {
#ifdef RUNTIME
	cout << "enumeration ";
	auto start = high_resolution_clock::now();
#endif

	dualEdges.reserve(cgcComp.getNumberOfCubes(2));
	value_t birth;
#ifdef USE_STABLE_SORT_OR_STABLE_PARTITION
	bool binaryInputs = true;
#endif
	for (index_t x = 0; x < cgcComp.shape[0]; ++x) {
		for (index_t y = 0; y < cgcComp.shape[1]; ++y) {
			for (index_t z = 0; z < cgcComp.shape[2]; ++z) {
				for (uint8_t type = 0; type < 3; ++type) {
					birth = cgcComp.getBirth(x, y, z, type, 2);
					if (birth < config.threshold) {
						dualEdges.push_back(Cube(birth, x, y, z, type));
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
		std::stable_partition(dualEdges.begin(), dualEdges.end(), [](Cube &cube) { return cube.birth == 0; });
	} else {
		std::stable_sort(dualEdges.begin(), dualEdges.end(), [](const Cube &cube1, const Cube &cube2) { return cube1.birth < cube2.birth; });
	}
#else
	std::sort(dualEdges.begin(), dualEdges.end(), CubeComparator());
#endif
#ifdef RUNTIME
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<milliseconds>(stop - start);
	cout << duration.count() << " ms, " << dualEdges.size() << " columns to reduce; ";
#endif
}


void Dimension2::computeInputAndImagePairs(vector<Cube>& dualEdges, const uint8_t& k) {
#ifdef RUNTIME
	cout << "barcodes ";
	auto start = high_resolution_clock::now();
#endif

	const CubicalGridComplex& cgc = (k == 0) ? cgc0 : cgc1;
	UnionFindDual& uf = (k == 0) ? uf0 : uf1; 
	vector<Pair>& pairs = (k == 0) ? pairs0 : pairs1;
	unordered_map<index_t, Pair>& matchMap = (k==0) ? matchMap0 : matchMap1;
	uf.reset();
	ufComp.reset();
	vector<index_t> boundaryIndices(2);
	index_t parentIdx0;
	index_t parentIdx1;
	index_t birthIdx;
	index_t birthIdxComp;
	value_t birth;
	vector<index_t> birthCoordinates(3);

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
			edge->index = NONE_INDEX;
		}
	}

	auto new_end = remove_if(dualEdges.begin(), dualEdges.end(), [](const Cube& cube){ return cube.index == NONE_INDEX; });
	dualEdges.erase(new_end, dualEdges.end());

#ifdef RUNTIME
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<milliseconds>(stop - start);
	cout << duration.count() << " ms";
#endif
}


void Dimension2::computeCompPairsAndMatch(vector<Cube>& dualEdges, vector<Cube>& ctrImage) {
#ifdef RUNTIME
	cout << "barcode and match ";
	auto start = high_resolution_clock::now();
#endif

	ufComp.reset();
	vector<index_t> boundaryIndices(2);
	index_t parentIdx0;
	index_t parentIdx1;
	index_t birthIdx;
	value_t birth;
	vector<index_t> birthCoordinates(3);
#ifdef USE_APPARENT_PAIRS_COMP
	BoundaryEnumerator enumeratorComp = BoundaryEnumerator(cgcComp);
	CoboundaryEnumerator coEnumeratorComp = CoboundaryEnumerator(cgcComp);
#ifdef RUNTIME
	size_t numApparentPairs = 0;
#endif
#endif

	for (auto edge = dualEdges.rbegin(), last = dualEdges.rend(); edge != last; ++edge) {
		boundaryIndices = ufComp.getBoundaryIndices(*edge);
		parentIdx0 = ufComp.find(boundaryIndices[0]);
		parentIdx1 = ufComp.find(boundaryIndices[1]);
		if (parentIdx0 != parentIdx1) {
			birthIdx = ufComp.link(parentIdx0, parentIdx1);
			birth = ufComp.getBirth(birthIdx);
			if (edge->birth != birth) {
#ifdef COMPUTE_COMPARISON
				birthCoordinates = ufComp.getCoordinates(birthIdx);
				pairsComp.push_back(Pair(*edge, Cube(birth, birthCoordinates[0], birthCoordinates[1], birthCoordinates[2], 0)));
#endif
				auto find0 = matchMap0.find(birthIdx);
				auto find1 = matchMap1.find(birthIdx);
				if (find0 != matchMap0.end() && find1 != matchMap1.end()) {
					matches.push_back(Match(find0->second, find1->second));
					isMatched0.emplace(find0->second.birth.index, true);
					isMatched1.emplace(find1->second.birth.index, true);
				}
			}
			edge->index = NONE_INDEX;
		} else {
#ifdef USE_APPARENT_PAIRS_COMP
			ctrImage.push_back(*edge);
			if (isApparentPair(*edge, enumeratorComp, coEnumeratorComp)) {
				edge->index = NONE_INDEX;
#ifdef RUNTIME
				++numApparentPairs;
#endif
			}
#endif
		}
	}

	auto new_end = std::remove_if(dualEdges.begin(), dualEdges.end(), [](const Cube& cube){ return cube.index == NONE_INDEX; });
	dualEdges.erase(new_end, dualEdges.end());
#ifdef USE_APPARENT_PAIRS_COMP
	reverse(ctrImage.begin(), ctrImage.end());
#endif

#ifdef RUNTIME
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<milliseconds>(stop - start);
	cout << duration.count() << " ms";
#ifdef USE_APPARENT_PAIRS_COMP
	cout << ", " << numApparentPairs << " apparent pairs";
#endif
#endif
}


#if defined(USE_APPARENT_PAIRS) or defined (USE_APPARENT_PAIRS_COMP)
bool Dimension2::isApparentPair(const Cube& dualEdge, BoundaryEnumerator& enumerator, CoboundaryEnumerator& coEnumerator) const {
	enumerator.setBoundaryEnumerator(dualEdge);
	while (enumerator.hasPreviousFace()) {
		if (enumerator.nextFace.birth == dualEdge.birth) {
			coEnumerator.setCoboundaryEnumerator(enumerator.nextFace);
			while (coEnumerator.hasNextCoface()) {
				if (coEnumerator.nextCoface == dualEdge) { return true; } 
				else if (coEnumerator.nextCoface.birth == dualEdge.birth) { return false; }
			}
		}
	}
	return false;
}
#endif