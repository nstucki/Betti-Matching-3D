#include "dimension_2.h"
#include "BettiMatching.h"
#include "data_structures.h"
#include "enumerators.h"

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <functional>
#include <iostream>
#include <map>
#include <set>
#include <stdexcept>
#include <tuple>
#include <unordered_map>

using namespace dim3;
using namespace std;
using namespace std::chrono;

Dimension2::Dimension2(const CubicalGridComplex &_cgc0,
                       const CubicalGridComplex &_cgc1,
                       const CubicalGridComplex &_cgcComp,
                       const Config &_config, vector<Pair> &_pairs0,
                       vector<Pair> &_pairs1, vector<Pair> &_pairsComp,
                       vector<Match> &_matches,
                       unordered_map<uint64_t, bool> &_isMatched0,
                       unordered_map<uint64_t, bool> &_isMatched1,
                       unordered_map<uint64_t, size_t> &_isMatchedWithIndexComp)
    : cgc0(_cgc0), cgc1(_cgc1), cgcComp(_cgcComp), config(_config),
      pairs0(_pairs0), pairs1(_pairs1), pairsComp(_pairsComp),
      matches(_matches), isMatched0(_isMatched0), isMatched1(_isMatched1),
      isMatchedWithIndexComp(_isMatchedWithIndexComp), uf0(UnionFindDual(cgc0)),
      uf1(UnionFindDual(cgc1)), ufComp(UnionFindDual(cgcComp)) {}

void Dimension2::computePairsAndMatch(vector<Cube> &ctr0, vector<Cube> &ctr1,
                                      vector<Cube> &ctrComp,
                                      vector<Cube> &ctrImage) {
    size_t actualDim = 3;
    for (index_t s : cgc0.shape) {
        if (s == 1) {
            --actualDim;
        }
    }
    bool needToCompute = (actualDim == 3);
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

void Dimension2::computeInput0Pairs(vector<Cube> &ctr0) {
    size_t actualDim = 3;
    for (index_t s : cgc0.shape) {
        if (s == 1) {
            --actualDim;
        }
    }
    bool needToCompute = (actualDim == 3);
    if (!needToCompute) {
        enumerateDualEdges(ctr0, cgc0);
    } else {
        enumerateDualEdges(ctr0, cgc0);
        computeInputAndImagePairs(ctr0, 0);
    }
}

void Dimension2::enumerateDualEdges(vector<Cube> &dualEdges,
                                    const CubicalGridComplex &cgc) const {
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
                        if (isApparentPair(dualEdge, enumerator,
                                           coEnumerator)) {
#ifdef RUNTIME
                            ++numApparentPairs;
#endif
                        } else {
                            dualEdges.push_back(dualEdge);
#ifdef USE_STABLE_SORT_OR_STABLE_PARTITION
                            if (binaryInputs && birth != 0 && birth != 1)
                                binaryInputs = false;
#endif
                        }
#else
                        dualEdges.push_back(Cube(birth, x, y, z, type));
#ifdef USE_STABLE_SORT_OR_STABLE_PARTITION
                        if (binaryInputs && birth != 0 && birth != 1)
                            binaryInputs = false;
#endif
#endif
                    }
                }
            }
        }
    }
#ifdef USE_STABLE_SORT_OR_STABLE_PARTITION
    if (binaryInputs) {
        std::stable_partition(dualEdges.begin(), dualEdges.end(),
                              [](Cube &cube) { return cube.birth == 0; });
    } else {
        std::stable_sort(dualEdges.begin(), dualEdges.end(),
                         [](const Cube &cube1, const Cube &cube2) {
                             return cube1.birth < cube2.birth;
                         });
    }
#else
    std::sort(dualEdges.begin(), dualEdges.end(), CubeComparator());
#endif

#ifdef RUNTIME
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);
    cout << duration.count() << " ms, " << dualEdges.size()
         << " columns to reduce";
#ifdef USE_APPARENT_PAIRS
    cout << ", " << numApparentPairs << " apparent pairs";
#endif
    cout << "; ";
#endif
}

void Dimension2::enumerateDualEdgesComp(vector<Cube> &dualEdges) const {
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
                        if (binaryInputs && birth != 0 && birth != 1)
                            binaryInputs = false;
#endif
                    }
                }
            }
        }
    }
#ifdef USE_STABLE_SORT_OR_STABLE_PARTITION
    if (binaryInputs) {
        std::stable_partition(dualEdges.begin(), dualEdges.end(),
                              [](Cube &cube) { return cube.birth == 0; });
    } else {
        std::stable_sort(dualEdges.begin(), dualEdges.end(),
                         [](const Cube &cube1, const Cube &cube2) {
                             return cube1.birth < cube2.birth;
                         });
    }
#else
    std::sort(dualEdges.begin(), dualEdges.end(), CubeComparator());
#endif
#ifdef RUNTIME
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);
    cout << duration.count() << " ms, " << dualEdges.size()
         << " columns to reduce; ";
#endif
}

void Dimension2::computeInputAndImagePairs(vector<Cube> &dualEdges,
                                           const uint8_t &k) {
#ifdef RUNTIME
    cout << "barcodes ";
    auto start = high_resolution_clock::now();
#endif

    UnionFindDual &uf = (k == 0) ? uf0 : uf1;
    vector<Pair> &pairs = (k == 0) ? pairs0 : pairs1;
    unordered_map<index_t, Pair> &matchMap = (k == 0) ? matchMap0 : matchMap1;
    uf.reset();
    ufComp.reset();
    vector<index_t> boundaryIndices(2);
    index_t parentIdx0;
    index_t parentIdx1;
    index_t birthIdx;
    index_t birthIdxComp;
    value_t birth;
    Coordinate birthCoordinates;

    for (auto edge = dualEdges.rbegin(), last = dualEdges.rend(); edge != last;
         ++edge) {
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
                pairs.push_back(
                    Pair(*edge, Cube(birth, std::get<0>(birthCoordinates),
                                     std::get<1>(birthCoordinates),
                                     std::get<2>(birthCoordinates), 0)));
                matchMap.emplace(birthIdxComp, pairs.back());
            }
            edge->index = NONE_INDEX;
        }
    }

    auto new_end =
        remove_if(dualEdges.begin(), dualEdges.end(),
                  [](const Cube &cube) { return cube.index == NONE_INDEX; });
    dualEdges.erase(new_end, dualEdges.end());

#ifdef RUNTIME
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);
    cout << duration.count() << " ms";
#endif
}

void Dimension2::computeCompPairsAndMatch(vector<Cube> &dualEdges,
                                          vector<Cube> &ctrImage) {
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
#ifdef USE_APPARENT_PAIRS_COMP
    BoundaryEnumerator enumeratorComp = BoundaryEnumerator(cgcComp);
    CoboundaryEnumerator coEnumeratorComp = CoboundaryEnumerator(cgcComp);
#ifdef RUNTIME
    size_t numApparentPairs = 0;
#endif
#endif

    for (auto edge = dualEdges.rbegin(), last = dualEdges.rend(); edge != last;
         ++edge) {
        boundaryIndices = ufComp.getBoundaryIndices(*edge);
        parentIdx0 = ufComp.find(boundaryIndices[0]);
        parentIdx1 = ufComp.find(boundaryIndices[1]);
        if (parentIdx0 != parentIdx1) {
            birthIdx = ufComp.link(parentIdx0, parentIdx1);
            birth = ufComp.getBirth(birthIdx);
            if (edge->birth != birth) {
#ifdef COMPUTE_COMPARISON
                auto birthCoordinates = ufComp.getCoordinates(birthIdx);
                pairsComp.push_back(
                    Pair(*edge, Cube(birth, std::get<0>(birthCoordinates),
                                     std::get<1>(birthCoordinates),
                                     std::get<2>(birthCoordinates), 0)));
#endif
                auto find0 = matchMap0.find(birthIdx);
                auto find1 = matchMap1.find(birthIdx);
                if (find0 != matchMap0.end() && find1 != matchMap1.end()) {
                    matches.push_back(Match(find0->second, find1->second));
                    isMatched0.emplace(find0->second.birth.index, true);
                    isMatched1.emplace(find1->second.birth.index, true);
#ifdef COMPUTE_COMPARISON
                    isMatchedWithIndexComp.emplace(edge->index,
                                                   matches.size() - 1);
#endif
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

    auto new_end = std::remove_if(
        dualEdges.begin(), dualEdges.end(),
        [](const Cube &cube) { return cube.index == NONE_INDEX; });
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

#if defined(USE_APPARENT_PAIRS) or defined(USE_APPARENT_PAIRS_COMP)
bool Dimension2::isApparentPair(const Cube &dualEdge,
                                BoundaryEnumerator &enumerator,
                                CoboundaryEnumerator &coEnumerator) const {
    enumerator.setBoundaryEnumerator(dualEdge);
    while (enumerator.hasPreviousFace()) {
        if (enumerator.nextFace.birth == dualEdge.birth) {
            coEnumerator.setCoboundaryEnumerator(enumerator.nextFace);
            while (coEnumerator.hasNextCoface()) {
                if (coEnumerator.nextCoface == dualEdge) {
                    return true;
                } else if (coEnumerator.nextCoface.birth == dualEdge.birth) {
                    return false;
                }
            }
        }
    }
    return false;
}
#endif

RepresentativeCycle
extractRepresentativeCycle(const Pair &pair, const CubicalGridComplex &cgc,
                           vector<dim3::Coordinate> &dualConnectedComponent) {
    map<Coordinate, size_t> boundaryVertices;
    for (const Coordinate &c : dualConnectedComponent) {
        for (uint8_t x = 0; x < 2; ++x) {
            for (uint8_t y = 0; y < 2; ++y) {
                for (uint8_t z = 0; z < 2; ++z) {
                    auto vertex =
                        Coordinate{std::get<0>(c) + x, std::get<1>(c) + y,
                                   std::get<2>(c) + z};
                    auto entry = boundaryVertices.find(vertex);
                    if (entry == boundaryVertices.end()) {
                        boundaryVertices[vertex] = 1;
                    } else {
                        entry->second++;
                    }
                }
            }
        }
    }

    vector<Coordinate> reprCycle;
    for (auto &entry : boundaryVertices) {
        if (entry.second < 8) {
            reprCycle.push_back(entry.first);
        }
    }
    reprCycle.push_back(cgc.getParentVoxel(pair.death, 3));

    return reprCycle;
}

vector<dim3::RepresentativeCycle> Dimension2::computeRepresentativeCycles(
    const int input,
    const std::vector<std::reference_wrapper<Pair>> &requestedPairs) {
    if (requestedPairs.size() == 0) {
        return {};
    }

    const CubicalGridComplex &cgc = (input == 0)   ? cgc0
                                    : (input == 1) ? cgc1
                                                   : cgcComp;
    UnionFindDual uf(cgc);

    vector<Cube> dualEdges;
    if (input == 0 || input == 1) {
        enumerateDualEdges(dualEdges, cgc);
    } else {
        enumerateDualEdgesComp(dualEdges);
    }

    unordered_map<uint64_t, RepresentativeCycle> cyclesByBirth;

    // Initialize singleton dual connected components
    vector<RepresentativeCycle> dualConnectedComponentsByBirthIdx(
        cgc.shape[0] * cgc.shape[1] * cgc.shape[2]);
    for (int birthIdx = 0; birthIdx < dualConnectedComponentsByBirthIdx.size();
         birthIdx++) {
        dualConnectedComponentsByBirthIdx[birthIdx].emplace_back(
            uf.getCoordinates(birthIdx));
    }

    // Gather the birth indices of requested pairs for use in the loop below
    // (the pairs are not necessarily ordered in the same order as they are
    // found, hence we need a set)
    unordered_map<uint64_t, std::reference_wrapper<Pair>> requestedPairsByBirth;
    for (auto &pair : requestedPairs) {
        requestedPairsByBirth.emplace(pair.get().birth.index, pair);
    }

    // Run the union-find algorithm on the dual edges and collect all requested
    // representative cycles on the way
    for (auto edge = dualEdges.rbegin(), last = dualEdges.rend();
         edge != last && cyclesByBirth.size() != requestedPairs.size();
         edge++) {
        // Compute the representative cycle if it was requested
        vector<index_t> boundaryIndices = uf.getBoundaryIndices(*edge);
        index_t parentIdx0 = uf.find(boundaryIndices[0]);
        index_t parentIdx1 = uf.find(boundaryIndices[1]);
        if (parentIdx0 != parentIdx1) {
            auto olderBirthIdx = uf.link(parentIdx0, parentIdx1);
            auto youngerBirthIdx =
                (parentIdx1 == olderBirthIdx) ? parentIdx0 : parentIdx1;
            dualConnectedComponentsByBirthIdx[youngerBirthIdx].insert(
                dualConnectedComponentsByBirthIdx[youngerBirthIdx].end(),
                dualConnectedComponentsByBirthIdx[olderBirthIdx].begin(),
                dualConnectedComponentsByBirthIdx[olderBirthIdx].end());

            auto maybeRequestedPair = requestedPairsByBirth.find(edge->index);
            if (maybeRequestedPair != requestedPairsByBirth.end()) {
                cyclesByBirth[edge->index] = extractRepresentativeCycle(
                    maybeRequestedPair->second, cgc,
                    dualConnectedComponentsByBirthIdx[olderBirthIdx]);
            }
            // The deceased dual connected component has either been converted
            // to its bordering 2-cycle, or the 2-cycle has not been requested.
            // Either way, we can delete it and save memory.
            dualConnectedComponentsByBirthIdx[olderBirthIdx].clear();
        }
    }

    if (cyclesByBirth.size() != requestedPairs.size()) {
        throw runtime_error(
            "Not all requested representative cycles were found");
    }

    vector<RepresentativeCycle> representativeCycles;
    representativeCycles.reserve(requestedPairs.size());
    for (auto &pair : requestedPairs) {
        representativeCycles.emplace_back(
            std::move(cyclesByBirth[pair.get().birth.index]));
    }
    return representativeCycles;
}
