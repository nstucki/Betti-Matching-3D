#include "dimension_1.h"
#include "data_structures.h"

#include <algorithm>
#include <chrono>
#include <iostream>
#include <map>
#include <set>

using namespace dim2;
using namespace std;
using namespace std::chrono;

Dimension1::Dimension1(const CubicalGridComplex &_cgc0,
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

void Dimension1::computePairsAndMatch(vector<Cube> &ctr0, vector<Cube> &ctr1,
                                      vector<Cube> &ctrComp) {
#ifdef RUNTIME
    cout << endl << "input & image 0: ";
#endif
    enumerateDualEdges(ctr0, cgc0);
    computeInputAndImagePairs(ctr0, 0);

#ifdef RUNTIME
    cout << endl << "input & image 1: ";
#endif
    enumerateDualEdges(ctr1, cgc1);
    ufComp.reset();
    computeInputAndImagePairs(ctr1, 1);

#ifdef RUNTIME
    cout << endl << "comparison & matching: ";
#endif
    enumerateDualEdges(ctrComp, cgcComp);
    ufComp.reset();
    computeCompPairsAndMatch(ctrComp);
}

void Dimension1::computeInput0Pairs(vector<Cube> &ctr0) {
    enumerateDualEdges(ctr0, cgc0);
    computeInputAndImagePairs(ctr0, 0);
}

void Dimension1::enumerateDualEdges(vector<Cube> &dualEdges,
                                    const CubicalGridComplex &cgc) const {
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
                if (birth < config.threshold) {
                    dualEdges.push_back(Cube(birth, x, y, type));
                }
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

void Dimension1::computeInputAndImagePairs(vector<Cube> &dualEdges,
                                           const uint8_t &k) {
#ifdef RUNTIME
    cout << "barcodes ";
    auto start = high_resolution_clock::now();
#endif

    UnionFindDual &uf = (k == 0) ? uf0 : uf1;
    vector<Pair> &pairs = (k == 0) ? pairs0 : pairs1;
    unordered_map<index_t, Pair> &matchMap = (k == 0) ? matchMap0 : matchMap1;
    vector<index_t> boundaryIndices;
    index_t parentIdx0;
    index_t parentIdx1;
    index_t birthIdx;
    index_t birthIdxComp;
    value_t birth;
    dim2::Coordinate birthCoordinates;
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
                                     std::get<1>(birthCoordinates), 0)));
                matchMap.emplace(birthIdxComp, pairs.back());
            }
#ifdef USE_CLEARING_DIM0
            edge->index = NONE_INDEX;
#endif
        }
    }

#ifdef USE_CLEARING_DIM0
    auto new_end =
        remove_if(dualEdges.begin(), dualEdges.end(),
                  [](const Cube &cube) { return cube.index == NONE_INDEX; });
    dualEdges.erase(new_end, dualEdges.end());
#endif

#ifdef RUNTIME
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);
    cout << duration.count() << " ms";
#endif
}

void Dimension1::computeCompPairsAndMatch(vector<Cube> &dualEdges) {
#ifdef RUNTIME
    cout << "barcode and matching ";
    auto start = high_resolution_clock::now();
#endif

    vector<index_t> boundaryIndices;
    index_t parentIdx0;
    index_t parentIdx1;
    index_t birthIdx;
    value_t birth;
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
                                     std::get<1>(birthCoordinates), 0)));
#endif
                auto find0 = matchMap0.find(birthIdx);
                auto find1 = matchMap1.find(birthIdx);
                if (find0 != matchMap0.end() && find1 != matchMap1.end()) {
                    matches.push_back(Match(find0->second, find1->second));
                    isMatched0.emplace(find0->second.birth.index, true);
                    isMatched1.emplace(find1->second.birth.index, true);
                    isMatchedWithIndexComp.emplace(edge->index,
                                                   matches.size() - 1);
                }
            }
#ifdef USE_CLEARING_DIM0
            edge->index = NONE_INDEX;
#endif
        }
    }

#ifdef USE_CLEARING_DIM0
    auto new_end = std::remove_if(
        dualEdges.begin(), dualEdges.end(),
        [](const Cube &cube) { return cube.index == NONE_INDEX; });
    dualEdges.erase(new_end, dualEdges.end());
#endif

#ifdef RUNTIME
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);
    cout << duration.count() << " ms";
#endif
}

RepresentativeCycle
extractRepresentativeCycle(const Pair &pair, const CubicalGridComplex &cgc,
                           vector<dim2::Coordinate> &dualConnectedComponent) {
    map<Coordinate, size_t> boundaryVertices;
    for (const Coordinate &c : dualConnectedComponent) {
        for (uint8_t x = 0; x < 2; ++x) {
            for (uint8_t y = 0; y < 2; ++y) {
                auto vertex =
                    Coordinate{std::get<0>(c) + x, std::get<1>(c) + y};
                auto entry = boundaryVertices.find(vertex);
                if (entry == boundaryVertices.end()) {
                    boundaryVertices[vertex] = 1;
                } else {
                    entry->second++;
                }
            }
        }
    }

    vector<Coordinate> reprCycle;
    for (auto &entry : boundaryVertices) {
        if (entry.second < 4) {
            reprCycle.push_back(entry.first);
        }
    }
    reprCycle.push_back(cgc.getParentVoxel(pair.death, 2));

    return reprCycle;
}

vector<dim2::RepresentativeCycle> Dimension1::computeRepresentativeCycles(
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
    enumerateDualEdges(dualEdges, cgc);

    unordered_map<uint64_t, RepresentativeCycle> cyclesByBirth;

    // Initialize singleton dual connected components
    vector<RepresentativeCycle> dualConnectedComponentsByBirthIdx(cgc.shape[0] *
                                                                  cgc.shape[1]);
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
