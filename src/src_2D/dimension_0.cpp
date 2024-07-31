#include "dimension_0.h"
#include "data_structures.h"

#include <algorithm>
#include <chrono>
#include <iostream>

using namespace dim2;
using namespace std::chrono;

Dimension0::Dimension0(const CubicalGridComplex &_cgc0,
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
      isMatchedWithIndexComp(_isMatchedWithIndexComp), uf0(UnionFind(cgc0)),
      uf1(UnionFind(cgc1)), ufComp(UnionFind(cgcComp)) {}

void Dimension0::computePairsAndMatch(vector<Cube> &ctr0, vector<Cube> &ctr1,
                                      vector<Cube> &ctrComp) {
#ifdef RUNTIME
    cout << endl << "input 0: ";
#endif
    computePairs(ctr0, 0);

#ifdef RUNTIME
    cout << endl << "input 1: ";
#endif
    computePairs(ctr1, 1);

#ifdef RUNTIME
    cout << endl << "comparison & image 0 & image 1 & matching: ";
#endif
    uf0.reset();
    uf1.reset();
    computeImagePairsAndMatch(ctrComp);
}

void Dimension0::computeInput0Pairs(vector<Cube> &ctr0) {
    computePairs(ctr0, 0);
}

void Dimension0::computePairs(vector<Cube> &edges, uint8_t k) {
#ifdef RUNTIME
    cout << "barcode ";
    auto start = high_resolution_clock::now();
#endif
    UnionFind &uf = (k == 0) ? uf0 : uf1;
    vector<Pair> &pairs = (k == 0) ? pairs0 : pairs1;
    unordered_map<index_t, Pair> &matchMap = (k == 0) ? matchMap0 : matchMap1;

    vector<index_t> boundaryIndices(2);
    index_t parentIdx0;
    index_t parentIdx1;
    index_t birthIdx;
    value_t birth;
    dim2::Coordinate birthCoordinates;
    for (Cube &edge : edges) {
        boundaryIndices = uf.getBoundaryIndices(edge);
        parentIdx0 = uf.find(boundaryIndices[0]);
        parentIdx1 = uf.find(boundaryIndices[1]);
        if (parentIdx0 != parentIdx1) {
            birthIdx = uf.link(parentIdx0, parentIdx1);
            birth = uf.getBirth(birthIdx);
            if (birth != edge.birth) {
                birthCoordinates = uf.getCoordinates(birthIdx);
                pairs.push_back(Pair(Cube(birth, std::get<0>(birthCoordinates),
                                          std::get<1>(birthCoordinates), 0),
                                     edge));
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

void Dimension0::computeImagePairsAndMatch(vector<Cube> &edges) {
#ifdef RUNTIME
    cout << "barcodes & matching ";
    auto start = high_resolution_clock::now();
#endif
    vector<index_t> boundaryIndices(2);
    index_t parentIdx0;
    index_t parentIdx1;
    index_t birthIdx0;
    index_t birthIdx1;
    index_t birthIdxComp;
    value_t birth;
    dim2::Coordinate birthCoordinates;
    for (Cube &edge : edges) {
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
                pairsComp.push_back(
                    Pair(Cube(birth, std::get<0>(birthCoordinates),
                              std::get<1>(birthCoordinates), 0),
                         edge));
#endif
                auto find0 = matchMap0.find(birthIdx0);
                auto find1 = matchMap1.find(birthIdx1);
                if (find0 != matchMap0.end() && find1 != matchMap1.end()) {
                    matches.push_back(Match(find0->second, find1->second));
                    isMatched0.emplace(find0->second.birth.index, true);
                    isMatched1.emplace(find1->second.birth.index, true);
#ifdef COMPUTE_COMPARISON
                    isMatchedWithIndexComp.emplace(pairsComp.back().birth.index,
                                                   matches.size() - 1);
#endif
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

void Dimension0::enumerateEdges(vector<Cube> &edges,
                                const CubicalGridComplex &cgc) const {
#ifdef RUNTIME
    cout << "enumeration ";
    auto start = high_resolution_clock::now();
#endif

    edges.reserve(cgc.getNumberOfCubes(1));
    value_t birth;
    for (index_t x = 0; x < cgc.shape[0]; ++x) {
        for (index_t y = 0; y < cgc.shape[1]; ++y) {
            for (uint8_t type = 0; type < 2; ++type) {
                birth = cgc.getBirth(x, y, type, 1);
                if (birth < config.threshold) {
                    edges.push_back(Cube(birth, x, y, type));
                }
            }
        }
    }

    sort(edges.begin(), edges.end(), CubeComparator());

#ifdef RUNTIME
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);
    cout << duration.count() << " ms, ";
#endif
}

vector<dim2::RepresentativeCycle> Dimension0::computeRepresentativeCycles(
    const int input,
    const std::vector<std::reference_wrapper<Pair>> &requestedPairs) {
    if (requestedPairs.size() == 0) {
        return {};
    }

    const CubicalGridComplex &cgc = (input == 0)   ? cgc0
                                    : (input == 1) ? cgc1
                                                   : cgcComp;
    UnionFind uf(cgc);

    // Map from cube indices to union find indices (to match cycles with
    // persistence pairs)
    unordered_map<uint64_t, index_t> unionFindIdxByDeath;

    // Gather the birth indices of requested pairs for use in the loop below
    // (the pairs are not necessarily ordered in the same order as they are
    // found, hence we need a set)
    unordered_map<uint64_t, std::reference_wrapper<Pair>> requestedPairsByDeath;
    for (auto &pair : requestedPairs) {
        requestedPairsByDeath.emplace(pair.get().death.index, pair);
    }

    // Initialize singleton representative cycles
    vector<RepresentativeCycle> cycleByBirthIdx(cgc.shape[0] * cgc.shape[1]);
    for (int birthIdx = 0; birthIdx < cycleByBirthIdx.size(); birthIdx++) {
        cycleByBirthIdx[birthIdx].emplace_back(uf.getCoordinates(birthIdx));
    }

    vector<Cube> edges;
    enumerateEdges(edges, cgc);

    int foundRequestedPairs = 0;
    // Follow the union-find algorithm:
    for (Cube &edge : edges) {
        vector<index_t> boundaryIndices = uf.getBoundaryIndices(edge);
        index_t parentIdx0 = uf.find(boundaryIndices[0]);
        index_t parentIdx1 = uf.find(boundaryIndices[1]);
        // When merging a younger component into an older one, extend the older
        // component's cycle by the younger component's cycle
        if (parentIdx0 != parentIdx1) {
            auto youngerBirthIdx = uf.link(parentIdx0, parentIdx1);
            auto olderBirthIdx =
                (parentIdx1 == youngerBirthIdx) ? parentIdx0 : parentIdx1;
            cycleByBirthIdx[olderBirthIdx].insert(
                cycleByBirthIdx[olderBirthIdx].end(),
                cycleByBirthIdx[youngerBirthIdx].begin(),
                cycleByBirthIdx[youngerBirthIdx].end());
            unionFindIdxByDeath[edge.index] = youngerBirthIdx;

            auto maybeRequestedPair = requestedPairsByDeath.find(edge.index);
            // If the deceased component does not correspond to a pair we'd like
            // to save, delete it to save memory:
            if (maybeRequestedPair == requestedPairsByDeath.end()) {
                cycleByBirthIdx[youngerBirthIdx].clear();
            }
            // Else, increment the counter of found requested pairs, and
            // terminate if we found all
            else {
                foundRequestedPairs++;
                if (foundRequestedPairs == requestedPairs.size()) {
                    break;
                }
            }
        }
    }

    if (foundRequestedPairs != requestedPairs.size()) {
        throw runtime_error(
            "Not all requested representative cycles were found");
    }

    vector<RepresentativeCycle> representativeCycles;
    representativeCycles.reserve(requestedPairs.size());
    for (auto &pair : requestedPairs) {
        auto unionFindIdx = unionFindIdxByDeath.find(pair.get().death.index);
        if (unionFindIdx == unionFindIdxByDeath.end()) {
            throw runtime_error("Union find index for pair cannot be found");
        }
        auto &representativeCycle = cycleByBirthIdx[unionFindIdx->second];
        representativeCycle.push_back(cgc.getParentVoxel(pair.get().death, 1));
        representativeCycles.emplace_back(std::move(representativeCycle));
    }
    return representativeCycles;
}
