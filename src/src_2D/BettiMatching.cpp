#include "BettiMatching.h"
#include "../data_structures.h"
#include "data_structures.h"
#include "dimension_0.h"
#include "dimension_1.h"

#include <chrono>
#include <iostream>
#include <numeric>
#include <unordered_map>

using namespace dim2;
using namespace std;
using namespace std::chrono;

BettiMatching::BettiMatching(vector<value_t> &&input0, vector<value_t> &&input1,
                             vector<value_t> &&comparison,
                             vector<index_t> &&shape, Config &&_config)
    : cgc0(input0, shape), cgc1(input1, shape), cgcComp(comparison, shape),
      config(_config) {
    pairs0 = vector<vector<Pair>>(2);
    pairs1 = vector<vector<Pair>>(2);
    pairsComp = vector<vector<Pair>>(2);
    matches = vector<vector<Match>>(2);
    isMatched0 = vector<unordered_map<uint64_t, bool>>(2);
    isMatched1 = vector<unordered_map<uint64_t, bool>>(2);
    isMatchedWithIndexComp = vector<unordered_map<uint64_t, size_t>>(2);
    _matched = vector<vector<VoxelMatch>>(2);
    _unmatched0 = vector<vector<VoxelPair>>(2);
    _unmatched1 = vector<vector<VoxelPair>>(2);
}

BettiMatching::BettiMatching(BettiMatching &&other)
    : cgc0(std::move(other.cgc0)), cgc1(std::move(other.cgc1)),
      cgcComp(std::move(other.cgcComp)), config(other.config),
      pairs0(other.pairs0), pairs1(other.pairs1), pairsComp(other.pairsComp),
      matches(other.matches), isMatched0(other.isMatched0),
      isMatched1(other.isMatched1),
      isMatchedWithIndexComp(other.isMatchedWithIndexComp),
      _matched(other.matched), _unmatched0(other.unmatched0),
      _unmatched1(other.unmatched1) {}

void BettiMatching::computeMatching() {
    vector<Cube> ctr0;
    vector<Cube> ctr1;
    vector<Cube> ctrComp;

    {
#ifdef RUNTIME
        cout << "dimension 1:";
        auto start = high_resolution_clock::now();
#endif
        Dimension1 dim1(cgc0, cgc1, cgcComp, config, pairs0[1], pairs1[1],
                        pairsComp[1], matches[1], isMatched0[1], isMatched1[1],
                        isMatchedWithIndexComp[1]);
        dim1.computePairsAndMatch(ctr0, ctr1, ctrComp);
#ifdef RUNTIME
        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<milliseconds>(stop - start);
        cout << endl << "total: " << duration.count() << " ms" << endl << endl;
#endif
    }

    {
#ifdef RUNTIME
        cout << "dimension 0:";
        auto start = high_resolution_clock::now();
#endif
        Dimension0 dim0(cgc0, cgc1, cgcComp, config, pairs0[0], pairs1[0],
                        pairsComp[0], matches[0], isMatched0[0], isMatched1[0],
                        isMatchedWithIndexComp[0]);
        dim0.computePairsAndMatch(ctr0, ctr1, ctrComp);
#ifdef RUNTIME
        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<milliseconds>(stop - start);
        cout << endl << "total: " << duration.count() << " ms" << endl << endl;
#endif
    }
}

void BettiMatching::computeVoxels() {
#ifdef RUNTIME
    cout << "computing voxels ... ";
    auto start = high_resolution_clock::now();
#endif

    for (uint8_t d = 0; d < 2; ++d) {
        for (auto &pair : pairs0[d]) {
            if (!isMatched0[d][pair.birth.index]) {
                _unmatched0[d].push_back(
                    VoxelPair(cgc0.getParentVoxel(pair.birth, d),
                              cgc0.getParentVoxel(pair.death, d + 1)));
            }
        }
        for (auto &pair : pairs1[d]) {
            if (!isMatched1[d][pair.birth.index]) {
                _unmatched1[d].push_back(
                    VoxelPair(cgc1.getParentVoxel(pair.birth, d),
                              cgc1.getParentVoxel(pair.death, d + 1)));
            }
        }
        for (auto &match : matches[d]) {
            _matched[d].push_back(VoxelMatch(
                VoxelPair(cgc0.getParentVoxel(match.pair0.birth, d),
                          cgc0.getParentVoxel(match.pair0.death, d + 1)),
                VoxelPair(cgc1.getParentVoxel(match.pair1.birth, d),
                          cgc1.getParentVoxel(match.pair1.death, d + 1))));
        }
    }

#ifdef RUNTIME
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);
    cout << duration.count() << " ms" << endl << endl;
#endif
}

vector<vector<VoxelPair>> BettiMatching::computePairsInput0() {
    vector<Cube> ctr0;

    Dimension1 dim1(cgc0, cgc1, cgcComp, config, pairs0[1], pairs1[1],
                    pairsComp[1], matches[1], isMatched0[1], isMatched1[1],
                    isMatchedWithIndexComp[1]);
    dim1.computeInput0Pairs(ctr0);

    Dimension0 dim0(cgc0, cgc1, cgcComp, config, pairs0[0], pairs1[0],
                    pairsComp[0], matches[0], isMatched0[0], isMatched1[0],
                    isMatchedWithIndexComp[0]);
    dim0.computeInput0Pairs(ctr0);

    vector<vector<VoxelPair>> voxelPairs(2);

    for (uint8_t d = 0; d < 2; ++d) {
        for (auto &pair : pairs0[d]) {
            voxelPairs[d].push_back(
                VoxelPair(cgc0.getParentVoxel(pair.birth, d),
                          cgc0.getParentVoxel(pair.death, d + 1)));
        }
    }
    return voxelPairs;
}

void BettiMatching::printResult() {
    index_t count;

    cout << "------------------------------------------------------------------"
            "---------------------------------------------"
         << endl;
    cout << "Input 0:" << endl << endl;
    if (cgc0.shape[0] < 10 && cgc0.shape[1] < 10) {
        cgc0.printImage();
    }
    for (uint8_t d = 0; d < 2; ++d) {
        cout << "dim " << unsigned(d) << ": ";
        count = pairs0[d].size();
        if (0 < count && count < 10) {
            cout << endl;
            for (auto &pair : pairs0[d]) {
                pair.print();
                cout << endl;
            }
        } else {
            cout << count << endl;
        }
    }
    cout << "------------------------------------------------------------------"
            "---------------------------------------------"
         << endl;
    cout << "Input 1" << endl << endl;
    if (cgc1.shape[0] < 10 && cgc1.shape[1] < 10) {
        cgc1.printImage();
    }
    for (uint8_t d = 0; d < 2; ++d) {
        cout << "dim " << unsigned(d) << ": ";
        count = pairs1[d].size();
        if (0 < count && count < 10) {
            cout << endl;
            for (auto &pair : pairs1[d]) {
                pair.print();
                cout << endl;
            }
        } else {
            cout << count << endl;
        }
    }

#ifdef COMPUTE_COMPARISON
    cout << "------------------------------------------------------------------"
            "---------------------------------------------"
         << endl;
    cout << "Comparison" << endl << endl;
    if (cgcComp.shape[0] < 10 && cgcComp.shape[1] < 10) {
        cgcComp.printImage();
    }
    for (uint8_t d = 0; d < 2; ++d) {
        cout << "dim " << unsigned(d) << ": ";
        count = pairsComp[d].size();
        if (0 < count && count < 10) {
            cout << endl;
            for (auto &pair : pairsComp[d]) {
                pair.print();
                cout << endl;
            }
        } else {
            cout << "number of pairs: " << count << endl;
        }
    }
#endif

    cout << "------------------------------------------------------------------"
            "---------------------------------------------"
         << endl;
    cout << "Betti Matching:" << endl << endl;
    cout << "matched: " << endl;
    for (uint8_t d = 0; d < 2; ++d) {
        cout << "dim " << unsigned(d) << ": ";
        count = matches[d].size();
        if (0 < count && count < 10) {
            cout << endl;
            std::vector<size_t> range(count);
            std::iota(range.begin(), range.end(), 0);
            vector<RepresentativeCycle> matchedCycles0 =
                std::get<0>(computeRepresentativeCycles(
                    0, d, range, {std::vector<size_t>()}));
            vector<RepresentativeCycle> matchedCycles1 =
                std::get<0>(computeRepresentativeCycles(
                    1, d, range, {std::vector<size_t>()}));
            for (size_t i = 0; i < matches[d].size(); i++) {
                matches[d][i].print();
                _matched[d][i].print();
                if (cgc0.shape[0] < 10 && cgc0.shape[1] < 10 &&
                    cgc0.shape[2] < 10) {
                    cgc0.printRepresentativeCycle(matchedCycles0[i]);
                    cout << endl;
                    cgc1.printRepresentativeCycle(matchedCycles1[i]);
                    cout << endl;
                }
            }
        } else {
            cout << count << endl;
        }
    }

    size_t counter;
    cout << endl << "unmatched in Input 0" << endl;
    for (uint8_t d = 0; d < 2; ++d) {
        cout << "dim " << unsigned(d) << ": ";
        count = 0;
        for (auto &pair : pairs0[d]) {
            if (!isMatched0[d][pair.birth.index]) {
                ++count;
            }
        }
        if (0 < count && count < 10) {
            cout << endl;
            std::vector<size_t> range(count);
            std::iota(range.begin(), range.end(), 0);
            vector<RepresentativeCycle> unmatchedCycles =
                std::get<1>(computeRepresentativeCycles(
                    0, d, {std::vector<size_t>()}, range));
            counter = 0;
            for (size_t i = 0; i < pairs0[d].size(); ++i) {
                if (!isMatched0[d][pairs0[d][i].birth.index]) {
                    pairs0[d][i].print();
                    cout << endl;
                    _unmatched0[d][counter].print();
                    cout << endl;
                    cgc0.printRepresentativeCycle(unmatchedCycles[counter]);
                    cout << endl;
                    ++counter;
                    if (counter == count) {
                        break;
                    }
                }
            }
        } else {
            cout << count << endl;
        }
    }

    cout << endl << "unmatched in Input 1" << endl;
    for (uint8_t d = 0; d < 2; ++d) {
        cout << "dim " << unsigned(d) << ": ";
        count = 0;
        for (auto &pair : pairs1[d]) {
            if (!isMatched1[d][pair.birth.index]) {
                ++count;
            }
        }
        if (0 < count && count < 10) {
            cout << endl;
            std::vector<size_t> range(count);
            std::iota(range.begin(), range.end(), 0);
            vector<RepresentativeCycle> unmatchedCycles =
                std::get<1>(computeRepresentativeCycles(
                    1, d, {std::vector<size_t>()}, range));
            counter = 0;
            for (size_t i = 0; i < pairs1[d].size(); ++i) {
                if (!isMatched1[d][pairs1[d][i].birth.index]) {
                    pairs1[d][i].print();
                    cout << endl;
                    _unmatched1[d][counter].print();
                    cout << endl;
                    cgc1.printRepresentativeCycle(unmatchedCycles[counter]);
                    cout << endl;
                    ++counter;
                    if (counter == count) {
                        break;
                    }
                }
            }
        } else {
            cout << count << endl;
        }
    }
}

tuple<vector<dim2::RepresentativeCycle>, vector<dim2::RepresentativeCycle>>
BettiMatching::computeRepresentativeCycles(
    const int input, const int dim,
    const optional<vector<size_t>> &matchedPairsIndices,
    const optional<vector<size_t>> &unmatchedPairsIndices) {
    if (dim < 0 || dim > 1) {
        throw runtime_error("Invalid value for dim");
    }

    // Assemble the list of requested pairs: First the matched pairs (all if
    // empty optional was passed, then the unmatched pairs (all if empty
    // optional was passed)
    vector<std::reference_wrapper<Pair>> requestedPairs;
    if (input == 0 || input == 1) {
        requestedPairs = assembleRequestedPairs(
            matchedPairsIndices, unmatchedPairsIndices,
            (input == 0 ? pairs0 : pairs1)[dim],
            (input == 0 ? isMatched0 : isMatched1)[dim], matches[dim], input);
    } else if (input == 2) {
        requestedPairs = assembleRequestedComparisonPairs<Pair>(
            matchedPairsIndices, pairsComp[dim], isMatchedWithIndexComp[dim]);
    } else {
        throw runtime_error("Invalid value for input");
    }

    // Hand over the representative cycle computation to the respective
    // dimension
    vector<RepresentativeCycle> representativeCycles;
    switch (dim) {
    case 0: {
        Dimension0 dim0(cgc0, cgc1, cgcComp, config, pairs0[0], pairs1[0],
                        pairsComp[0], matches[0], isMatched0[0], isMatched1[0],
                        isMatchedWithIndexComp[0]);
        representativeCycles =
            dim0.computeRepresentativeCycles(input, requestedPairs);
        break;
    }
    case 1: {
        Dimension1 dim1(cgc0, cgc1, cgcComp, config, pairs0[1], pairs1[1],
                        pairsComp[1], matches[1], isMatched0[1], isMatched1[1],
                        isMatchedWithIndexComp[1]);
        representativeCycles =
            dim1.computeRepresentativeCycles(input, requestedPairs);
        break;
    }
    }

    // Split the returned representative cycles into a matched and unmatched
    // portion
    vector<RepresentativeCycle> matchedRepresentativeCycles;
    vector<RepresentativeCycle> unmatchedRepresentativeCycles;

    int numMatchedRequested = matchedPairsIndices.has_value()
                                  ? matchedPairsIndices->size()
                                  : matches[dim].size();
    matchedRepresentativeCycles.insert(
        matchedRepresentativeCycles.end(),
        std::make_move_iterator(representativeCycles.begin()),
        std::make_move_iterator(representativeCycles.begin()) +
            numMatchedRequested);
    unmatchedRepresentativeCycles.insert(
        unmatchedRepresentativeCycles.end(),
        std::make_move_iterator(representativeCycles.begin() +
                                numMatchedRequested),
        std::make_move_iterator(representativeCycles.end()));

    return {matchedRepresentativeCycles, unmatchedRepresentativeCycles};
}
