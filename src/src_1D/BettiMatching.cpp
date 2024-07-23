#include "BettiMatching.h"
#include "../data_structures.h"
#include "data_structures.h"
#include "dimension_0.h"

#include <chrono>
#include <iostream>

using namespace dim1;
using namespace std;
using namespace std::chrono;

BettiMatching::BettiMatching(vector<value_t> &&input0, vector<value_t> &&input1,
                             vector<value_t> &&comparison,
                             vector<index_t> &&shape, Config &&_config)
    : cgc0(input0, shape), cgc1(input1, shape), cgcComp(comparison, shape),
      config(_config) {
    pairs0 = vector<vector<Pair>>(1);
    pairs1 = vector<vector<Pair>>(1);
    pairsComp = vector<vector<Pair>>(1);
    matches = vector<vector<Match>>(1);
    isMatched0 = vector<unordered_map<uint64_t, bool>>(1);
    isMatched1 = vector<unordered_map<uint64_t, bool>>(1);
}

BettiMatching::BettiMatching(BettiMatching &&other)
    : cgc0(std::move(other.cgc0)), cgc1(std::move(other.cgc1)),
      cgcComp(std::move(other.cgcComp)), config(other.config),
      pairs0(other.pairs0), pairs1(other.pairs1), pairsComp(other.pairsComp),
      matches(other.matches), isMatched0(other.isMatched0),
      isMatched1(other.isMatched1), _matched(other.matched),
      _unmatched0(other.unmatched0), _unmatched1(other.unmatched1) {}

void BettiMatching::computeMatching() {
    vector<Cube> ctr0;
    vector<Cube> ctr1;
    vector<Cube> ctrComp;

#ifdef RUNTIME
    cout << "dimension 0:";
    auto start = high_resolution_clock::now();
#endif

    Dimension0 dim0(cgc0, cgc1, cgcComp, config, pairs0[0], pairs1[0],
                    pairsComp[0], matches[0], isMatched0[0], isMatched1[0]);
    dim0.computePairsAndMatch(ctr0, ctr1, ctrComp);

#ifdef RUNTIME
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);
    cout << endl << "total: " << duration.count() << " ms" << endl << endl;
#endif
}

void BettiMatching::computeVoxels() {
#ifdef RUNTIME
    cout << "computing voxels ... ";
    auto start = high_resolution_clock::now();
#endif

    for (auto &pair : pairs0[0]) {
        if (!isMatched0[0][pair.birth.index]) {
            _unmatched0.push_back(
                VoxelPair(std::vector{cgc0.getParentVoxel(pair.birth, 0)},
                          std::vector{cgc0.getParentVoxel(pair.death, 1)}));
        }
    }
    for (auto &pair : pairs1[0]) {
        if (!isMatched1[0][pair.birth.index]) {
            _unmatched1.push_back(
                VoxelPair(std::vector{cgc1.getParentVoxel(pair.birth, 0)},
                          std::vector{cgc1.getParentVoxel(pair.death, 1)}));
        }
    }
    for (auto &match : matches[0]) {
        _matched.push_back(VoxelMatch(
            VoxelPair(std::vector{cgc0.getParentVoxel(match.pair0.birth, 0)},
                      std::vector{cgc0.getParentVoxel(match.pair0.death, 1)}),
            VoxelPair(std::vector{cgc1.getParentVoxel(match.pair1.birth, 0)},
                      std::vector{cgc1.getParentVoxel(match.pair1.death, 1)})));
    }

#ifdef RUNTIME
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);
    cout << duration.count() << " ms" << endl << endl;
#endif
}

vector<vector<VoxelPair>> BettiMatching::computePairsInput0() {
    vector<Cube> ctr0;

    Dimension0 dim0(cgc0, cgc1, cgcComp, config, pairs0[0], pairs1[0],
                    pairsComp[0], matches[0], isMatched0[0], isMatched1[0]);
    dim0.computeInput0Pairs(ctr0);

    vector<vector<VoxelPair>> voxelPairs(1);

    for (uint8_t d = 0; d < 1; ++d) {
        for (auto &pair : pairs0[d]) {
            voxelPairs[d].push_back(
                VoxelPair({cgc0.getParentVoxel(pair.birth, d)},
                          {cgc0.getParentVoxel(pair.death, d + 1)}));
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
    if (cgc0.shape[0] < 10) {
        cgc0.printImage();
    }
    cout << "dim " << 0 << ": ";
    count = pairs0[0].size();
    if (0 < count && count < 10) {
        cout << endl;
        for (auto &pair : pairs0[0]) {
            pair.print();
            cout << endl;
        }
    } else {
        cout << count << endl;
    }

    cout << "------------------------------------------------------------------"
            "---------------------------------------------"
         << endl;
    cout << "Input 1" << endl << endl;
    if (cgc1.shape[0] < 10) {
        cgc1.printImage();
    }
    cout << "dim " << 0 << ": ";
    count = pairs1[0].size();
    if (0 < count && count < 10) {
        cout << endl;
        for (auto &pair : pairs1[0]) {
            pair.print();
            cout << endl;
        }
    } else {
        cout << count << endl;
    }

#ifdef COMPUTE_COMPARISON
    cout << "------------------------------------------------------------------"
            "---------------------------------------------"
         << endl;
    cout << "Comparison" << endl << endl;
    if (cgcComp.shape[0] < 10) {
        cgcComp.printImage();
    }
    cout << "dim " << 0 << ": ";
    count = pairsComp[0].size();
    if (0 < count && count < 10) {
        cout << endl;
        for (auto &pair : pairsComp[0]) {
            pair.print();
            cout << endl;
        }
    } else {
        cout << "number of pairs: " << count << endl;
    }
#endif

    cout << "------------------------------------------------------------------"
            "---------------------------------------------"
         << endl;
    cout << "Betti Matching:" << endl << endl;
    cout << "matched: " << endl;
    cout << "dim " << 0 << ": ";
    count = matches[0].size();
    if (0 < count && count < 10) {
        cout << endl;
        for (size_t i = 0; i < matches[0].size(); i++) {
            matches[0][i].print();
            _matched[i].print();
            if (cgc0.shape[0] < 10 && cgc0.shape[1] < 10 &&
                cgc0.shape[2] < 10) {
                pair<dim1::RepresentativeCycle, dim1::RepresentativeCycle>
                    reprCycles = getMatchedRepresentativeCycles(0, i);
                cgc0.printRepresentativeCycle(get<0>(reprCycles));
                cout << endl;
                cgc1.printRepresentativeCycle(get<1>(reprCycles));
                cout << endl;
            }
        }
    } else {
        cout << count << endl;
    }

    size_t counter;
    cout << endl << "unmatched in Input 0" << endl;
    cout << "dim " << 0 << ": ";
    count = 0;
    for (auto &pair : pairs0[0]) {
        if (!isMatched0[0][pair.birth.index]) {
            ++count;
        }
    }
    if (0 < count && count < 10) {
        cout << endl;
        counter = 0;
        for (size_t i = 0; i < pairs0[0].size(); ++i) {
            if (!isMatched0[0][pairs0[0][i].birth.index]) {
                pairs0[0][i].print();
                cout << endl;
                _unmatched0[counter].print();
                cout << endl;
                dim1::RepresentativeCycle reprCycle =
                    getUnmatchedRepresentativeCycle(0, 0, counter);
                cgc0.printRepresentativeCycle(reprCycle);
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

    cout << endl << "unmatched in Input 1" << endl;
    cout << "dim " << 0 << ": ";
    count = 0;
    for (auto &pair : pairs1[0]) {
        if (!isMatched1[0][pair.birth.index]) {
            ++count;
        }
    }
    if (0 < count && count < 10) {
        cout << endl;
        counter = 0;
        for (size_t i = 0; i < pairs1[0].size(); ++i) {
            if (!isMatched1[0][pairs1[0][i].birth.index]) {
                pairs1[0][i].print();
                cout << endl;
                _unmatched1[counter].print();
                cout << endl;
                dim1::RepresentativeCycle reprCycle =
                    getUnmatchedRepresentativeCycle(1, 0, counter);
                cgc1.printRepresentativeCycle(reprCycle);
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

pair<dim1::RepresentativeCycle, dim1::RepresentativeCycle>
BettiMatching::getMatchedRepresentativeCycles(const uint8_t &dim,
                                              const size_t &index) {
    if (dim >= 1) {
        throw runtime_error("Invalid value for dim");
    }

    pair<dim1::RepresentativeCycle, dim1::RepresentativeCycle> reprCycles;

    if (index >= matches[dim].size()) {
        throw runtime_error("Cycle index " + std::to_string((int)index) +
                            " out of range: There are only " +
                            std::to_string((int)matches[dim].size()) +
                            " matched cycles in dimension " +
                            std::to_string((int)dim));
    }

    Dimension0 dim0(cgc0, cgc1, cgcComp, config, pairs0[0], pairs1[0],
                    pairsComp[0], matches[0], isMatched0[0], isMatched1[0]);
    get<0>(reprCycles) =
        dim0.getRepresentativeCycle(matches[0][index].pair0, cgc0);
    get<1>(reprCycles) =
        dim0.getRepresentativeCycle(matches[0][index].pair1, cgc1);

    return reprCycles;
}

dim1::RepresentativeCycle BettiMatching::getUnmatchedRepresentativeCycle(
    const uint8_t &input, const size_t &dim, const size_t &index) {
    if (dim >= 1) {
        throw runtime_error("Invalid value for dim");
    }
    if (input >= 2) {
        throw runtime_error("Invalid value for input: must be 0 or 1");
    }

    const CubicalGridComplex &cgc = (input == 0) ? cgc0 : cgc1;
    vector<vector<Pair>> &pairs = (input == 0) ? pairs0 : pairs1;
    vector<unordered_map<uint64_t, bool>> isMatched =
        (input == 0) ? isMatched0 : isMatched1;

    dim1::RepresentativeCycle reprCycle;

    size_t numPairs = pairs[dim].size();
    if (index >= numPairs) {
        throw runtime_error(
            "Cycle index " + std::to_string((int)index) +
            " out of range: There are only " + std::to_string((int)numPairs) +
            " unmatched cycles in dimension " + std::to_string((int)dim));
    }

    size_t count = 0;
    for (Pair &pair : pairs[dim]) {
        if (!isMatched[dim][pair.birth.index]) {
            if (count == index) {
                Dimension0 dim0(cgc0, cgc1, cgcComp, config, pairs0[0],
                                pairs1[0], pairsComp[0], matches[0],
                                isMatched0[0], isMatched1[0]);
                reprCycle = dim0.getRepresentativeCycle(pair, cgc);
                break;
            }
            ++count;
        }
    }

    return reprCycle;
}
