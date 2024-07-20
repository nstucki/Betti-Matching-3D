#include "BettiMatching.h"
#include "dimension_0.h"
#include "dimension_1.h"
#include "dimension_2.h"

#include <iostream>
#include <chrono>

using namespace dim3;
using namespace std;
using namespace std::chrono;


BettiMatching::BettiMatching(vector<value_t>&& input0, vector<value_t>&& input1, vector<value_t>&& comparison, vector<index_t>&& shape,
                                Config&& _config) : cgc0(input0, shape), cgc1(input1, shape), cgcComp(comparison, shape), config(_config) {
    pairs0 = vector<vector<Pair>>(3);
    pairs1 = vector<vector<Pair>>(3);
    pairsComp = vector<vector<Pair>>(3);
    matches = vector<vector<Match>>(3);
    isMatched0 = vector<unordered_map<uint64_t, bool>>(3);
    isMatched1 = vector<unordered_map<uint64_t, bool>>(3);
    _matched = vector<vector<VoxelMatch>>(3);
    _unmatched0 = vector<vector<VoxelPair>>(3);
    _unmatched1 = vector<vector<VoxelPair>>(3);
}


 BettiMatching::BettiMatching(BettiMatching&& other) : 
    cgc0(std::move(other.cgc0)), cgc1(std::move(other.cgc1)), cgcComp(std::move(other.cgcComp)),
    config(other.config), pairs0(other.pairs0), pairs1(other.pairs1), pairsComp(other.pairsComp),
    matches(other.matches), isMatched0(other.isMatched0), isMatched1(other.isMatched1),
    _matched(other.matched), _unmatched0(other.unmatched0), _unmatched1(other.unmatched1) {}



void BettiMatching::computeMatching() {
    vector<Cube> ctr0;
    vector<Cube> ctr1;
    vector<Cube> ctrComp;
    vector<Cube> ctrImage;

    {
#ifdef RUNTIME
        cout << "dimension 2:";
        auto start = high_resolution_clock::now();
#endif

        Dimension2 dim2(cgc0, cgc1, cgcComp,  config, pairs0[2], pairs1[2], pairsComp[2], matches[2], isMatched0[2], isMatched1[2]);       
        dim2.computePairsAndMatch(ctr0, ctr1, ctrComp, ctrImage);

#ifdef RUNTIME
        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<milliseconds>(stop - start);
        cout << endl << "total: " << duration.count() << " ms" << endl << endl;
#endif
    }

    {
#ifdef RUNTIME
        cout << "dimension 1:";
        auto start = high_resolution_clock::now();
#endif

        Dimension1 dim1(cgc0, cgc1, cgcComp,  config, pairs0[1], pairs1[1], pairsComp[1], matches[1], isMatched0[1], isMatched1[1]);       
        dim1.computePairsAndMatch(ctr0, ctr1, ctrComp, ctrImage);

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

        Dimension0 dim0(cgc0, cgc1, cgcComp,  config, pairs0[0], pairs1[0], pairsComp[0], matches[0], isMatched0[0], isMatched1[0]);       
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

    for (uint8_t d = 0; d < 3; ++d) {
        for (auto& pair : pairs0[d]) {
            if (!isMatched0[d][pair.birth.index]) {
                _unmatched0[d].push_back(VoxelPair(cgc0.getParentVoxel(pair.birth, d), cgc0.getParentVoxel(pair.death, d+1)));
            }
        }
        for (auto& pair : pairs1[d]) {
            if (!isMatched1[d][pair.birth.index]) {
                _unmatched1[d].push_back(VoxelPair(cgc1.getParentVoxel(pair.birth, d), cgc1.getParentVoxel(pair.death, d+1)));
            }
        }
        for (auto& match : matches[d]) {
            _matched[d].push_back(VoxelMatch(VoxelPair(cgc0.getParentVoxel(match.pair0.birth, d),
                                                        cgc0.getParentVoxel(match.pair0.death, d+1)), 
                                                VoxelPair(cgc1.getParentVoxel(match.pair1.birth, d),
                                                            cgc1.getParentVoxel(match.pair1.death, d+1))));
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

    Dimension2 dim2(cgc0, cgc1, cgcComp,  config, pairs0[2], pairs1[2], pairsComp[2], matches[2], isMatched0[2], isMatched1[2]);
    dim2.computeInput0Pairs(ctr0);

    Dimension1 dim1(cgc0, cgc1, cgcComp,  config, pairs0[1], pairs1[1], pairsComp[1], matches[1], isMatched0[1], isMatched1[1]);
    dim1.computeInput0Pairs(ctr0);

    Dimension0 dim0(cgc0, cgc1, cgcComp,  config, pairs0[0], pairs1[0], pairsComp[0], matches[0], isMatched0[0], isMatched1[0]);
    dim0.computeInput0Pairs(ctr0);

    vector<vector<VoxelPair>> voxelPairs(3);

    for (uint8_t d = 0; d < 3; ++d) {
        for (auto& pair : pairs0[d]) {
            voxelPairs[d].push_back(VoxelPair(cgc0.getParentVoxel(pair.birth, d), cgc0.getParentVoxel(pair.death, d+1)));
        }
    }
    return voxelPairs;
}


void BettiMatching::printResult() {
    index_t count;

    cout << "---------------------------------------------------------------------------------------------------------------" << endl;
    cout << "Input 0:" << endl;
    if (cgc0.shape[0] < 10 && cgc0.shape[1] < 10 && cgc0.shape[2] < 10) { cgc0.printImage(); }
    for (uint8_t d = 0; d < 3; ++d) {
        cout << "dim " << unsigned(d) << ": ";
        count = pairs0[d].size();
        if (0 < count && count < 10) {
            cout << endl;
            for (auto &pair : pairs0[d]) { pair.print(); cout << endl; }
        } else { cout << count << endl; }
    }

    cout << "---------------------------------------------------------------------------------------------------------------" << endl;
    cout << "Input 1:" << endl; 
    if (cgc1.shape[0] < 10 && cgc1.shape[1] < 10 && cgc1.shape[2] < 10) { cgc1.printImage(); }
    for (uint8_t d = 0; d < 3; ++d) {
        cout << "dim " << unsigned(d) << ": ";
        count = pairs1[d].size();
        if (0 < count && count < 10) {
            cout << endl;
            for (auto &pair : pairs1[d]) { pair.print(); cout << endl; }
        } else { cout << count << endl; }
    }

#ifdef COMPUTE_COMPARISON
    cout << "---------------------------------------------------------------------------------------------------------------" << endl;
    cout << "Comparison:" << endl; 
    if (cgcComp.shape[0] < 10 && cgcComp.shape[1] < 10 && cgcComp.shape[2] < 10) { cgcComp.printImage(); }
    for (uint8_t d = 0; d < 3; ++d) {
        cout << "dim " << unsigned(d) << ": ";
        count = pairsComp[d].size();
        if (0 < count && count < 10) {
            cout << endl;
            for (auto &pair : pairsComp[d]) { pair.print(); cout << endl; }
        } else { cout << count << endl; }
    }
#endif

    cout << "---------------------------------------------------------------------------------------------------------------" << endl;
    cout << "Betti Matching:" << endl << endl;
    cout << "matched: " << endl;
    for (uint8_t d = 0; d < 3; ++d) {
        cout << "dim " << unsigned(d) << ": ";
        count = matches[d].size();
        if (0 < count && count < 10) {
            cout << endl;
            for (size_t i = 0; i < matches[d].size(); i++) {
                matches[d][i].print();
                _matched[d][i].print();
                if (cgc0.shape[0] < 10 && cgc0.shape[1] < 10 && cgc0.shape[2] < 10) {
                    pair<vector<vector<index_t>>, vector<vector<index_t>>> reprCycles = getMatchedRepresentativeCycles(d, i);
                    cgc0.printRepresentativeCycle(get<0>(reprCycles));
                    cout << endl;
                    cgc1.printRepresentativeCycle(get<1>(reprCycles));
                    cout << endl;
                }
            }
        } else { cout << count << endl; }
    }

    size_t counter;
    cout << endl << "unmatched in Input 0:" << endl;
    for (uint8_t d = 0; d < 3; ++d) {
        cout << "dim " << unsigned(d) << ": ";
        count = 0;
        for (auto &pair : pairs0[d]) { if (!isMatched0[d][pair.birth.index]) { ++count; } }
        if (0 < count && count < 10) {
            cout << endl;
            counter = 0;
            for (size_t i = 0; i < pairs0[d].size(); ++i) { 
                if (!isMatched0[d][pairs0[d][i].birth.index]) {
                    pairs0[d][i].print(); cout << endl;
                    _unmatched0[d][counter].print(); cout << endl;
                    vector<vector<index_t>> reprCycle = getUnmatchedRepresentativeCycle(0, d, counter);
                    cgc0.printRepresentativeCycle(reprCycle);
                    cout << endl;
                    ++counter;
                    if (counter == count) { break; }
                }
            }
        } else { cout << count << endl; }
    }

    cout << endl << "unmatched in Input 1:" << endl;
    for (uint8_t d = 0; d < 3; ++d) {
        cout << "dim " << unsigned(d) << ": ";
        count = 0;
        for (auto &pair : pairs1[d]) { if (!isMatched1[d][pair.birth.index]) { ++count; } }
        if (0 < count && count < 10) {
            cout << endl;
            counter = 0;
            for (size_t i = 0; i < pairs1[d].size(); ++i) { 
                if (!isMatched1[d][pairs1[d][i].birth.index]) {
                    pairs1[d][i].print(); cout << endl;
                    _unmatched1[d][counter].print(); cout << endl;
                    vector<vector<index_t>> reprCycle = getUnmatchedRepresentativeCycle(1, d, counter);
                    cgc1.printRepresentativeCycle(reprCycle);
                    cout << endl;
                    ++counter;
                    if (counter == count) { break; }
                 }
            }
        } else { cout << count << endl; }
    }
}


pair<vector<vector<index_t>>, vector<vector<index_t>>> BettiMatching::getMatchedRepresentativeCycles(const uint8_t& dim, const size_t& index) {
    pair<vector<vector<index_t>>, vector<vector<index_t>>> reprCycles;

    if (index >= matches[dim].size()) { return reprCycles; }

    switch(dim) {
        case 0: {
            Dimension0 dim0(cgc0, cgc1, cgcComp, config, pairs0[0], pairs1[0], pairsComp[0], matches[0], isMatched0[0], isMatched1[0]);
            get<0>(reprCycles) = dim0.getRepresentativeCycle(matches[0][index].pair0, cgc0);
            get<1>(reprCycles) = dim0.getRepresentativeCycle(matches[0][index].pair1, cgc1); 
            break;
        }

        case 1: {
            Dimension1 dim1(cgc0, cgc1, cgcComp, config, pairs0[1], pairs1[1], pairsComp[1], matches[1], isMatched0[1], isMatched1[1]);
            get<0>(reprCycles) = dim1.getRepresentativeCycle(matches[1][index].pair0, cgc0);
            get<1>(reprCycles) = dim1.getRepresentativeCycle(matches[1][index].pair1, cgc1);
            break;
        }

        case 2: {
            Dimension2 dim2(cgc0, cgc1, cgcComp, config, pairs0[2], pairs1[2], pairsComp[2], matches[2], isMatched0[2], isMatched1[2]);
            get<0>(reprCycles) = dim2.getRepresentativeCycle(matches[2][index].pair0, cgc0);
            get<1>(reprCycles) = dim2.getRepresentativeCycle(matches[2][index].pair1, cgc1);
            break;
        }
    }
    return reprCycles;
}


vector<vector<index_t>> BettiMatching::getUnmatchedRepresentativeCycle(const uint8_t& input, const uint8_t& dim, const size_t& index) {
    const CubicalGridComplex& cgc = (input == 0) ? cgc0 : cgc1;
    vector<vector<Pair>>& pairs = (input == 0) ? pairs0 : pairs1;
    vector<unordered_map<uint64_t, bool>> isMatched = (input == 0) ? isMatched0 : isMatched1;

    vector<vector<index_t>> reprCycle;

    size_t numPairs = pairs[dim].size();
    if (index >= numPairs) { return reprCycle; }

    size_t count = 0;
    for (Pair& pair : pairs[dim]) {
        if (!isMatched[dim][pair.birth.index]) {
            if (count == index) {
                switch(dim) {
                    case 0: {
                        Dimension0 dim0(cgc0, cgc1, cgcComp, config, pairs0[0], pairs1[0], pairsComp[0], 
                                        matches[0], isMatched0[0], isMatched1[0]);
                        reprCycle = dim0.getRepresentativeCycle(pair, cgc);
                        break;
                    }

                    case 1: {
                        Dimension1 dim1(cgc0, cgc1, cgcComp, config, pairs0[1], pairs1[1], pairsComp[1], 
                                        matches[1], isMatched0[1], isMatched1[1]);
                        reprCycle = dim1.getRepresentativeCycle(pair, cgc);
                        break;
                    }

                    case 2: {
                        Dimension2 dim2(cgc0, cgc1, cgcComp, config, pairs0[2], pairs1[2], pairsComp[2], 
                                        matches[2], isMatched0[2], isMatched1[2]);
                        reprCycle = dim2.getRepresentativeCycle(pair, cgc);
                        break;
                    }
                }
                break;
            }
            ++count;
        }
    }
    
    return reprCycle;
}