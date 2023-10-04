#include "BettiMatching.h"
#include "dimension_0.h"
#include "dimension_1.h"
#include "dimension_2.h"

#include <iostream>
#include <chrono>

using namespace dim3;
using namespace std;
using namespace std::chrono;


BettiMatching::BettiMatching(vector<value_t> input0, vector<value_t> input1, vector<value_t> comparison, vector<index_t> shape,
                                Config& _config) : 
    cgc0(input0, shape), cgc1(input1, shape), cgcComp(comparison, shape), config(_config) {
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

BettiMatching::BettiMatching(BettiMatching &&other) : cgc0(std::move(other.cgc0)), cgc1(std::move(other.cgc1)), cgcComp(std::move(other.cgcComp)),
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

void BettiMatching::printResult() {
    index_t count;
    cout << "---------------------------------------------------------------------------------------------------------------" << endl;
    cout << "Input 0:" << endl << endl;
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
    cout << "Input 1:" << endl << endl; 
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
    cout << "Comparison:" << endl << endl; 
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
    cout << "matched cubes: " << endl;
    for (uint8_t d = 0; d < 3; ++d) {
        cout << "dim " << unsigned(d) << ": ";
        count = matches[d].size();
        if (0 < count && count < 10) {
            cout << endl;
            for (auto &match : matches[d]) { match.print(); }
        } else { cout << count << endl; }
    }
    cout << endl;
    cout << "unmatched cubes in Input 0:" << endl;
    for (uint8_t d = 0; d < 3; ++d) {
        cout << "dim " << unsigned(d) << ": ";
        count = 0;
        for (auto &pair : pairs0[d]) { if (!isMatched0[d][pair.birth.index]) { ++count; } }
        if (0 < count && count < 10) {
            cout << endl;
            for (auto &pair : pairs0[d]) { if (!isMatched0[d][pair.birth.index]) { pair.print(); cout << endl; } }
        } else { cout << count << endl; }
    }
    cout << endl;
    cout << "unmatched cubes in Input 1:" << endl;
    for (uint8_t d = 0; d < 3; ++d) {
        cout << "dim " << unsigned(d) << ": ";
        count = 0;
        for (auto &pair : pairs1[d]) {if (!isMatched1[d][pair.birth.index]) { ++count; } }
        if (0 < count && count < 10) {
            cout << endl;
            for (auto &pair : pairs1[d]) { if (!isMatched1[d][pair.birth.index]) { pair.print(); cout << endl; } }
        } else { cout << count << endl; }
    }
    cout << endl;
    cout << "matched voxels:" << endl;
    for (uint8_t d = 0; d < 3; ++d) {
        cout << "dim " << unsigned(d) << ": ";
        count = _matched[d].size();
        if (0 < count && count < 10) {
            cout << endl;
            for (auto& match : _matched[d]) { match.print(); }
        } else { cout << count << endl; }
    }
    cout << endl;
    cout << "unmatched voxels in Input 0:" << endl;
    for (uint8_t d = 0; d < 3; ++d) {
        cout << "dim " << unsigned(d) << ": ";
        count = _unmatched0[d].size();
        if (0 < count && count < 10) {
            cout << endl;
            for (auto& pair : _unmatched0[d]) { pair.print(); cout << endl; }
        } else { cout << count << endl; }
    }
    cout << endl;
    cout << "unmatched voxels in Input 1:" << endl;
    for (uint8_t d = 0; d < 3; ++d) {
        cout << "dim " << unsigned(d) << ": ";
        count = _unmatched1[d].size();
        if (0 < count && count < 10) {
            cout << endl;
            for (auto& pair : _unmatched1[d]) { pair.print(); cout << endl; }
        } else { cout << count << endl; }
    }
    cout << endl;
}