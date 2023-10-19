#include "BettiMatching.h"
#include "../data_structures.h"
#include "dimension_0.h"

#include <iostream>
#include <chrono>

using namespace dim1;
using namespace std;
using namespace std::chrono;



BettiMatching::BettiMatching(vector<value_t>&& input0, vector<value_t>&& input1, vector<value_t>&& comparison, vector<index_t>&& shape,
                                Config&& _config) : cgc0(input0, shape), cgc1(input1, shape), cgcComp(comparison, shape), config(_config) {}


BettiMatching::BettiMatching(BettiMatching&& other) : 
    cgc0(std::move(other.cgc0)), cgc1(std::move(other.cgc1)), cgcComp(std::move(other.cgcComp)),
    config(other.config), pairs0(other.pairs0), pairs1(other.pairs1), pairsComp(other.pairsComp),
    matches(other.matches), isMatched0(other.isMatched0), isMatched1(other.isMatched1),
    _matched(other.matched), _unmatched0(other.unmatched0), _unmatched1(other.unmatched1) {}


void BettiMatching::computeMatching() {
    vector<Cube> ctr0;
    vector<Cube> ctr1;
    vector<Cube> ctrComp;

#ifdef RUNTIME
        cout << "dimension 0:";
        auto start = high_resolution_clock::now();
#endif
        Dimension0 dim0(cgc0, cgc1, cgcComp,  config, pairs0, pairs1, pairsComp, matches, isMatched0, isMatched1);       
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

    for (auto& pair : pairs0) {
        if (!isMatched0[pair.birth.index]) {
            _unmatched0.push_back(VoxelPair(std::vector{cgc0.getParentVoxel(pair.birth, 0)}, std::vector{cgc0.getParentVoxel(pair.death, 1)}));
        }
    }
    for (auto& pair : pairs1) {
        if (!isMatched1[pair.birth.index]) {
            _unmatched1.push_back(VoxelPair(std::vector{cgc1.getParentVoxel(pair.birth, 0)}, std::vector{cgc1.getParentVoxel(pair.death, 1)}));
        }
    }
    for (auto& match : matches) {
        _matched.push_back(VoxelMatch(VoxelPair(std::vector{cgc0.getParentVoxel(match.pair0.birth, 0)}, 
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


void BettiMatching::printResult() {
    index_t count;
    cout << "---------------------------------------------------------------------------------------------------------------" << endl;
    cout << "Input 0:" << endl << endl;
    if (cgc0.shape[0] < 10) { cgc0.printImage(); }
    cout << "dim " << 0 << ": ";
    count = pairs0.size();
    if (0 < count && count < 10) {
        cout << endl;
        for (auto &pair : pairs0) { pair.print(); cout << endl; }
    } else { cout << count << endl; }
    cout << "---------------------------------------------------------------------------------------------------------------" << endl;
    cout << "Input 1" << endl << endl; 
    if (cgc1.shape[0] < 10) { cgc1.printImage(); }
    cout << "dim " << 0 << ": ";
    count = pairs1.size();
    if (0 < count && count < 10) {
        cout << endl;
        for (auto &pair : pairs1) { pair.print(); cout << endl; }
    } else { cout << count << endl; }

#ifdef COMPUTE_COMPARISON
    cout << "---------------------------------------------------------------------------------------------------------------" << endl;
    cout << "Comparison" << endl << endl; 
    if (cgcComp.shape[0] < 10) { cgcComp.printImage(); }
    cout << "dim " << 0 << ": ";
    count = pairsComp.size();
    if (0 < count && count < 10) {
        cout << endl;
        for (auto &pair : pairsComp) { pair.print(); cout << endl; }
    } else { cout << "number of pairs: " << count << endl; }
#endif

    cout << "---------------------------------------------------------------------------------------------------------------" << endl;
    cout << "Betti Matching:" << endl << endl;
    cout << "matched cubes: " << endl;
    cout << "dim " << 0 << ": ";
    count = matches.size();
    if (0 < count && count < 10) {
        cout << endl;
        for (auto &match : matches) { match.print(); }
    } else { cout << count << endl; }
    cout << endl;
    cout << "unmatched cubes in Input 0" << endl;
    cout << "dim " << 0 << ": ";
    count = 0;
    for (auto &pair : pairs0) { if (!isMatched0[pair.birth.index]) { ++count; } }
    if (0 < count && count < 10) {
        cout << endl;
        for (auto &pair : pairs0) { if (!isMatched0[pair.birth.index]) { pair.print(); cout << endl; } }
    } else { cout << count << endl; }
    cout << endl;
    cout << "unmatched cubes in Input 1" << endl;
    cout << "dim " << 0 << ": ";
    count = 0;
    for (auto &pair : pairs1) {if (!isMatched1[pair.birth.index]) { ++count; } }
    if (0 < count && count < 10) {
        cout << endl;
        for (auto &pair : pairs1) { if (!isMatched1[pair.birth.index]) { pair.print(); cout << endl; } }
    } else { cout << count << endl; }
    cout << endl;
    cout << "matched voxels: " << endl;
    cout << "dim " << 0 << ": ";
    count = _matched.size();
    if (0 < count && count < 10) {
        cout << endl;
        for (auto& match : _matched) { match.print(); }
    } else { cout << count << endl; }
    cout << endl;
    cout << "unmatched voxels in Input 0" << endl;
    cout << "dim " << 0 << ": ";
    count = _unmatched0.size();
    if (0 < count && count < 10) {
        cout << endl;
        for (auto& pair : _unmatched0) { pair.print(); cout << endl; }
    } else { cout << count << endl; }
    cout << endl;
    cout << "unmatched voxels in Input 1" << endl;
    cout << "dim " << 0 << ": ";
    count = _unmatched1.size();
    if (0 < count && count < 10) {
        cout << endl;
        for (auto& pair : _unmatched1) { pair.print(); cout << endl; }
    } else { cout << count << endl; }
    cout << endl;
}