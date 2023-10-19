#include "BettiMatching.h"
#include "dimension_0.h"
#include "inter_dimensions.h"
#include "top_dimension.h"

#include <iostream>
#include <chrono>

using namespace dimN;
using namespace std;
using namespace std::chrono;



BettiMatching::BettiMatching(vector<value_t>&& input0, vector<value_t>&& input1, vector<value_t>&& comparison, vector<index_t>&& shape,
                                Config&& _config) : cgc0(input0, shape), cgc1(input1, shape), cgcComp(comparison, shape), config(_config) {
    dim = shape.size();
    pairs0 = vector<vector<Pair>>(dim);
    pairs1 = vector<vector<Pair>>(dim);
    pairsComp = vector<vector<Pair>>(dim);
    matches = vector<vector<Match>>(dim);
    _matched = vector<vector<VoxelMatch>>(dim);
    _unmatched0 = vector<vector<VoxelPair>>(dim);
    _unmatched1 = vector<vector<VoxelPair>>(dim);
}


void BettiMatching::computeMatching() {
    vector<Cube> ctr0;
    vector<Cube> ctr1;
    vector<Cube> ctrComp;

    {
#ifdef RUNTIME
        cout << "dimension " << dim-1 << ":";
        auto start = high_resolution_clock::now();
#endif
        TopDimension topDim(cgc0, cgc1, cgcComp, config, pairs0[dim-1], pairs1[dim-1], pairsComp[dim-1], matches[dim-1],
                            isMatched0, isMatched1);
        topDim.computePairsAndMatch(ctr0, ctr1, ctrComp);
#ifdef RUNTIME
        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<milliseconds>(stop - start);
        cout << endl << "total: " << duration.count() << " ms" << endl << endl;
#endif
    }

    {
#ifdef RUNTIME
        auto start = high_resolution_clock::now();
#endif
        InterDimensions interDim(cgc0, cgc1, cgcComp, config, pairs0, pairs1, pairsComp, matches, isMatched0, isMatched1);
        interDim.computePairsAndMatch(ctr0, ctr1, ctrComp);
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
        Dimension0 dim0(cgc0, cgc1, cgcComp,  config, pairs0[0], pairs1[0], pairsComp[0], matches[0], isMatched0, isMatched1);       
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

    for (uint8_t d = 0; d < dim; ++d) {
        for (auto& pair : pairs0[d]) {
            if (!isMatched0[cgc0.getCubeIndex(pair.birth)]) {
                _unmatched0[d].push_back(VoxelPair(cgc0.getParentVoxel(pair.birth), cgc0.getParentVoxel(pair.death)));
            }
        }
        for (auto& pair : pairs1[d]) {
            if (!isMatched1[cgc1.getCubeIndex(pair.birth)]) {
                _unmatched1[d].push_back(VoxelPair(cgc1.getParentVoxel(pair.birth), cgc1.getParentVoxel(pair.death)));
            }
        }
        for (auto& match : matches[d]) {
            _matched[d].push_back(VoxelMatch(VoxelPair(cgc0.getParentVoxel(match.pair0.birth),
                                                        cgc0.getParentVoxel(match.pair0.death)), 
                                            VoxelPair(cgc1.getParentVoxel(match.pair1.birth),
                                                        cgc1.getParentVoxel(match.pair1.death))));
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
    for (uint8_t d = 0; d < dim; ++d) {
        cout << "dim " << unsigned(d) << ": ";
        count = pairs0[d].size();
        if (0 < count && count < 10) {
            cout << endl;
            for (auto &pair : pairs0[d]) { pair.print(); cout << endl; }
        } else { cout << count << endl; }
    }
    cout << "---------------------------------------------------------------------------------------------------------------" << endl;
    cout << "Input 1" << endl << endl; 
    for (uint8_t d = 0; d < dim; ++d) {
        cout << "dim " << unsigned(d) << ": ";
        count = pairs1[d].size();
        if (0 < count && count < 10) {
            cout << endl;
            for (auto &pair : pairs1[d]) { pair.print(); cout << endl; }
        } else { cout << count << endl; }
    }

#ifdef COMPUTE_COMPARISON
    cout << "---------------------------------------------------------------------------------------------------------------" << endl;
    cout << "Comparison" << endl << endl; 
    for (uint8_t d = 0; d < dim; ++d) {
        cout << "dim " << unsigned(d) << ": ";
        count = pairsComp[d].size();
        if (0 < count && count < 10) {
            cout << endl;
            for (auto &pair : pairsComp[d]) { pair.print(); cout << endl; }
        } else { cout << "number of pairs: " << count << endl; }
    }
#endif

    cout << "---------------------------------------------------------------------------------------------------------------" << endl;
    cout << "Betti Matching:" << endl << endl;
    cout << "matched cubes: " << endl;
    for (uint8_t d = 0; d < dim; ++d) {
        cout << "dim " << unsigned(d) << ": ";
        count = matches[d].size();
        if (0 < count && count < 10) {
            cout << endl;
            for (auto &match : matches[d]) { match.print(); }
        } else { cout << count << endl; }
    }
    cout << endl;
    cout << "unmatched cubes in Input 0" << endl;
    for (uint8_t d = 0; d < dim; ++d) {
        cout << "dim " << unsigned(d) << ": ";
        count = 0;
        for (auto &pair : pairs0[d]) { if (!isMatched0[cgc0.getCubeIndex(pair.birth)]) { ++count; } }
        if (0 < count && count < 10) {
            cout << endl;
            for (auto &pair : pairs0[d]) { if (!isMatched0[cgc0.getCubeIndex(pair.birth)]) { pair.print(); cout << endl; } }
        } else { cout << count << endl; }
    }
    cout << endl;
    cout << "unmatched cubes in Input 1" << endl;
    for (uint8_t d = 0; d < dim; ++d) {
        cout << "dim " << unsigned(d) << ": ";
        count = 0;
        for (auto &pair : pairs1[d]) {if (!isMatched1[cgc1.getCubeIndex(pair.birth)]) { ++count; } }
        if (0 < count && count < 10) {
            cout << endl;
            for (auto &pair : pairs1[d]) { if (!isMatched1[cgc1.getCubeIndex(pair.birth)]) { pair.print(); cout << endl; } }
        } else { cout << count << endl; }
    }
    cout << endl;
    cout << "matched voxels: " << endl;
    for (uint8_t d = 0; d < dim; ++d) {
        cout << "dim " << unsigned(d) << ": ";
        count = _matched[d].size();
        if (0 < count && count < 10) {
            cout << endl;
            for (auto& match : _matched[d]) { match.print(); }
        } else { cout << count << endl; }
    }
    cout << endl;
    cout << "unmatched voxels in Input 0" << endl;
    for (uint8_t d = 0; d < dim; ++d) {
        cout << "dim " << unsigned(d) << ": ";
        count = _unmatched0[d].size();
        if (0 < count && count < 10) {
            cout << endl;
            for (auto& pair : _unmatched0[d]) { pair.print(); cout << endl; }
        } else { cout << count << endl; }
    }
    cout << endl;
    cout << "unmatched voxels in Input 1" << endl;
    for (uint8_t d = 0; d < dim; ++d) {
        cout << "dim " << unsigned(d) << ": ";
        count = _unmatched1[d].size();
        if (0 < count && count < 10) {
            cout << endl;
            for (auto& pair : _unmatched1[d]) { pair.print(); cout << endl; }
        } else { cout << count << endl; }
    }
    cout << endl;
}