#pragma once

#include "data_structures.h"

#include <cstdint>
#include <queue>
#include <unordered_map>

namespace dimN {
typedef priority_queue<Cube, vector<Cube>, CubeComparator> CubeQueue;

class InterDimensions {
  public:
    InterDimensions(const CubicalGridComplex &cgc0,
                    const CubicalGridComplex &cgc1,
                    const CubicalGridComplex &cgcComp, const Config &config,
                    vector<vector<Pair>> &pairs0, vector<vector<Pair>> &pairs1,
                    vector<vector<Pair>> &pairsComp,
                    vector<vector<Match>> &matches,
                    unordered_map<index_t, bool> &isMatched0,
                    unordered_map<index_t, bool> &isMatched1);
    void computePairsAndMatch(vector<Cube> &ctr0, vector<Cube> &ctr1,
                              vector<Cube> &ctrComp);
    void computeInput0Pairs(vector<Cube> &ctr0);

  private:
    const CubicalGridComplex &cgc0;
    const CubicalGridComplex &cgc1;
    const CubicalGridComplex &cgcComp;
    const Config &config;
    vector<vector<Pair>> &pairs0;
    vector<vector<Pair>> &pairs1;
    vector<vector<Pair>> &pairsComp;
    vector<vector<Match>> &matches;
    unordered_map<index_t, bool> &isMatched0;
    unordered_map<index_t, bool> &isMatched1;
    unordered_map<index_t, bool> isMatchedComp;
    unordered_map<index_t, Pair> matchMap0;
    unordered_map<index_t, Pair> matchMap1;
    unordered_map<index_t, Cube> matchMapIm0;
    unordered_map<index_t, Cube> matchMapIm1;
    unordered_map<index_t, index_t> pivotColumnIndex;
    unordered_map<index_t, CubeQueue> cache;
    index_t computeDim;

    void computePairs(const vector<Cube> &ctr, uint8_t k);
    void computePairsComp(vector<Cube> &ctr);
    void computeImagePairs(vector<Cube> &ctr, uint8_t k);
    void computeMatching();
    void assembleColumnsToReduce(const CubicalGridComplex &cgc,
                                 vector<Cube> &ctr) const;
    Cube popPivot(CubeQueue &column) const;
    Cube getPivot(CubeQueue &column) const;
    void addCache(const index_t &i, CubeQueue &working_boundary);
};
} // namespace dimN
