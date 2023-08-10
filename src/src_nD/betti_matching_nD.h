#pragma once

#include "data_structures.h"
#include "config.h"

#include <unordered_map>

class BettiMatchingND
{
public:
    BettiMatchingND(
        const CubicalGridComplex &cgc0,
        const CubicalGridComplex &cgc1,
        const CubicalGridComplex &cgcComp,
        const Config &config);
    void computeMatching();

    const int dim;
    const CubicalGridComplex &cgc0;
    const CubicalGridComplex &cgc1;
    const CubicalGridComplex &cgcComp;

    // Read-only properties
    const vector<vector<Pair>> &pairs0 = _pairs0;
    const vector<vector<Pair>> &pairs1 = _pairs1;
    const vector<vector<Pair>> &pairsComp = _pairsComp;
    const vector<vector<Match>> &matches = _matches;
    const unordered_map<index_t, bool> &isMatched0 = _isMatched0;
    const unordered_map<index_t, bool> &isMatched1 = _isMatched1;
    const vector<vector<Pair>> &unmatchedPairs0 = _unmatchedPairs0;
    const vector<vector<Pair>> &unmatchedPairs1 = _unmatchedPairs1;

private:
    void computeUnmatched();

    vector<vector<Pair>> _pairs0;
    vector<vector<Pair>> _pairs1;
    vector<vector<Pair>> _pairsComp;
    vector<vector<Match>> _matches;
    unordered_map<index_t, bool> _isMatched0;
    unordered_map<index_t, bool> _isMatched1;
    vector<vector<Pair>> _unmatchedPairs0;
    vector<vector<Pair>> _unmatchedPairs1;

    const Config config;
    vector<Cube> ctr0;
    vector<Cube> ctr1;
    vector<Cube> ctrComp;
};
