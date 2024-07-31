#pragma once

#include "data_structures.h"

#include <unordered_map>

namespace dim2 {
class Dimension1 {
  public:
    Dimension1(const CubicalGridComplex &cgc0, const CubicalGridComplex &cgc1,
               const CubicalGridComplex &cgcComp, const Config &config,
               vector<Pair> &pairs0, vector<Pair> &pairs1,
               vector<Pair> &pairsComp, vector<Match> &matches,
               unordered_map<uint64_t, bool> &isMatched0,
               unordered_map<uint64_t, bool> &isMatched1,
               unordered_map<uint64_t, size_t> &isMatchedWithIndexComp);
    void computePairsAndMatch(vector<Cube> &ctr0, vector<Cube> &ctr1,
                              vector<Cube> &ctrComp);
    void computeInput0Pairs(vector<Cube> &ctr0);
    vector<dim2::RepresentativeCycle> computeRepresentativeCycles(
        const int input,
        const std::vector<std::reference_wrapper<Pair>> &requestedPairs);

  private:
    const CubicalGridComplex &cgc0;
    const CubicalGridComplex &cgc1;
    const CubicalGridComplex &cgcComp;
    const Config &config;
    vector<Pair> &pairs0;
    vector<Pair> &pairs1;
    vector<Pair> &pairsComp;
    vector<Match> &matches;
    unordered_map<uint64_t, bool> &isMatched0;
    unordered_map<uint64_t, bool> &isMatched1;
    unordered_map<uint64_t, size_t> &isMatchedWithIndexComp;
    unordered_map<index_t, Pair> matchMap0;
    unordered_map<index_t, Pair> matchMap1;
    UnionFindDual uf0;
    UnionFindDual uf1;
    UnionFindDual ufComp;

    void enumerateDualEdges(vector<Cube> &dualEdges,
                            const CubicalGridComplex &cgc) const;
    void computeInputAndImagePairs(vector<Cube> &dualEdges, const uint8_t &k);
    void computeCompPairsAndMatch(vector<Cube> &dualEdges);
};
} // namespace dim2
