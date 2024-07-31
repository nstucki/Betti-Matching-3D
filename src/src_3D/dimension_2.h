#pragma once

#include "../config.h"
#include "BettiMatching.h"
#include "data_structures.h"
#include "enumerators.h"

#include <unordered_map>

namespace dim3 {
class Dimension2 {
  public:
    Dimension2(const CubicalGridComplex &cgc0, const CubicalGridComplex &cgc1,
               const CubicalGridComplex &cgcComp, const Config &config,
               vector<Pair> &pairs0, vector<Pair> &pairs1,
               vector<Pair> &pairsComp, vector<Match> &matches,
               unordered_map<uint64_t, bool> &isMatched0,
               unordered_map<uint64_t, bool> &isMatched1,
               unordered_map<uint64_t, size_t> &isMatchedWithIndexComp);
    void computePairsAndMatch(vector<Cube> &ctr0, vector<Cube> &ctr1,
                              vector<Cube> &ctrComp, vector<Cube> &ctrImage);
    void computeInput0Pairs(vector<Cube> &ctr0);
    RepresentativeCycle
    computeRepresentativeCycle(const Pair &pair,
                               const CubicalGridComplex &cgc) const;
    vector<dim3::RepresentativeCycle> computeRepresentativeCycles(
        const int input,
        const std::vector<std::reference_wrapper<Pair>> &requestedPairs);

  private:
    void enumerateDualEdges(vector<Cube> &dualEdges,
                            const CubicalGridComplex &cgc) const;
    void enumerateDualEdgesComp(vector<Cube> &dualEdges) const;
    void computeInputAndImagePairs(vector<Cube> &dualEdges, const uint8_t &k);
    void computeCompPairsAndMatch(vector<Cube> &dualEdges,
                                  vector<Cube> &ctrImage);
#if defined(USE_APPARENT_PAIRS) or defined(USE_APPARENT_PAIRS_COMP)
    bool isApparentPair(const Cube &dualEdge, BoundaryEnumerator &enumerator,
                        CoboundaryEnumerator &coEnumerator) const;
#endif
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
};
} // namespace dim3
