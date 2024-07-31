#include "../data_structures.h"
#include "data_structures.h"
#include <unordered_map>

namespace dim1 {
class BettiMatching {
  public:
    BettiMatching(vector<value_t> &&input0, vector<value_t> &&input1,
                  vector<value_t> &&comparison, vector<index_t> &&shape,
                  Config &&config);
    BettiMatching(BettiMatching &&other);
    void computeMatching();
    void computeVoxels();
    vector<vector<VoxelPair>> computePairsInput0();
    void printResult();
    tuple<vector<dim1::RepresentativeCycle>, vector<dim1::RepresentativeCycle>>
    computeRepresentativeCycles(
        const int input, const optional<vector<size_t>> &matchedPairsIndices,
        const optional<vector<size_t>> &unmatchedPairsIndices);
    const vector<VoxelMatch> &matched = _matched;
    const vector<VoxelPair> &unmatched0 = _unmatched0;
    const vector<VoxelPair> &unmatched1 = _unmatched1;

  private:
    CubicalGridComplex cgc0;
    CubicalGridComplex cgc1;
    CubicalGridComplex cgcComp;
    vector<vector<Pair>> pairs0;
    vector<vector<Pair>> pairs1;
    vector<vector<Pair>> pairsComp;
    vector<vector<Match>> matches;
    vector<unordered_map<uint64_t, bool>> isMatched0;
    vector<unordered_map<uint64_t, bool>> isMatched1;
    vector<unordered_map<uint64_t, size_t>> isMatchedWithIndexComp;
    vector<VoxelMatch> _matched;
    vector<VoxelPair> _unmatched0;
    vector<VoxelPair> _unmatched1;
    Config config;
};
} // namespace dim1
