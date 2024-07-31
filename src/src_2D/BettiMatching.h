#include "../data_structures.h"
#include "data_structures.h"
#include <unordered_map>

namespace dim2 {
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
    const vector<vector<VoxelMatch>> &matched = _matched;
    const vector<vector<VoxelPair>> &unmatched0 = _unmatched0;
    const vector<vector<VoxelPair>> &unmatched1 = _unmatched1;
    tuple<vector<dim2::RepresentativeCycle>, vector<dim2::RepresentativeCycle>>
    computeRepresentativeCycles(
        const int input, const int dim,
        const optional<vector<size_t>> &matchedPairsIndices,
        const optional<vector<size_t>> &unmatchedPairsIndices);

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
    vector<vector<VoxelMatch>> _matched;
    vector<vector<VoxelPair>> _unmatched0;
    vector<vector<VoxelPair>> _unmatched1;
    Config config;
};
} // namespace dim2
