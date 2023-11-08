#pragma once

#include "config.h"
#include "src_1D/BettiMatching.h"
#include "src_2D/BettiMatching.h"
#include "src_3D/BettiMatching.h"
#include "src_nD/BettiMatching.h"

#include <variant>
#include <optional>

using namespace std;

class BettiMatching {
public:
    BettiMatching(vector<value_t>&& input0, vector<value_t>&& input1, vector<index_t>&& shape, Config&& config);
    BettiMatching(BettiMatching&& other);
    void computeMatching();
    vector<vector<VoxelPair>> computePairsInput0();
    void printResult();
    tuple<vector<vector<VoxelMatch>>, vector<vector<VoxelPair>>, vector<vector<VoxelPair>>> getMatching();
    tuple<vector<vector<index_t>>, vector<vector<index_t>>> getMatchedRepresentativeCycles(const size_t& dim, const size_t& index);
    vector<vector<index_t>> getUnmatchedRepresentativeCycle(const size_t& dim, const size_t& index, const uint8_t& input);
    vector<index_t> shape;


private:
    optional<variant<dim1::BettiMatching, dim2::BettiMatching, dim3::BettiMatching, dimN::BettiMatching>> dimensionSpecificBettiMatching;
    size_t dimension;
    bool computed;
    bool computedPairsInput0;
};