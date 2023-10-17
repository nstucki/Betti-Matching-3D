#pragma once

#include <variant>
#include <optional>
#include "src_1D/BettiMatching.h"
#include "src_2D/BettiMatching.h"
#include "src_3D/BettiMatching.h"
#include "src_nD/BettiMatching.h"

using namespace std;

class BettiMatching {
public:
    BettiMatching(vector<value_t> input0, vector<value_t> input1, vector<index_t> shape, Config &config);
    BettiMatching(BettiMatching &&other);
    tuple<vector<vector<VoxelMatch>>, vector<vector<VoxelPair>>, vector<vector<VoxelPair>>> computeMatchingWithVoxels();

private:
    optional<variant<dim1::BettiMatching, dim2::BettiMatching, dim3::BettiMatching, dimN::BettiMatching>> dimensionSpecificBettiMatching;
    int dimension;
};