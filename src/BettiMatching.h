#pragma once

#include "config.h"
#include "src_1D/BettiMatching.h"
#include "src_2D/BettiMatching.h"
#include "src_3D/BettiMatching.h"
#include "src_3D/data_structures.h"
#include "src_nD/BettiMatching.h"

#include <optional>
#include <variant>

using namespace std;

class BettiMatching {
  public:
    BettiMatching(vector<value_t> &&input0, vector<value_t> &&input1,
                  vector<index_t> &&shape, Config &&config);
    BettiMatching(BettiMatching &&other);
    void computeMatching();
    vector<vector<VoxelPair>> computePairsInput0();
    void printResult();
    tuple<vector<vector<VoxelMatch>>, vector<vector<VoxelPair>>,
          vector<vector<VoxelPair>>>
    getMatching();
    variant<tuple<vector<dim1::RepresentativeCycle>,
                  vector<dim1::RepresentativeCycle>>,
            tuple<vector<dim2::RepresentativeCycle>,
                  vector<dim2::RepresentativeCycle>>,
            tuple<vector<dim3::RepresentativeCycle>,
                  vector<dim3::RepresentativeCycle>>>
    computeRepresentativeCycles(
        const int input, const int dim,
        const optional<vector<size_t>> &matchedPairsIndices,
        const optional<vector<size_t>> &unmatchedPairsIndices);
    vector<index_t> shape;

  private:
    optional<variant<dim1::BettiMatching, dim2::BettiMatching,
                     dim3::BettiMatching, dimN::BettiMatching>>
        dimensionSpecificBettiMatching;
    size_t dimension;
    bool computed;
    bool computedPairsInput0;
};
