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
    variant<pair<dim1::RepresentativeCycle, dim1::RepresentativeCycle>,
            pair<dim2::RepresentativeCycle, dim2::RepresentativeCycle>,
            pair<dim3::RepresentativeCycle, dim3::RepresentativeCycle>>
    getMatchedRepresentativeCycles(const size_t &dim, const size_t &index);
    variant<dim1::RepresentativeCycle, dim2::RepresentativeCycle,
            dim3::RepresentativeCycle>
    getUnmatchedRepresentativeCycle(const uint8_t &input, const size_t &dim,
                                    const size_t &index);
    tuple<vector<dim3::RepresentativeCycle>, vector<dim3::RepresentativeCycle>>
    computeAllRepresentativeCycles(const int input, const int dim,
                                   bool computeMatchedCycles,
                                   bool computeUnmatchedCycles);
    vector<index_t> shape;

  private:
    optional<variant<dim1::BettiMatching, dim2::BettiMatching,
                     dim3::BettiMatching, dimN::BettiMatching>>
        dimensionSpecificBettiMatching;
    size_t dimension;
    bool computed;
    bool computedPairsInput0;
};
