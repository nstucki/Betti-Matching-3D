#include "BettiMatching.h"
#include "data_structures.h"
#include <vector>
#include "utils.h"
#include "config.h"
#include <optional>
#include <stdexcept>

BettiMatching::BettiMatching(vector<value_t> input0, vector<value_t> input1, vector<value_t> comparison, vector<index_t> shape, Config &config) : dimension(shape.size()), dimensionSpecificBettiMatching(dim1::BettiMatching(input0, input1, comparison, shape, config))
{
    switch (dimension)
    {
    case 0:
        throw std::invalid_argument("Input must be at least 1-dimensional.");
    case 1:
        dimensionSpecificBettiMatching.emplace(dim1::BettiMatching(input0, input1, comparison, shape, config));
        break;
    case 2:
        dimensionSpecificBettiMatching.emplace(dim2::BettiMatching(input0, input1, comparison, shape, config));
        break;
    case 3:
        dimensionSpecificBettiMatching.emplace(dim3::BettiMatching(input0, input1, comparison, shape, config));
        break;
    default:
        dimensionSpecificBettiMatching.emplace(dimN::BettiMatching(input0, input1, comparison, shape, config));
    }
}
BettiMatching::BettiMatching(BettiMatching &&other) : dimension(other.dimension), dimensionSpecificBettiMatching(std::move(other.dimensionSpecificBettiMatching)) {}

std::tuple<vector<vector<VoxelMatch>>, vector<vector<VoxelPair>>, vector<vector<VoxelPair>>> BettiMatching::computeMatchingWithVoxels()
{
    switch (dimension)
    {
    case 1:
    {
        dim1::BettiMatching &bettiMatching = std::get<dim1::BettiMatching>(dimensionSpecificBettiMatching.value());
        bettiMatching.computeMatching();
        bettiMatching.computeVoxels();
        return {
            std::vector<std::vector<VoxelMatch>>{bettiMatching.matched},
            std::vector<std::vector<VoxelPair>>{bettiMatching.unmatched0},
            std::vector<std::vector<VoxelPair>>{bettiMatching.unmatched1}};
    }
    break;
    case 2:
    {
        dim2::BettiMatching &bettiMatching = std::get<dim2::BettiMatching>(dimensionSpecificBettiMatching.value());
        bettiMatching.computeMatching();
        bettiMatching.computeVoxels();
        return {bettiMatching.matched, bettiMatching.unmatched0, bettiMatching.unmatched1};
    }
    break;
    case 3:
    {
        dim3::BettiMatching &bettiMatching = std::get<dim3::BettiMatching>(dimensionSpecificBettiMatching.value());
        bettiMatching.computeMatching();
        bettiMatching.computeVoxels();
        return {bettiMatching.matched, bettiMatching.unmatched0, bettiMatching.unmatched1};
    }
    break;
    default:
    {
        dimN::BettiMatching bettiMatching = std::get<dimN::BettiMatching>(dimensionSpecificBettiMatching.value());
        bettiMatching.computeMatching();
        bettiMatching.computeVoxels();
        return {bettiMatching.matched, bettiMatching.unmatched0, bettiMatching.unmatched1};
    }
    }
}