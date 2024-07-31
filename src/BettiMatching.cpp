#include "BettiMatching.h"
#include "config.h"
#include "data_structures.h"
#include "src_3D/data_structures.h"
#include "utils.h"

#include <algorithm>
#include <chrono>
#include <iostream>
#include <optional>
#include <stdexcept>
#include <vector>

using namespace std;
using namespace std::chrono;

BettiMatching::BettiMatching(vector<value_t> &&input0, vector<value_t> &&input1,
                             vector<index_t> &&_shape, Config &&config)
    : shape(_shape), dimension(shape.size()), computed(false),
      computedPairsInput0(false) {
#ifdef RUNTIME
    cout << "initializing BettiMatching ... ";
    auto start = high_resolution_clock::now();
#endif

    vector<value_t> comparison;
    transform(input0.begin(), input0.end(), input1.begin(),
              back_inserter(comparison),
              [](value_t a, value_t b) { return min(a, b); });

    switch (dimension) {
    case 0:
        throw std::invalid_argument("Input must be at least 1-dimensional.");

    case 1:
        dimensionSpecificBettiMatching.emplace(dim1::BettiMatching(
            std::move(input0), std::move(input1), std::move(comparison),
            std::move(shape), std::move(config)));
        break;

    case 2:
        dimensionSpecificBettiMatching.emplace(dim2::BettiMatching(
            std::move(input0), std::move(input1), std::move(comparison),
            std::move(shape), std::move(config)));
        break;

    case 3:
        dimensionSpecificBettiMatching.emplace(dim3::BettiMatching(
            std::move(input0), std::move(input1), std::move(comparison),
            std::move(shape), std::move(config)));
        break;

    default:
        dimensionSpecificBettiMatching.emplace(dimN::BettiMatching(
            std::move(input0), std::move(input1), std::move(comparison),
            std::move(shape), std::move(config)));
    }

#ifdef RUNTIME
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);
    cout << duration.count() << " ms" << endl << endl;
#endif
}

BettiMatching::BettiMatching(BettiMatching &&other)
    : shape(other.shape), dimensionSpecificBettiMatching(
                              std::move(other.dimensionSpecificBettiMatching)),
      dimension(other.dimension), computed(other.computed),
      computedPairsInput0(other.computedPairsInput0) {}

void BettiMatching::computeMatching() {
    if (computed) {
#ifdef RUNTIME
        cout << "already computed ..." << endl;
#endif
        return;
    }
    if (computedPairsInput0) {
        throw runtime_error("computeMatching cannot be called after "
                            "computePairsInput0 has been called");
    }

#ifdef RUNTIME
    cout << "computing Betti Matching ..." << endl;
    auto start = high_resolution_clock::now();
#endif

    switch (dimension) {
    case 1: {
        dim1::BettiMatching &bettiMatching = std::get<dim1::BettiMatching>(
            dimensionSpecificBettiMatching.value());
        bettiMatching.computeMatching();
        bettiMatching.computeVoxels();
        break;
    }

    case 2: {
        dim2::BettiMatching &bettiMatching = std::get<dim2::BettiMatching>(
            dimensionSpecificBettiMatching.value());
        bettiMatching.computeMatching();
        bettiMatching.computeVoxels();
        break;
    }

    case 3: {
        dim3::BettiMatching &bettiMatching = std::get<dim3::BettiMatching>(
            dimensionSpecificBettiMatching.value());
        bettiMatching.computeMatching();
        bettiMatching.computeVoxels();
        break;
    }

    default: {
        dimN::BettiMatching &bettiMatching = std::get<dimN::BettiMatching>(
            dimensionSpecificBettiMatching.value());
        bettiMatching.computeMatching();
        bettiMatching.computeVoxels();
        break;
    }
    }

    computed = true;

#ifdef RUNTIME
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);
    cout << "total runtime: " << duration.count() << " ms" << endl << endl;
#endif
}

vector<vector<VoxelPair>> BettiMatching::computePairsInput0() {
    if (computed || computedPairsInput0) {
        throw runtime_error("computePairsInput0 can only be called once and if "
                            "computeMatching has not been called before");
    }

    vector<vector<VoxelPair>> pairsInput0;
    switch (dimension) {
    case 1: {
        dim1::BettiMatching &bettiMatching = std::get<dim1::BettiMatching>(
            dimensionSpecificBettiMatching.value());
        pairsInput0 = bettiMatching.computePairsInput0();
        break;
    }

    case 2: {
        dim2::BettiMatching &bettiMatching = std::get<dim2::BettiMatching>(
            dimensionSpecificBettiMatching.value());
        pairsInput0 = bettiMatching.computePairsInput0();
        break;
    }

    case 3: {
        dim3::BettiMatching &bettiMatching = std::get<dim3::BettiMatching>(
            dimensionSpecificBettiMatching.value());
        pairsInput0 = bettiMatching.computePairsInput0();
        break;
    }

    default: {
        dimN::BettiMatching &bettiMatching = std::get<dimN::BettiMatching>(
            dimensionSpecificBettiMatching.value());
        pairsInput0 = bettiMatching.computePairsInput0();
        break;
    }
    }

    computedPairsInput0 = true;

    return pairsInput0;
}

void BettiMatching::printResult() {
    if (!computed) {
        throw std::runtime_error("Betti Matching not computed yet");
    }
    switch (dimension) {
    case 1: {
        dim1::BettiMatching &bettiMatching = std::get<dim1::BettiMatching>(
            dimensionSpecificBettiMatching.value());
        bettiMatching.printResult();
        break;
    }

    case 2: {
        dim2::BettiMatching &bettiMatching = std::get<dim2::BettiMatching>(
            dimensionSpecificBettiMatching.value());
        bettiMatching.printResult();
        break;
    }

    case 3: {
        dim3::BettiMatching &bettiMatching = std::get<dim3::BettiMatching>(
            dimensionSpecificBettiMatching.value());
        bettiMatching.printResult();
        break;
    }

    default: {
        dimN::BettiMatching &bettiMatching = std::get<dimN::BettiMatching>(
            dimensionSpecificBettiMatching.value());
        bettiMatching.printResult();
        break;
    }
    }
}

std::tuple<vector<vector<VoxelMatch>>, vector<vector<VoxelPair>>,
           vector<vector<VoxelPair>>>
BettiMatching::getMatching() {
    if (!computed) {
        throw std::runtime_error("Betti Matching not computed yet");
    }
    switch (dimension) {
    case 1: {
        dim1::BettiMatching &bettiMatching = std::get<dim1::BettiMatching>(
            dimensionSpecificBettiMatching.value());
        return {vector<vector<VoxelMatch>>{bettiMatching.matched},
                vector<vector<VoxelPair>>{bettiMatching.unmatched0},
                vector<vector<VoxelPair>>{bettiMatching.unmatched1}};
    }

    case 2: {
        dim2::BettiMatching &bettiMatching = std::get<dim2::BettiMatching>(
            dimensionSpecificBettiMatching.value());
        return {bettiMatching.matched, bettiMatching.unmatched0,
                bettiMatching.unmatched1};
    }

    case 3: {
        dim3::BettiMatching &bettiMatching = std::get<dim3::BettiMatching>(
            dimensionSpecificBettiMatching.value());
        return {bettiMatching.matched, bettiMatching.unmatched0,
                bettiMatching.unmatched1};
    }

    default: {
        dimN::BettiMatching &bettiMatching = std::get<dimN::BettiMatching>(
            dimensionSpecificBettiMatching.value());
        return {bettiMatching.matched, bettiMatching.unmatched0,
                bettiMatching.unmatched1};
    }
    }
}

variant<
    tuple<vector<dim1::RepresentativeCycle>, vector<dim1::RepresentativeCycle>>,
    tuple<vector<dim2::RepresentativeCycle>, vector<dim2::RepresentativeCycle>>,
    tuple<vector<dim3::RepresentativeCycle>, vector<dim3::RepresentativeCycle>>>
BettiMatching::computeRepresentativeCycles(
    const int input, const int dim,
    const optional<vector<size_t>> &matchedPairsIndices,
    const optional<vector<size_t>> &unmatchedPairsIndices) {
    if (!computed) {
        throw std::runtime_error("Betti Matching not computed yet");
    }

#ifndef COMPUTE_COMPARISON
    if (input == 2) {
        throw runtime_error(
            "Computing representative cycles for the comparison image requires "
            "the COMPUTE_COMPARISON compilation option to be set");
    }
#endif

    switch (dimension) {
    case 1: {
        dim1::BettiMatching &bettiMatching = std::get<dim1::BettiMatching>(
            dimensionSpecificBettiMatching.value());
        return bettiMatching.computeRepresentativeCycles(
            input, matchedPairsIndices, unmatchedPairsIndices);
    }

    case 2: {
        dim2::BettiMatching &bettiMatching = std::get<dim2::BettiMatching>(
            dimensionSpecificBettiMatching.value());
        return bettiMatching.computeRepresentativeCycles(
            input, dim, matchedPairsIndices, unmatchedPairsIndices);
    }

    case 3: {
        dim3::BettiMatching &bettiMatching = std::get<dim3::BettiMatching>(
            dimensionSpecificBettiMatching.value());
        return bettiMatching.computeRepresentativeCycles(
            input, dim, matchedPairsIndices, unmatchedPairsIndices);
    }

    default: {
        throw runtime_error("Computing representative cycles is only supported "
                            "for 1D, 2D and 3D volumes");
    }
    }
}
