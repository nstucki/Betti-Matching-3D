#include "BettiMatching.h"
#include "data_structures.h"
#include "utils.h"
#include "config.h"

#include <vector>
#include <optional>
#include <stdexcept>
#include <iostream>

using namespace std;
using namespace std::chrono;



BettiMatching::BettiMatching(vector<value_t>&& input0, vector<value_t>&& input1, vector<index_t>&& shape, Config&& config) : 
    dimension(shape.size()), computed(false) {
#ifdef RUNTIME
    cout << "initializing BettiMatching ... ";
	auto start = high_resolution_clock::now();
#endif

    vector<value_t> comparison;
    transform(input0.begin(), input0.end(), input1.begin(), back_inserter(comparison), [](value_t a, value_t b) { return min(a, b); });
    
    switch (dimension) {
        case 0:
            throw std::invalid_argument("Input must be at least 1-dimensional.");

        case 1:
            dimensionSpecificBettiMatching.emplace(dim1::BettiMatching(std::move(input0), std::move(input1), std::move(comparison), 
                                                                        std::move(shape), std::move(config)));
            break;

        case 2:
            dimensionSpecificBettiMatching.emplace(dim2::BettiMatching(std::move(input0), std::move(input1), std::move(comparison), 
                                                                        std::move(shape), std::move(config)));
            break;

        case 3:
            dimensionSpecificBettiMatching.emplace(dim3::BettiMatching(std::move(input0), std::move(input1), std::move(comparison),
                                                                        std::move(shape), std::move(config)));
            break;

        default:
            dimensionSpecificBettiMatching.emplace(dimN::BettiMatching(std::move(input0), std::move(input1), std::move(comparison),
                                                                        std::move(shape), std::move(config)));
    }

#ifdef RUNTIME
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<milliseconds>(stop - start);
	cout << duration.count() << " ms" << endl << endl;
#endif
}


//BettiMatching::BettiMatching(BettiMatching&& other) : 
//    dimension(other.dimension), dimensionSpecificBettiMatching(std::move(other.dimensionSpecificBettiMatching)) {}


void BettiMatching::computeMatching() {
    if (computed) {
#ifdef RUNTIME
        cout << "already computed ..." << endl;
#endif        
        return;
    }

#ifdef RUNTIME
	cout << "computing Betti Matching ..." << endl;
    auto start = high_resolution_clock::now();
#endif

    switch (dimension) {
        case 1: {
            dim1::BettiMatching &bettiMatching = std::get<dim1::BettiMatching>(dimensionSpecificBettiMatching.value());
            bettiMatching.computeMatching();
            bettiMatching.computeVoxels();
            break;
        }

        case 2: {
            dim2::BettiMatching &bettiMatching = std::get<dim2::BettiMatching>(dimensionSpecificBettiMatching.value());
            bettiMatching.computeMatching();
            bettiMatching.computeVoxels();
            break;
        }

        case 3: {
            dim3::BettiMatching &bettiMatching = std::get<dim3::BettiMatching>(dimensionSpecificBettiMatching.value());
            bettiMatching.computeMatching();
            bettiMatching.computeVoxels();
            break;
        }

        default: {
            dimN::BettiMatching &bettiMatching = std::get<dimN::BettiMatching>(dimensionSpecificBettiMatching.value());
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


std::tuple<vector<vector<VoxelMatch>>, vector<vector<VoxelPair>>, vector<vector<VoxelPair>>> BettiMatching::getMatching() {
    switch (dimension) {
        case 1: {
            dim1::BettiMatching& bettiMatching = std::get<dim1::BettiMatching>(dimensionSpecificBettiMatching.value());
            return { vector<vector<VoxelMatch>> {bettiMatching.matched}, 
                        vector<vector<VoxelPair>> {bettiMatching.unmatched0},
                        vector<vector<VoxelPair>> {bettiMatching.unmatched1} };
        }

        case 2: {
            dim2::BettiMatching& bettiMatching = std::get<dim2::BettiMatching>(dimensionSpecificBettiMatching.value());
            return { bettiMatching.matched, bettiMatching.unmatched0, bettiMatching.unmatched1 };
        }

        case 3: {
            dim3::BettiMatching& bettiMatching = std::get<dim3::BettiMatching>(dimensionSpecificBettiMatching.value());
            return { bettiMatching.matched, bettiMatching.unmatched0, bettiMatching.unmatched1 };
        }

        default: {
            dimN::BettiMatching& bettiMatching = std::get<dimN::BettiMatching>(dimensionSpecificBettiMatching.value());
            return { bettiMatching.matched, bettiMatching.unmatched0, bettiMatching.unmatched1 };
        }
    }
}


tuple<set<vector<index_t>>, set<vector<index_t>>> BettiMatching::getMatchedRepresentativeCycle(const size_t& dim, const size_t& index) {
    tuple<set<vector<index_t>>, set<vector<index_t>>> repCycles;
    switch (dimension) {
        case 1: {
            return repCycles;
        }

        case 2: {
            return repCycles;
        }

        case 3: {
            dim3::BettiMatching& bettiMatching = std::get<dim3::BettiMatching>(dimensionSpecificBettiMatching.value());
            repCycles = bettiMatching.getMatchedRepresentativeCycle(dim, index);
            return repCycles;
        }

        default: {
            return repCycles;
        }
    }
}


set<vector<index_t>> BettiMatching::getUnmatchedRepresentativeCycle(const size_t& dim, const size_t& index, const uint8_t& input) {
    set<vector<index_t>> repCycle;
    switch (dimension) {
        case 1: {
            return repCycle;
        }

        case 2: {
            return repCycle;
        }

        case 3: {
            dim3::BettiMatching& bettiMatching = std::get<dim3::BettiMatching>(dimensionSpecificBettiMatching.value());
            repCycle = bettiMatching.getUnmatchedRepresentativeCycle(dim, index, input);
            return repCycle;
        }

        default: {
            return repCycle;
        }
    }
}


void BettiMatching::printResult() {
    switch (dimension) {
        case 1: {
            dim1::BettiMatching& bettiMatching = std::get<dim1::BettiMatching>(dimensionSpecificBettiMatching.value());
            bettiMatching.printResult();
            break;
        }

        case 2: {
            dim2::BettiMatching& bettiMatching = std::get<dim2::BettiMatching>(dimensionSpecificBettiMatching.value());
            bettiMatching.printResult();
            break;
        }

        case 3: {
            dim3::BettiMatching& bettiMatching = std::get<dim3::BettiMatching>(dimensionSpecificBettiMatching.value());
            bettiMatching.printResult();
            break;
        }

        default: {
            dimN::BettiMatching& bettiMatching = std::get<dimN::BettiMatching>(dimensionSpecificBettiMatching.value());
            bettiMatching.printResult();
            break;
        }
    }
}