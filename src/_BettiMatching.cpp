#include "utils.h"
#include "config.h"
#include "data_structures.h"
#include "BettiMatching.h"

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <future>
#include <iostream>
#include <iterator>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <string>
#include <tuple>
#include <utility>

using namespace std;
namespace py = pybind11;

struct BettiMatchingResult {
    py::array_t<int64_t> predictionMatchesBirthCoordinates;
    py::array_t<int64_t> predictionMatchesDeathCoordinates;
    py::array_t<int64_t> targetMatchesBirthCoordinates;
    py::array_t<int64_t> targetMatchesDeathCoordinates;
    py::array_t<int64_t> predictionUnmatchedBirthCoordinates;
    py::array_t<int64_t> predictionUnmatchedDeathCoordinates;
    py::array_t<int64_t> numMatchesByDim;
    py::array_t<int64_t> numUnmatchedPredictionByDim;
};

string repr_vector(const vector<index_t> shape, std::tuple<string, string> parentheses = make_tuple("(", ")"), string separator = ", ") {
    stringstream out_stream;
    out_stream << std::get<0>(parentheses);
    for_each(shape.begin(), std::prev(shape.end()), [&out_stream, &shape, &separator](auto size)
                { out_stream << std::to_string(size) << separator; });
    out_stream << shape[shape.size() - 1] << std::get<1>(parentheses);
    return out_stream.str();
};

std::tuple<vector<vector<VoxelMatch>>, vector<vector<VoxelPair>>, vector<vector<VoxelPair>>> computeMatchingFromNumpyArrays(
    py::array_t<value_t> &input0, py::array_t<value_t> &input1)
{
    vector<index_t> shape0(input0.shape(), input0.shape() + input0.ndim());
    vector<index_t> shape1(input1.shape(), input1.shape() + input1.ndim());
    if (shape0 != shape1) {
        throw invalid_argument("The shapes of the input volumes must agree. Got " + repr_vector(shape0) + " and " + repr_vector(shape1));
    }
    vector<value_t> input0Vector(input0.mutable_data(), input0.mutable_data() + input0.size());
    vector<value_t> input1Vector(input1.mutable_data(), input1.mutable_data() + input1.size());

    Config config;
    BettiMatching bettiMatching(std::move(input0Vector), std::move(input1Vector), std::move(shape0), std::move(config));
    bettiMatching.computeMatching();
    return bettiMatching.getMatching();
};

PYBIND11_MODULE(betti_matching, m) {
    py::class_<BettiMatching>(m, "BettiMatching")
        .def(py::init([](py::array_t<value_t>& input0, py::array_t<value_t>& input1) {
                vector<index_t> shape0(input0.shape(), input0.shape() + input0.ndim());
                vector<index_t> shape1(input1.shape(), input1.shape() + input1.ndim());
                if (shape0 != shape1) {
                    throw invalid_argument("The shapes of the input volumes must agree. Got " + repr_vector(shape0)
                                            + " and " + repr_vector(shape1));
                }
                vector<value_t> input0Vector(input0.mutable_data(), input0.mutable_data() + input0.size());
                vector<value_t> input1Vector(input1.mutable_data(), input1.mutable_data() + input1.size());
                Config config;
                return BettiMatching(std::move(input0Vector), std::move(input1Vector), std::move(shape0), std::move(config));
            }))

        .def(py::init([](string input0_path, string input1_path) {
                vector<value_t> input0Vector;
                vector<value_t> input1Vector;
                vector<index_t> shape0;
                vector<index_t> shape1;
                readImage(input0_path, NUMPY, input0Vector, shape0);
                readImage(input1_path, NUMPY, input1Vector, shape1);
                if (shape0 != shape1) {
                    throw invalid_argument("The shapes of the tree input volumes must agree. Got " + repr_vector(shape0)
                                            + " and " + repr_vector(shape1));
                }
                Config config;
                return BettiMatching(std::move(input0Vector), std::move(input1Vector), std::move(shape0), std::move(config));
            }))

        .def_readwrite("shape", &BettiMatching::shape)

        .def("compute_matching", &BettiMatching::computeMatching)

        .def("print", &BettiMatching::printResult)

        .def("get_matching", &BettiMatching::getMatching)

        .def("get_matched_cycles", &BettiMatching::getMatchedRepresentativeCycles)

        .def("get_unmatched_cycle", &BettiMatching::getUnmatchedRepresentativeCycle);


    m.def("compute_matching", &computeMatchingFromNumpyArrays);

    m.def("compute_matching", [](
        vector<py::array_t<value_t>> inputs0,
        vector<py::array_t<value_t>> inputs1
    )
    {
        if (inputs0.size() != inputs1.size()) {
            throw invalid_argument("Different numbers of inputs where provided: " + std::to_string(inputs0.size()) + " (inputs0) and " + std::to_string(inputs1.size()) + " (inputs1)");
        }
        size_t batchSize = inputs0.size();
        size_t numDimensions = inputs0[0].ndim();

        // Create one asynchronous task for Betti matching computation per input in the batch
        std::vector<std::future<std::tuple<vector<vector<VoxelMatch>>, vector<vector<VoxelPair>>, vector<vector<VoxelPair>>>>> resultFutures;
        for (int i = 0; i < inputs0.size(); i++) {
            resultFutures.push_back(std::async(computeMatchingFromNumpyArrays, std::ref(inputs0[i]), std::ref(inputs1[i])));
        }

        // Now block on all of them one at a time.
        std::vector<std::tuple<vector<vector<VoxelMatch>>, vector<vector<VoxelPair>>, vector<vector<VoxelPair>>>> pairResults;
        for (auto& future : resultFutures) {
            pairResults.push_back(future.get());
        }

        vector<BettiMatchingResult> arrayResults;

        // Write the results into numpy arrays
        for (auto &result : pairResults) {
            auto matchesByDimension = std::get<0>(result);
            auto predictionUnmatchedByDimension = std::get<1>(result);

            size_t numMatches = 0;
            size_t numPredictionUnmatched = 0;
            for (size_t d = 0; d < numDimensions; d++) {
                numMatches += matchesByDimension[d].size();
                numPredictionUnmatched += predictionUnmatchedByDimension[d].size();
            }

            // Shape (num matched in total, num dimensions)
            py::array_t<int64_t> predictionMatchesBirthCoordinates({numMatches, numDimensions}, {numDimensions * sizeof(int64_t), sizeof(int64_t)});
            py::array_t<int64_t> predictionMatchesDeathCoordinates({numMatches, numDimensions}, {numDimensions * sizeof(int64_t), sizeof(int64_t)});
            py::array_t<int64_t> targetMatchesBirthCoordinates({numMatches, numDimensions}, {numDimensions * sizeof(int64_t), sizeof(int64_t)});
            py::array_t<int64_t> targetMatchesDeathCoordinates({numMatches, numDimensions}, {numDimensions * sizeof(int64_t), sizeof(int64_t)});
            // Shape (num unmatched (prediction) in total, num dimensions)
            py::array_t<int64_t> predictionUnmatchedBirthCoordinates({numPredictionUnmatched, numDimensions}, {numDimensions * sizeof(int64_t), sizeof(int64_t)});
            py::array_t<int64_t> predictionUnmatchedDeathCoordinates({numPredictionUnmatched, numDimensions}, {numDimensions * sizeof(int64_t), sizeof(int64_t)});
            // Shape (num dimensions,)
            py::array_t<int64_t> numMatchesByDim({numDimensions}, {sizeof(int64_t)});
            py::array_t<int64_t> numUnmatchedPredictionByDim({numDimensions}, {sizeof(int64_t)});

            auto predictionMatchesBirthCoordinatesView = predictionMatchesBirthCoordinates.mutable_unchecked();
            auto predictionMatchesDeathCoordinatesView = predictionMatchesDeathCoordinates.mutable_unchecked();
            auto targetMatchesBirthCoordinatesView = targetMatchesBirthCoordinates.mutable_unchecked();
            auto targetMatchesDeathCoordinatesView = targetMatchesDeathCoordinates.mutable_unchecked();
            auto predictionUnmatchedBirthCoordinatesView = predictionUnmatchedBirthCoordinates.mutable_unchecked();
            auto predictionUnmatchedDeathCoordinatesView = predictionUnmatchedDeathCoordinates.mutable_unchecked();
            auto numMatchesByDimView = numMatchesByDim.mutable_unchecked();
            auto numUnmatchedPredictionByDimView = numUnmatchedPredictionByDim.mutable_unchecked();

            int i = 0;
            int currentDimension = 0;
            for (auto &matchesInDimension : matchesByDimension) {
                for (auto &match : matchesInDimension) {
                    for (int d = 0; d < numDimensions; d++) {
                        predictionMatchesBirthCoordinatesView(i, d) = match.pair0.birth[d];
                        predictionMatchesDeathCoordinatesView(i, d) = match.pair0.death[d];
                        targetMatchesBirthCoordinatesView(i, d) = match.pair1.birth[d];
                        targetMatchesDeathCoordinatesView(i, d) = match.pair1.death[d];
                    }
                    i += 1;
                }
                numMatchesByDimView(currentDimension) = matchesInDimension.size();
                currentDimension++;
            }
            i = 0;
            currentDimension = 0;
            for (auto &unmatchedInDimension : predictionUnmatchedByDimension) {
                for (auto &unmatched : unmatchedInDimension) {
                    for (int d = 0; d < numDimensions; d++) {
                        predictionUnmatchedBirthCoordinatesView(i, d) = unmatched.birth[d];
                        predictionUnmatchedDeathCoordinatesView(i, d) = unmatched.death[d];
                    }
                    i += 1;
                }
                numUnmatchedPredictionByDimView(currentDimension) = unmatchedInDimension.size();
                currentDimension += 1;
            }

            arrayResults.emplace_back(BettiMatchingResult{
                std::move(predictionMatchesBirthCoordinates),
                std::move(predictionMatchesDeathCoordinates),
                std::move(targetMatchesBirthCoordinates),
                std::move(targetMatchesDeathCoordinates),
                std::move(predictionUnmatchedBirthCoordinates),
                std::move(predictionUnmatchedDeathCoordinates),
                std::move(numMatchesByDim),
                std::move(numUnmatchedPredictionByDim)});
        }
        return arrayResults;
    });



    py::class_<VoxelMatch>(m, "VoxelMatch")
        .def_readonly("pair0", &VoxelMatch::pair0)
        .def_readonly("pair1", &VoxelMatch::pair1)
        .def("__repr__", [](VoxelMatch &self) {
            return "VoxelMatch(pair0=(birth=" + repr_vector(self.pair0.birth) + ", death=" + repr_vector(self.pair0.death)
                    + "), pair1=(birth=" + repr_vector(self.pair1.birth) + ", death=" + repr_vector(self.pair1.death) + ")"; });


    py::class_<VoxelPair>(m, "VoxelPair")
        .def_readonly("birth", &VoxelPair::birth)
        .def_readonly("death", &VoxelPair::death)
        .def("__repr__", [](VoxelPair &self) {
          return "VoxelPair(birth=" + repr_vector(self.birth) +
                 ", death=" + repr_vector(self.death) + ")";
        });

    auto resultTypesModule = m.def_submodule("return_types", "Return types for betti_matching");

    py::class_<BettiMatchingResult>(resultTypesModule, "BettiMatchingResult")
        .def_readonly("prediction_matches_birth_coordinates", &BettiMatchingResult::predictionMatchesBirthCoordinates)
        .def_readonly("prediction_matches_death_coordinates", &BettiMatchingResult::predictionMatchesDeathCoordinates)
        .def_readonly("target_matches_birth_coordinates", &BettiMatchingResult::targetMatchesBirthCoordinates)
        .def_readonly("target_matches_death_coordinates", &BettiMatchingResult::targetMatchesDeathCoordinates)
        .def_readonly("prediction_unmatched_birth_coordinates", &BettiMatchingResult::predictionUnmatchedBirthCoordinates)
        .def_readonly("prediction_unmatched_death_coordinates", &BettiMatchingResult::predictionUnmatchedDeathCoordinates)
        .def_readonly("num_matches_by_dim", &BettiMatchingResult::numMatchesByDim)
        .def_readonly("num_unmatched_prediction_by_dim", &BettiMatchingResult::numUnmatchedPredictionByDim)
        .def("__repr__", [](BettiMatchingResult &self) {
            auto reprMemberArray = [](string name, py::array_t<int64_t>& array) {
                return name + "=" + repr_vector(std::vector<index_t>(array.shape(), array.shape() + array.ndim()), make_tuple("[", "]"));
            };
            return "BettiMatchingResult(" + (
                reprMemberArray("prediction_matches_birth_coordinates", self.predictionMatchesBirthCoordinates) + ", " +
                reprMemberArray("prediction_matches_death_coordinates", self.predictionMatchesDeathCoordinates) + ", " +
                reprMemberArray("target_matches_birth_coordinates", self.targetMatchesBirthCoordinates) + ", " +
                reprMemberArray("target_matches_death_coordinates", self.targetMatchesDeathCoordinates) + ", " +
                reprMemberArray("prediction_unmatched_birth_coordinates", self.predictionUnmatchedBirthCoordinates) + ", " +
                reprMemberArray("prediction_unmatched_death_coordinates", self.predictionUnmatchedDeathCoordinates) + ", " +
                reprMemberArray("num_matches_by_dim", self.numMatchesByDim) + ", " +
                reprMemberArray("num_unmatched_prediction_by_dim", self.numUnmatchedPredictionByDim)
            ) + ")";
        });
}
