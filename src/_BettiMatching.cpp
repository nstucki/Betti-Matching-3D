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

typedef py::array InputVolume;
typedef py::array_t<value_t, py::array::c_style> TypedInputVolume;
typedef std::tuple<vector<vector<VoxelMatch>>, vector<vector<VoxelPair>>, vector<vector<VoxelPair>>> BettiMatchingPairsResult;

struct BettiMatchingResult {
    py::array_t<int64_t> predictionMatchesBirthCoordinates;
    py::array_t<int64_t> predictionMatchesDeathCoordinates;
    py::array_t<int64_t> targetMatchesBirthCoordinates;
    py::array_t<int64_t> targetMatchesDeathCoordinates;
    py::array_t<int64_t> predictionUnmatchedBirthCoordinates;
    py::array_t<int64_t> predictionUnmatchedDeathCoordinates;
    py::array_t<int64_t> numMatchesByDim;
    py::array_t<int64_t> numUnmatchedPredictionByDim;
    // Optional return values if include_target_unmatched_pairs flag is set (not needed for Betti matching training loss)
    std::optional<py::array_t<int64_t>> targetUnmatchedBirthCoordinates;
    std::optional<py::array_t<int64_t>> targetUnmatchedDeathCoordinates;
    std::optional<py::array_t<int64_t>> numUnmatchedTargetByDim;
};
struct BarcodeResult {
    py::array_t<int64_t> birthCoordinates;
    py::array_t<int64_t> deathCoordinates;
    py::array_t<int64_t> numPairsByDim;
};

string repr_vector(const vector<index_t> shape, std::tuple<string, string> parentheses = make_tuple("(", ")"), string separator = ", ") {
    stringstream out_stream;
    out_stream << std::get<0>(parentheses);
    for_each(shape.begin(), std::prev(shape.end()), [&out_stream, &shape, &separator](auto size)
                { out_stream << std::to_string(size) << separator; });
    out_stream << shape[shape.size() - 1] << std::get<1>(parentheses);
    return out_stream.str();
};

BettiMatchingPairsResult computeMatchingFromNumpyArrays(TypedInputVolume &input0, TypedInputVolume &input1)
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

vector<vector<VoxelPair>> computeBarcodeFromNumpyArrays(TypedInputVolume &input0)
{
    vector<index_t> shape0(input0.shape(), input0.shape() + input0.ndim());
    vector<value_t> input0Vector(input0.mutable_data(), input0.mutable_data() + input0.size());
    vector<value_t> input0VectorCopy(input0.mutable_data(), input0.mutable_data() + input0.size());

    Config config;
    BettiMatching bettiMatching(std::move(input0Vector), std::move(input0VectorCopy), std::move(shape0), std::move(config));
    return bettiMatching.computePairsInput0();
};

BettiMatchingResult convertPairResultToArrayResult(BettiMatchingPairsResult &result, size_t numDimensions, bool includeTargetUnmatchedPairs);
BarcodeResult convertPairResultToBarcodeResult(vector<vector<VoxelPair>>& pairsByDimension);

vector<BettiMatchingResult> computeMatchingFromInputs(
        vector<InputVolume> &untypedInputs0,
        vector<InputVolume> &untypedInputs1,
        bool includeTargetUnmatchedPairs)
{
    // Validate and convert inputs
    vector<TypedInputVolume> inputs0;
    vector<TypedInputVolume> inputs1;
    if (untypedInputs0.size() != untypedInputs1.size()) {
        throw invalid_argument("Different numbers of inputs where provided: " + std::to_string(untypedInputs0.size()) + " (inputs0) and " + std::to_string(untypedInputs1.size()) + " (inputs1)");
    }
    for (size_t i = 0; i < untypedInputs0.size(); i++) {
        vector<index_t> shape0(untypedInputs0[i].shape(), untypedInputs0[i].shape() + untypedInputs0[i].ndim());
        vector<index_t> shape1(untypedInputs1[i].shape(), untypedInputs1[i].shape() + untypedInputs1[i].ndim());
        if (shape0 != shape1) {
            throw invalid_argument("The shapes of the input volumes must agree. Got " + repr_vector(shape0) + " and " + repr_vector(shape1));
        }
        inputs0.push_back(untypedInputs0[i].cast<TypedInputVolume>());
        inputs1.push_back(untypedInputs1[i].cast<TypedInputVolume>());
    }
    size_t batchSize = inputs0.size();
    size_t numDimensions = inputs0[0].ndim();

    // Create one asynchronous task for Betti matching computation per input in the batch
    std::vector<std::future<BettiMatchingPairsResult>> resultFutures;
    for (int i = 0; i < inputs0.size(); i++) {
        resultFutures.push_back(std::async(computeMatchingFromNumpyArrays, std::ref(inputs0[i]), std::ref(inputs1[i])));
    }

    // Now block on all of them one at a time.
    std::vector<BettiMatchingPairsResult> pairResults;
    for (auto& future : resultFutures) {
        pairResults.push_back(future.get());
    }

    vector<BettiMatchingResult> arrayResults;

    // Write the results into numpy arrays
    for (auto &result : pairResults) {
        auto arrayResult = convertPairResultToArrayResult(result, numDimensions, includeTargetUnmatchedPairs);
        arrayResults.emplace_back(arrayResult);
    }
    return arrayResults;
}

vector<BarcodeResult> computeBarcodeFromInputs(vector<InputVolume> &untypedInputs) {
    vector<TypedInputVolume> inputs;
    for (auto &untypedInput : untypedInputs) {
        inputs.push_back(untypedInput.cast<TypedInputVolume>());
    }
    size_t batchSize = inputs.size();
    size_t numDimensions = inputs[0].ndim();

    std::vector<vector<vector<VoxelPair>>> pairResults;
    // Create one asynchronous task for barcode computation per input in the batch
    std::vector<std::future<vector<vector<VoxelPair>>>> resultFutures;
    for (int i = 0; i < inputs.size(); i++) {
        resultFutures.push_back(std::async(computeBarcodeFromNumpyArrays, std::ref(inputs[i])));
    }

    // Now block on all of them one at a time.
    for (auto& future : resultFutures) {
        pairResults.push_back(future.get());
    }

    // Write the results into numpy arrays
    vector<BarcodeResult> arrayResults;
    for (auto &result : pairResults) {
        arrayResults.emplace_back(convertPairResultToBarcodeResult(result));
    }
    return arrayResults;
}

BettiMatchingResult convertPairResultToArrayResult(BettiMatchingPairsResult &result, size_t numDimensions, bool includeTargetUnmatchedPairs) {
    auto& matchesByDimension = std::get<0>(result);
    auto& predictionUnmatchedByDimension = std::get<1>(result);
    auto& targetUnmatchedByDimension = std::get<2>(result);

    size_t numMatches = 0;
    size_t numPredictionUnmatched = 0;
    size_t numTargetUnmatched = 0;
    for (size_t d = 0; d < numDimensions; d++) {
        numMatches += matchesByDimension[d].size();
        numPredictionUnmatched += predictionUnmatchedByDimension[d].size();
        numTargetUnmatched += targetUnmatchedByDimension[d].size();
    }

    // Shape: (num matched in total, num dimensions)
    py::array_t<int64_t> predictionMatchesBirthCoordinates({numMatches, numDimensions}, {numDimensions * sizeof(int64_t), sizeof(int64_t)});
    py::array_t<int64_t> predictionMatchesDeathCoordinates({numMatches, numDimensions}, {numDimensions * sizeof(int64_t), sizeof(int64_t)});
    py::array_t<int64_t> targetMatchesBirthCoordinates({numMatches, numDimensions}, {numDimensions * sizeof(int64_t), sizeof(int64_t)});
    py::array_t<int64_t> targetMatchesDeathCoordinates({numMatches, numDimensions}, {numDimensions * sizeof(int64_t), sizeof(int64_t)});
    // Shape: (num unmatched (prediction) in total, num dimensions)
    py::array_t<int64_t> predictionUnmatchedBirthCoordinates({numPredictionUnmatched, numDimensions}, {numDimensions * sizeof(int64_t), sizeof(int64_t)});
    py::array_t<int64_t> predictionUnmatchedDeathCoordinates({numPredictionUnmatched, numDimensions}, {numDimensions * sizeof(int64_t), sizeof(int64_t)});
    // Shape: (num unmatched (target) in total, num dimensions)
    py::array_t<int64_t> targetUnmatchedBirthCoordinates;
    py::array_t<int64_t> targetUnmatchedDeathCoordinates;
    if (includeTargetUnmatchedPairs) {
        targetUnmatchedBirthCoordinates = py::array_t<int64_t>({numTargetUnmatched, numDimensions}, {numDimensions * sizeof(int64_t), sizeof(int64_t)});
        targetUnmatchedDeathCoordinates = py::array_t<int64_t>({numTargetUnmatched, numDimensions}, {numDimensions * sizeof(int64_t), sizeof(int64_t)});
    }
    // Shape: (num dimensions,)
    py::array_t<int64_t> numMatchesByDim({numDimensions}, {sizeof(int64_t)});
    py::array_t<int64_t> numUnmatchedPredictionByDim({numDimensions}, {sizeof(int64_t)});
    py::array_t<int64_t> numUnmatchedTargetByDim;
    if (includeTargetUnmatchedPairs) {
        numUnmatchedTargetByDim = py::array_t<int64_t>({numDimensions}, {sizeof(int64_t)});
    }

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
            i++;
        }
        numMatchesByDimView(currentDimension++) = matchesInDimension.size();
    }
    i = 0;
    currentDimension = 0;
    for (auto &unmatchedInDimension : predictionUnmatchedByDimension) {
        for (auto &unmatched : unmatchedInDimension) {
            for (int d = 0; d < numDimensions; d++) {
                predictionUnmatchedBirthCoordinatesView(i, d) = unmatched.birth[d];
                predictionUnmatchedDeathCoordinatesView(i, d) = unmatched.death[d];
            }
            i++;
        }
        numUnmatchedPredictionByDimView(currentDimension++) = unmatchedInDimension.size();
    }
    if (includeTargetUnmatchedPairs) {
        i = 0;
        currentDimension = 0;
        auto targetUnmatchedBirthCoordinatesView = targetUnmatchedBirthCoordinates.mutable_unchecked();
        auto targetUnmatchedDeathCoordinatesView = targetUnmatchedDeathCoordinates.mutable_unchecked();
        auto numUnmatchedTargetByDimView = numUnmatchedTargetByDim.mutable_unchecked();
        for (auto &unmatchedInDimension : targetUnmatchedByDimension) {
            for (auto &unmatched : unmatchedInDimension) {
                for (int d = 0; d < numDimensions; d++) {
                    targetUnmatchedBirthCoordinatesView(i, d) = unmatched.birth[d];
                    targetUnmatchedDeathCoordinatesView(i, d) = unmatched.death[d];
                }
                i++;
            }
            numUnmatchedTargetByDimView(currentDimension++) = unmatchedInDimension.size();
        }
    }

    auto arrayResult = BettiMatchingResult{
        predictionMatchesBirthCoordinates=std::move(predictionMatchesBirthCoordinates),
        predictionMatchesDeathCoordinates=std::move(predictionMatchesDeathCoordinates),
        targetMatchesBirthCoordinates=std::move(targetMatchesBirthCoordinates),
        targetMatchesDeathCoordinates=std::move(targetMatchesDeathCoordinates),
        predictionUnmatchedBirthCoordinates=std::move(predictionUnmatchedBirthCoordinates),
        predictionUnmatchedDeathCoordinates=std::move(predictionUnmatchedDeathCoordinates),
        numMatchesByDim=std::move(numMatchesByDim),
        numUnmatchedPredictionByDim=std::move(numUnmatchedPredictionByDim)};
    if (includeTargetUnmatchedPairs) {
        arrayResult.targetUnmatchedBirthCoordinates = std::move(targetUnmatchedBirthCoordinates);
        arrayResult.targetUnmatchedDeathCoordinates = std::move(targetUnmatchedDeathCoordinates);
        arrayResult.numUnmatchedTargetByDim = std::move(numUnmatchedTargetByDim);
    }
    return arrayResult;
}

BarcodeResult convertPairResultToBarcodeResult(vector<vector<VoxelPair>>& pairsByDimension) {
    auto numDimensions = pairsByDimension.size();

    size_t numPairs = 0;
    for (size_t d = 0; d < numDimensions; d++) {
        numPairs += pairsByDimension[d].size();
    }

    // Shape: (num pairs, num dimensions)
    py::array_t<int64_t> birthCoordinates({numPairs, numDimensions}, {numDimensions * sizeof(int64_t), sizeof(int64_t)});
    py::array_t<int64_t> deathCoordinates({numPairs, numDimensions}, {numDimensions * sizeof(int64_t), sizeof(int64_t)});

    // Shape: (num dimensions,)
    py::array_t<int64_t> numPairsByDim({numDimensions}, {sizeof(int64_t)});

    auto birthCoordinatesView = birthCoordinates.mutable_unchecked();
    auto deathCoordinatesView = deathCoordinates.mutable_unchecked();
    auto numPairsByDimView = numPairsByDim.mutable_unchecked();

    int i = 0;
    int currentDimension = 0;
    for (auto &pairsInDimension : pairsByDimension) {
        for (auto &pair : pairsInDimension) {
            for (int d = 0; d < numDimensions; d++) {
                birthCoordinatesView(i, d) = pair.birth[d];
                deathCoordinatesView(i, d) = pair.death[d];
            }
            i++;
        }
        numPairsByDimView(currentDimension++) = pairsInDimension.size();
    }

    auto arrayResult = BarcodeResult{
        birthCoordinates=birthCoordinates,
        deathCoordinates=deathCoordinates,
        numPairsByDim=numPairsByDim
    };
    return arrayResult;
}

py::array_t<int64_t> coordinateListToArray(vector<tuple<uint64_t, uint64_t, uint64_t>>& coordinateList) {
    py::array_t<int64_t> coordinateArray({(int)coordinateList.size(), 3}, {3 * sizeof(int64_t), sizeof(int64_t)});
    auto coordinateArrayView = coordinateArray.mutable_unchecked();
    for (int i = 0; i < coordinateList.size(); i++) {
        auto& coordinates = coordinateList[i];
        coordinateArrayView(i, 0) = std::get<0>(coordinates);
        coordinateArrayView(i, 1) = std::get<1>(coordinates);
        coordinateArrayView(i, 2) = std::get<2>(coordinates);
    }
    return coordinateArray;
}

py::array_t<int64_t> coordinateListToArray(vector<vector<index_t>>& coordinateList, size_t numDimensions) {
    py::array_t<int64_t> coordinateArray({(int)coordinateList.size(), (int)numDimensions}, {numDimensions * sizeof(int64_t), sizeof(int64_t)});
    auto coordinateArrayView = coordinateArray.mutable_unchecked();
    for (int i = 0; i < coordinateList.size(); i++) {
        for (int j = 0; j < numDimensions; j++) {
            coordinateArrayView(i, j) = coordinateList[i][j];
        }
    }
    return coordinateArray;
}

PYBIND11_MODULE(betti_matching, m) {
    m.def("compute_matching",
        [](
            InputVolume &untypedInput0,
            InputVolume &untypedInput1,
            bool includeTargetUnmatchedPairs
        ) {
            vector<InputVolume> untypedInputs0 = {untypedInput0};
            vector<InputVolume> untypedInputs1 = {untypedInput1};
            return computeMatchingFromInputs(untypedInputs0, untypedInputs1, includeTargetUnmatchedPairs)[0];
        },
        py::arg("input0").noconvert(),
        py::arg("input1").noconvert(),
        py::arg("include_target_unmatched_pairs") = true,
        R"(
        compute_matching(input0, input1, include_target_unmatched_pairs=True)

        Compute the Betti matching between two input volumes.

        Parameters
        ----------
        input0 : np.ndarray
            The first input volume (the "prediction" in the machine
            learning context).
        input1 : np.ndarray
            The second input volume (the "target" in the machine
            learning context).
        include_target_unmatched_pairs : bool, optional
            Whether to include the unmatched pairs in the input1 (target) volume
            in the result. Default is True. Can be deactivated when the target
            unmatched pairs are not needed, such as in training with the Betti
            matching loss, where they do not contribute to the gradient.

        Returns
        -------
        result : BettiMatchingResult
            The result of the Betti matching between the two input volumes.
            See the documentation of `betti_matching.return_types.BettiMatchingResult`
            for details about the contained data.

        Example
        -------
        ```python
        a, b = np.random.rand(10, 10, 10), np.random.rand(10, 10, 10)
        result = betti_matching.compute_matching(a, b)
        a_dim0_matches_birth_coordinates = (
            result.prediction_matches_birth_coordinates[:result.num_matches_by_dim[0]])
        a_dim0_matches_death_coordinates = (
            result.prediction_matches_death_coordinates[:result.num_matches_by_dim[0]])
        a_dim0_matched_bars_lengths = (a[*a_dim0_matches_death_coordinates.T]
            - a[*a_dim0_matches_birth_coordinates.T])
        ```
        )"
    );

    m.def("compute_matching",
        &computeMatchingFromInputs,
        py::arg("inputs0").noconvert(),
        py::arg("inputs1").noconvert(),
        py::arg("include_target_unmatched_pairs") = true,
        R"(
        compute_matching(inputs0, inputs1, include_target_unmatched_pairs=True)

        Compute the Betti matching between two batches of input volumes in
        parallel. The Betti matching computation are parallelized using
        std::async.

        Parameters
        ----------
        inputs0 : list of np.ndarray
            The batch of first input volumes (the "predictions" in the machine
            learning context).
        inputs1 : list of np.ndarray
            The batch of second input volumes (the "targets" in the machine
            learning context).
        include_target_unmatched_pairs : bool, optional
            Whether to include the unmatched pairs in the input1 (target) volumes
            in the result. Default is True. Can be deactivated when the target
            unmatched pairs are not needed, such as in training with the Betti
            matching loss, where they do not contribute to the gradient.

        Returns
        -------
        results : list of BettiMatchingResult
            The results of the Betti matching between the corresponding pairs of
            input volumes.
            See the documentation of`betti_matching.return_types.BettiMatchingResult`
            for details about the contained data.

        Example
        -------
        ```python
        a1, b1 = np.random.rand(10, 10, 10), np.random.rand(10, 10, 10)
        a2, b2 = np.random.rand(8, 8, 8), np.random.rand(8, 8, 8)
        results = betti_matching.compute_matching([a1, a2], [b1, b2])
        num_matches_a1_b1 = results[0].num_matches_by_dim.sum()
        num_matches_a2_b2 = results[1].num_matches_by_dim.sum()
        ```
        )"
    );

    m.def("compute_barcode",
        [](
            InputVolume &untypedInput
        ) {
            vector<InputVolume> untypedInputs = {untypedInput};
            return computeBarcodeFromInputs(untypedInputs)[0];
        },
        py::arg("input").noconvert()
    );

    m.def("compute_barcode",
        &computeBarcodeFromInputs,
        py::arg("inputs").noconvert()
    );

    py::class_<BettiMatching>(m, "BettiMatching",
        R"(
        BettiMatching(input0, input1)

        Class for computing the Betti matching between two input volumes, and computing
        representative cycles. In comparison to `betti_matching.compute_matching()`,
        this class allows to cache the computed matching and use it to compute
        representative cycles for the matched and unmatched features.

        Representative cycles have the following format:
        - Each cycle is represented as a NumPy array of shape (len_cycle+1, d).
        - The entries are voxel coordinates which outline the cycle boundary.
            - Dimension 0: vertex (0-cube) coordinates of the cycle.
            - Dimension 1: endpoint vertex coordinates of edges (1-cubes) in the cycle
            - Dimension 2: corner vertex coordinates of faces (2-cubes) in the cycle
        - The last index coordinate is the voxel coordinate corresponding to the death
          cube of the persistence pair (a `dim+1`-dimensional cube).
        - Cycles in dimension 1 may have duplicate voxel coordinates.

        Parameters
        ----------
        input0 : np.ndarray
            The first input volume (the "prediction" in the machine learning context).
        input1 : np.ndarray
            The second input volume (the "target" in the machine learning context).

        Example
        -------
        ```python
        a, b = np.random.rand(10, 10, 10), np.random.rand(10, 10, 10)
        bm = betti_matching.BettiMatching(a, b)
        bm.compute_matching()
        result = bm.get_matching()
        ```
        )"
    )
        .def(py::init([](InputVolume& untypedInput0, InputVolume& untypedInput1) {
                auto input0 = untypedInput0.cast<TypedInputVolume>();
                auto input1 = untypedInput1.cast<TypedInputVolume>();
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
            }),
            py::arg("input0").noconvert(),
            py::arg("input1").noconvert()
        )

        .def_readwrite("shape", &BettiMatching::shape)

        .def("compute_matching", &BettiMatching::computeMatching,
            R"(
            compute_matching()

            Compute the Betti matching between the two input volumes.
            Skips the computation if it has already been performed previously on
            the same `BettiMatching` instance.
            )"
        )

        .def("get_matching", [](BettiMatching &self, bool includeTargetUnmatchedPairs)
            {
                BettiMatchingPairsResult result = self.getMatching();
                return convertPairResultToArrayResult(result, self.shape.size(), includeTargetUnmatchedPairs);
            },
            py::arg("include_target_unmatched_pairs") = true,
            R"(
            get_matching()

            Retrieve the computed matching. `compute_matching()` must be called
            before calling this method.

            Returns
            -------
            result : BettiMatchingResult
                The result of the Betti matching between the two input volumes.
                See the documentation of `betti_matching.return_types.BettiMatchingResult`
                for details about the contained data.
            )"
        )

        .def("print", &BettiMatching::printResult,
            R"(
            print()

            Print a summary of the computed matching. `compute_matching()` must be called
            before calling this method.
            )"
        )

        .def("compute_matched_representative_cycles",
            [](BettiMatching &self, size_t dim, size_t index) {
                auto result = self.getMatchedRepresentativeCycles(dim, index);

                py::array_t<int64_t> matchedCycle0 = coordinateListToArray(std::get<0>(result), self.shape.size());
                py::array_t<int64_t> matchedCycle1 = coordinateListToArray(std::get<1>(result), self.shape.size());

                return std::make_tuple(matchedCycle0, matchedCycle1);
            },
            py::arg("dim"),
            py::arg("index"),
            R"(
            compute_matched_representative_cycles(dim, index)

            Compute a representative cycle in each of the two input volumes for a given pair
            of matched persistence pairs.

            Computation time is on the order of computing the barcodes in dimension `dim` of
            the input volumes (depending on `index`, since small indices allow the
            computation to exit early).

            Parameters
            ----------
            dim : int
                The dimension of the persistence pair.
            index : int
                The index of the persistence pair, relative to the dimension.

            Returns
            -------
            matched_cycles : tuple of (np.ndarray, np.ndarray)
                The representative cycles in each input volume for the given pair of matched
                persistence pairs. The first array contains the representative cycle in the
                first input volume, and the second array contains the representative cycle
                in the second input volume.
                See the documentation of the `BettiMatching` class for details about the
                representative cycle format.

            Example
            -------
            ```python
            bm = betti_matching.BettiMatching(prediction, target)
            bm.compute_matching()
            first_matched_component_in_prediction, _ = (
                bm.compute_matched_representative_cycles(dim=0, index=0))
            cycle_image = np.zeros(prediction.shape)
            cycle_image[*first_matched_component_prediction.T] = 1
            ```
            )"
        )

        .def("compute_unmatched_representative_cycle",
            [](BettiMatching &self, uint8_t input, size_t dim, size_t index) {
                auto unmatchedCycle = self.getUnmatchedRepresentativeCycle(input, dim, index);
                return coordinateListToArray(unmatchedCycle, self.shape.size());
            },
            py::arg("input"),
            py::arg("dim"),
            py::arg("index"),
            R"(
            compute_unmatched_representative_cycle(input, dim, index)

            Compute a representative cycle in one of the input volumes for a given unmatched
            persistence pair.

            Computation time is on the order of computing the barcode in dimension `dim` of
            the input volume (depending on `index`, since small indices allow the
            computation to exit early).

            Parameters
            ----------
            input : int
                The index of the input volume for which to compute the unmatched cycle.
            dim : int
                The dimension of the persistence pair.
            index : int
                The index of the persistence pair, relative to the dimension.
            
            Returns
            -------
            unmatched_cycle : np.ndarray
                The representative cycle in the input volume for the given unmatched
                persistence pair.
                See the documentation of the `BettiMatching` class for details about the
                representative cycle format.
            
            Example
            -------
            ```python
            bm = betti_matching.BettiMatching(prediction, target)
            bm.compute_matching()
            first_unmatched_component_in_prediction = (
                bm.compute_unmatched_representative_cycle(input=0, dim=0, index=0))
            cycle_image = np.zeros(prediction.shape)
            cycle_image[*first_unmatched_component_prediction.T] = 1
            ```
            )"
        )

        .def("compute_all_representative_cycles",
            [](BettiMatching &self, int input, int dim, bool computeMatchedCycles, bool computeUnmatchedCycles) {
                auto matchedAndUnmatchedCycles = self.computeAllRepresentativeCycles(input, dim, computeMatchedCycles, computeUnmatchedCycles);
                auto matchedCycles = std::get<0>(matchedAndUnmatchedCycles);
                auto unmatchedCycles = std::get<1>(matchedAndUnmatchedCycles);

                std::optional<vector<py::array_t<int64_t>>> matchedCyclesArrays;
                std::optional<vector<py::array_t<int64_t>>> unmatchedCyclesArrays;
                
                if (computeMatchedCycles) {
                    matchedCyclesArrays = {{}};
                    std::transform(matchedCycles.begin(), matchedCycles.end(), back_inserter(*matchedCyclesArrays),
                        [](auto coords){return coordinateListToArray(coords);});
                }
                if (computeUnmatchedCycles) {
                    unmatchedCyclesArrays = {{}};
                    std::transform(unmatchedCycles.begin(), unmatchedCycles.end(), back_inserter(*unmatchedCyclesArrays),
                        [](auto coords){return coordinateListToArray(coords);});
                }

                return std::make_tuple(matchedCyclesArrays, unmatchedCyclesArrays);
            },
            py::arg("input"),
            py::arg("dim"),
            py::arg("compute_matched_cycles") = true,
            py::arg("compute_unmatched_cycles") = true,
            R"(
            compute_all_representative_cycles(input,
                                              dim,
                                              compute_matched_cycles=True,
                                              compute_unmatched_cycles=True)

            Currently only supported for 3D input volumes.

            Compute all (matched and unmatched) representative cycles in one of the input
            volumes for a given dimension. The method is more efficient than calling
            `compute_matched_representative_cycles` and
            `compute_unmatched_representative_cycle` for each index.

            The computation can be restricted to not include matched cycles, or not include
            unmatched cycles.

            Computation time:
            - for `dim==0` and `dim==2`: at least on the order of computing the barcode in
                dimension `dim` of the input volume. Additional overhead for saving the
                cycles applies.
            - for `dim==1`: utilizes the cache of the previous Betti matching computation,
                hence the computation time only comes from saving the cycles to arrays.

            Parameters
            ----------
            input : int
                The index of the input volume for which to compute the cycles.
            dim : int
                The dimension of the persistence pairs.
            compute_matched_cycles : bool, optional
                Whether to compute the matched cycles. Default is True.
            compute_unmatched_cycles : bool, optional
                Whether to compute the unmatched cycles. Default is True.

            Returns
            -------
            matched_and_unmatched_cycles : tuple of (
                    optional list of np.ndarray, optional list of np.ndarray)
                The matched and unmatched representative cycles in the input volume for the
                given dimension. The first list contains the matched cycles, and the second
                list contains the unmatched cycles.
                See the documentation of the `BettiMatching` class for details about the
                representative cycle format.
                If `compute_matched_cycles == False` or `compute_unmatched_cycles == False`,
                `None` is returned for the respective entry.
            
            Example
            -------
            ```python
            bm = betti_matching.BettiMatching(prediction, target)
            bm.compute_matching()
            matched_components_prediction, _ = bm.compute_all_representative_cycles(
                input=0, dim=0, compute_unmatched_cycles=False)
            cycle_image = np.zeros(prediction.shape)
            for i in reversed(range(len(matched_components_prediction))):
                cycle_image[*matched_components_prediction[i].T] = i+1
            ```
            )"
        );

    py::class_<dim3::Cube>(m, "Cube")
        .def("x", &dim3::Cube::x)
        .def("y", &dim3::Cube::y)
        .def("z", &dim3::Cube::z)
        .def("type", &dim3::Cube::type)
        .def("__repr__", [](dim3::Cube& self) {
          return "dim3::Cube(x=" + std::to_string(self.x()) + ", y=" + std::to_string(self.y()) + ", z=" + std::to_string(self.z()) + ")";
        });

    auto resultTypesModule = m.def_submodule("return_types", "Return types for betti_matching");

    py::class_<BettiMatchingResult>(resultTypesModule, "BettiMatchingResult",
        R"(
        BettiMatchingResult

        Holds the result of the Betti matching between two arrays (called "prediction"
        and "target"). The result contains the coordinates of the birth and death voxels
        for each matched and unmatched persistence pair, represented as NumPy arrays.

        Each coordinate array contains the birth or death coordinates of persistence
        pairs in all dimensions: first of the 0-dim. pairs, then of the 1-dim. pairs,
        and so on. The starting index of the persistence pairs birth/death coordinates
        in dimension i can be recovered using `np.cumsum()` on `num_matches_by_dim`,
        `num_unmatched_prediction_by_dim`, and `num_unmatched_target_by_dim`,
        respectively.

        Attributes
        ----------
        prediction_matches_birth_coordinates : numpy.ndarray
        prediction_matches_death_coordinates : numpy.ndarray
        target_matches_birth_coordinates : numpy.ndarray
        target_matches_death_coordinates : numpy.ndarray
            Arrays of shape (n_matched, d). The birth and death coordinates of matched
            persistence pairs in the prediction and target, respectively.
        prediction_unmatched_birth_coordinates : numpy.ndarray
        prediction_unmatched_death_coordinates : numpy.ndarray
            Arrays of shape (n_unmatched_prediction, d). The birth and death coordinates
            of unmatched persistence pairs in the prediction.
        target_unmatched_birth_coordinates : numpy.ndarray, optional
        target_unmatched_death_coordinates : numpy.ndarray, optional
            Arrays of shape (n_unmatched_target, d). The birth and death
            coordinates of unmatched persistence pairs in the target. Only present if
            the `include_target_unmatched_pairs` flag was set to `True` in the
            Betti matching computation (can be turned off since the target unmatched
            pairs do not contribute to the gradient when training with the Betti
            matching loss).
        num_matches_by_dim : numpy.ndarray
            Array of shape (d,). The number of matched persistence pairs by dimension.
        num_unmatched_prediction_by_dim : numpy.ndarray
        num_unmatched_target_by_dim : numpy.ndarray
            Arrays of shape (d,). The number of unmatched persistence pairs in the
            prediction and target, respectively, by dimension.
        )")
        .def_readonly("prediction_matches_birth_coordinates", &BettiMatchingResult::predictionMatchesBirthCoordinates)
        .def_readonly("prediction_matches_death_coordinates", &BettiMatchingResult::predictionMatchesDeathCoordinates)
        .def_readonly("target_matches_birth_coordinates", &BettiMatchingResult::targetMatchesBirthCoordinates)
        .def_readonly("target_matches_death_coordinates", &BettiMatchingResult::targetMatchesDeathCoordinates)
        .def_readonly("prediction_unmatched_birth_coordinates", &BettiMatchingResult::predictionUnmatchedBirthCoordinates)
        .def_readonly("prediction_unmatched_death_coordinates", &BettiMatchingResult::predictionUnmatchedDeathCoordinates)
        .def_readonly("target_unmatched_birth_coordinates", &BettiMatchingResult::targetUnmatchedBirthCoordinates)
        .def_readonly("target_unmatched_death_coordinates", &BettiMatchingResult::targetUnmatchedDeathCoordinates)
        .def_readonly("num_matches_by_dim", &BettiMatchingResult::numMatchesByDim)
        .def_readonly("num_unmatched_prediction_by_dim", &BettiMatchingResult::numUnmatchedPredictionByDim)
        .def_readonly("num_unmatched_target_by_dim", &BettiMatchingResult::numUnmatchedTargetByDim)
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
                (self.targetUnmatchedBirthCoordinates.has_value() ? (reprMemberArray("target_unmatched_birth_coordinates", *self.targetUnmatchedBirthCoordinates) + ", ") : "") +
                (self.targetUnmatchedDeathCoordinates.has_value() ? (reprMemberArray("target_unmatched_death_coordinates", *self.targetUnmatchedDeathCoordinates) + ", ") : "") +
                reprMemberArray("num_matches_by_dim", self.numMatchesByDim) + ", " +
                reprMemberArray("num_unmatched_prediction_by_dim", self.numUnmatchedPredictionByDim) +
                (self.numUnmatchedTargetByDim.has_value() ? (", " + reprMemberArray("num_unmatched_target_by_dim", *self.numUnmatchedTargetByDim) + ", ") : "")
            ) + ")";
        });

    py::class_<BarcodeResult>(resultTypesModule, "BarcodeResult",
        R"(
        BarcodeResult

        Holds the result of the barcode computation for a single input volume.

        Each coordinate array contains the birth or death coordinates of persistence
        pairs in all dimensions: first of the 0-dim. pairs, then of the 1-dim. pairs,
        and so on. The starting index of the persistence pairs birth/death coordinates
        in dimension i can be recovered using `np.cumsum()` on `num_matches_by_dim`,
        `num_unmatched_prediction_by_dim`, and `num_unmatched_target_by_dim`,
        respectively.

        Attributes
        ----------
        birth_coordinates : numpy.ndarray
        death_coordinates : numpy.ndarray
            Arrays of shape (n_pairs, d). The birth and death coordinates of
            all persistence pairs in the input volume.
        num_pairs_by_dim : numpy.ndarray
            Array of shape (d,). The number of persistence pairs by dimension.
        )"
        )
        .def_readonly("birth_coordinates", &BarcodeResult::birthCoordinates)
        .def_readonly("death_coordinates", &BarcodeResult::deathCoordinates)
        .def_readonly("num_pairs_by_dim", &BarcodeResult::numPairsByDim)
        .def("__repr__", [](BarcodeResult &self) {
            auto reprMemberArray = [](string name, py::array_t<int64_t>& array) {
                return name + "=" + repr_vector(std::vector<index_t>(array.shape(), array.shape() + array.ndim()), make_tuple("[", "]"));
            };
            return "BarcodeResult(" + (
                reprMemberArray("birth_coordinates", self.birthCoordinates) + ", " +
                reprMemberArray("death_coordinates", self.deathCoordinates) + ", " +
                reprMemberArray("num_pairs_by_dim", self.numPairsByDim)
            ) + ")";
        });

    m.doc() = R"(
        Betti-matching-3D
        =================

        Provides
          1. A fast C++ implementation of the Betti matching algorithm for
            1D, 2D, 3D and n-D images.

        How to compute the Betti matching
        ---------------------------------

        The Betti matching between two NumPy arrays `a`, `b` can be computed using
        the `compute_matching` function:

        ```python
        result = betti_matching.compute_matching(a, b)
        ```

        The result is a `BettiMatchingResult` object that contains the birth and death
        voxel coordinates belonging to matched features and unmatched features in the
        two input arrays, encoded in respective NumPy arrays. 

        Parallel computation of the Betti matching for batches of inputs is supported:

        ```python
        # results[0] constains the matching result between a1 and b1, results[1]
        # contains the matching result between a2 and b2, and so on.
        results = betti_matching.compute_matching([a1, a2, a3], [b1, b2, b3])
        ```

        The Betti matching can also be computed in more fine-grained steps using the
        `BettiMatching` class:

        ```python
        betti_matching = BettiMatching(a, b)
        betti_matching.compute_matching()
        result = betti_matching.get_matching()
        ```
    )";
}
