#include "src_1D/data_structures.h"
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
#include <stdexcept>
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

struct RepresentativeCycleResult {
    optional<vector<py::array_t<int64_t>>> matchedCycles;
    optional<vector<py::array_t<int64_t>>> unmatchedCycles;
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

template<typename... Types>
py::array_t<int64_t> coordinateListToArray(typename vector<std::tuple<Types...>>::iterator coordinates1, typename vector<std::tuple<Types...>>::iterator coordinates2)
{
    size_t size = coordinates2 - coordinates1;
    auto dim = (int)std::tuple_size<std::tuple<Types...>>::value;
    py::array_t<int64_t> coordinateArray({(int)size, dim}, {dim * sizeof(int64_t), sizeof(int64_t)});
    auto coordinateArrayView = coordinateArray.mutable_unchecked();
    for (int i = 0; coordinates1 < coordinates2; i++, coordinates1++) {
        auto& coordinates = *coordinates1;
        std::size_t j = 0;
        std::apply([&](const Types&... args) mutable {
            // Using a fold expression to compute the sum of elements and their indices
            ((coordinateArrayView(i, j++) = args), ...);
        }, coordinates);
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
        input0 : numpy.ndarray
            The first input volume (the "prediction" in the machine
            learning context).
        input1 : numpy.ndarray
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
        inputs0 : list of numpy.ndarray
            The batch of first input volumes (the "predictions" in the machine
            learning context).
        inputs1 : list of numpy.ndarray
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

        Parameters
        ----------
        input0 : numpy.ndarray
            The first input volume (the "prediction" in the machine learning context).
        input1 : numpy.ndarray
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

        .def("compute_representative_cycles",
            [](BettiMatching &self, size_t input, size_t dim, variant<size_t, vector<size_t>, string> &matchedPairsIndices, variant<size_t, vector<size_t>, string> &unmatchedPairsIndices, bool includeDeathVoxel) {
                if (input != 0 && input != 1) {
                    throw invalid_argument("Invalid value for input: " + std::to_string(input));
                }

                // Convert the input arguments to lists of indices (or the empty optional for "all")
                bool matchedNoneRequested = false, unmatchedNoneRequested = false;
                optional<vector<size_t>> matchedPairsIndicesVector, unmatchedPairsIndicesVector;
                if (matchedPairsIndices.index() == 0) {
                    matchedPairsIndicesVector = {std::get<size_t>(matchedPairsIndices)};
                } else if (matchedPairsIndices.index() == 1) {
                    matchedPairsIndicesVector = std::get<vector<size_t>>(matchedPairsIndices);
                } else{
                    auto keyword = std::get<string>(matchedPairsIndices);
                    if (keyword == "all") {} // Nothing to do: matchedPairsIndicesVector is already the empty optional
                    else if (keyword == "none") {
                        matchedNoneRequested = true;
                        matchedPairsIndicesVector = {{}};
                    } else {
                        throw invalid_argument(R"(Invalid value for matched ")" + keyword + R"(": Valid keywords are "all" and "none".)");
                    }
                }
                if (unmatchedPairsIndices.index() == 0) {
                    unmatchedPairsIndicesVector = {std::get<size_t>(unmatchedPairsIndices)};
                } else if (unmatchedPairsIndices.index() == 1) {
                    unmatchedPairsIndicesVector = std::get<vector<size_t>>(unmatchedPairsIndices);
                } else{
                    auto keyword = std::get<string>(unmatchedPairsIndices);
                    if (keyword == "all") {} // Nothing to do: unmatchedPairsIndicesVector is already the empty optional
                    else if (keyword == "none") {
                        unmatchedNoneRequested = true;
                        unmatchedPairsIndicesVector = {{}};
                    } else {
                        throw invalid_argument(R"(Invalid value for unmatched ")" + keyword + R"(": Valid keywords are "all" and "none".)");
                    }
                }
                if (matchedNoneRequested && unmatchedNoneRequested) {
                    throw invalid_argument("No matched or unmatched cycles were requested.");
                }

                // Compute representative cycles
                auto result = self.computeRepresentativeCycles(input, dim, matchedPairsIndicesVector, unmatchedPairsIndicesVector);
                
                // Convert the computed cycles to NumPy arrays
                vector<py::array_t<int64_t>> matchedCyclesArrays;
                vector<py::array_t<int64_t>> unmatchedCyclesArrays;
                switch (self.shape.size()) {
                    case 1: {
                        auto &resultDim1 = std::get<tuple<vector<dim1::RepresentativeCycle>, vector<dim1::RepresentativeCycle>>>(result);
                        for (auto &matchedCycle : std::get<0>(resultDim1)) {
                            matchedCyclesArrays.emplace_back(coordinateListToArray<index_t>(matchedCycle.begin(), matchedCycle.end() - (!includeDeathVoxel)));
                        }
                        for (auto &unmatchedCycle : std::get<1>(resultDim1)) {
                            unmatchedCyclesArrays.emplace_back(coordinateListToArray<index_t>(unmatchedCycle.begin(), unmatchedCycle.end() - (!includeDeathVoxel)));
                        }
                        break;
                    }
                    case 2: {
                        auto &resultDim2 = std::get<tuple<vector<dim2::RepresentativeCycle>, vector<dim2::RepresentativeCycle>>>(result);
                        for (auto &matchedCycle : std::get<0>(resultDim2)) {
                            matchedCyclesArrays.emplace_back(coordinateListToArray<index_t, index_t>(matchedCycle.begin(), matchedCycle.end() - (!includeDeathVoxel)));
                        }
                        for (auto &unmatchedCycle : std::get<1>(resultDim2)) {
                            unmatchedCyclesArrays.emplace_back(coordinateListToArray<index_t, index_t>(unmatchedCycle.begin(), unmatchedCycle.end() - (!includeDeathVoxel)));
                        }
                        break;
                    }
                    case 3: {
                        auto &resultDim3 = std::get<tuple<vector<dim3::RepresentativeCycle>, vector<dim3::RepresentativeCycle>>>(result);
                        for (auto &matchedCycle : std::get<0>(resultDim3)) {
                            matchedCyclesArrays.emplace_back(coordinateListToArray<index_t, index_t, index_t>(matchedCycle.begin(), matchedCycle.end() - (!includeDeathVoxel)));
                        }
                        for (auto &unmatchedCycle : std::get<1>(resultDim3)) {
                            unmatchedCyclesArrays.emplace_back(coordinateListToArray<index_t, index_t, index_t>(unmatchedCycle.begin(), unmatchedCycle.end() - (!includeDeathVoxel)));
                        }
                        break;
                    }
                    default: {
                        throw runtime_error("Invalid value for dim: " + std::to_string(self.shape.size()));
                    }
                }
                return RepresentativeCycleResult{
                    matchedNoneRequested ? std::optional<vector<py::array_t<int64_t>>>{} : matchedCyclesArrays,
                    unmatchedNoneRequested ? std::optional<vector<py::array_t<int64_t>>>{} : unmatchedCyclesArrays
                };
            },
            py::arg("input"),
            py::arg("dim"),
            py::arg("matched") = "none",
            py::arg("unmatched") = "none",
            py::arg("include_death_voxels") = false,
            R"(
            compute_representative_cycles(
                input,
                dim,
                matched="none",
                unmatched="none",
                include_death_voxels=False
            )

            Compute representative cycles for matched and unmatched persistence pairs.
            Lists of matched/unmatched persistence pair indices can be specified for
            which cycles should be computed.
            Computing several cycles in one pass is typically significantly more efficient
            than computing them individually.

            Returned representative cycles are coordinate arrays of shape (len_cycle, d)
            and are voxelized (i.e. coarsened) representations of representatives of
            homology classes belonging to the respective persistence pairs.
            See Notes for more details.

            Parameters
            ----------
            input : int
                The index of the input volume for which to compute the cycles (0 or 1).
            dim : int
                The homology dimension in which to compute the representative cycles
                (for example, for 3D volumes, valid values are 0, 1, and 2).
            matched : int, list of int, "all" or "none", optional
                The indices of the matched pairs for which to compute the representative
                cycles (index i means the i-th matched pair in the specified homology
                dimension). A single index can be passed if only one cycle is requested.
                If "all", all matched cycles in this dimension are computed.
                If "none", no matched cycles are computed (default).
            unmatched : int, list of int, "all" or "none", optional
                The indices of the unmatched pairs for which to compute the representative
                cycles (index i means the i-th unmatched pair in the specified homology
                dimension).  A single index can be passed if only one cycle is requested.
                If "all", all unmatched cycles in this dimension are computed.
                If "none", no unmatched cycles are computed (default).
            include_death_voxels : bool, optional
                Whether to include the death voxel of the persistence pair at the last
                index of the respective representative cycle output array.
                Default is False.

            Returns
            -------
            result : RepresentativeCycleResult
                The representative cycles for the matched and unmatched persistence pairs,
                contained in `result.matched_cycles` and `result.unmatched_cycles` (lists
                of NumPy arrays).

            Format of returned cycles:
            - Each cycle is represented as a NumPy array of shape (len_cycle, d).
            - The entries are voxel coordinates which outline the cycle boundary.
            - If `include_death_voxel == True` is specified, the last coordinate in the
                representative cycle array is the voxel coordinate corresponding to the
                death cube of the persistence pair (a `dim+1`-dimensional cube), which
                can be useful for some visualizations.
            - Cycles for dim=1 in 3D volumes may have duplicate voxel coordinates.

            Notes
            -----
            Choice of representative cycles:
            Representatives cycles of the persistence pairs are not unique, and are chosen
            based on being algorithmically natural. Other algorithms may determine different
            representatives both by choosing a different sublevel set threshold, and by
            choosing a different representative within the homology class belonging to a
            persistence pair at the given threshold.
            - Bottom dimension (dim=0 for 1D, 2D, 3D volumes):
                Representative cycles are connected components just before their death,
                i.e. just before being merged with a component that was born earlier.
            - Inter dimension (dim=1 for 3D volumes):
                Representative cycles are derived from fully reduced working boundaries in
                the matrix reduction algorithm for barcode computation.
                These cycles are the least well-behaved in terms of uniqueness: they may
                be unnecessarily long, and circle unnecessarily many holes.
            - Top dimension (dim=1 for 2D volumes, dim=2 for 3D volumes):
                Representative cycles are determined as boundaries of connected components
                in the dual, at the time of death in the dual (i.e. the birth in the primary
                complex). In other words, they outline loops (2D volumes) or cavities
                (3D volumes) at the time of their creation.
                A voxel counts as being on the boundary if the corresponding vertex is
                a corner vertex of at least one 2-cube (2D) / 3-cube (3D), and at most
                of three 2-cubes (2D) / seven 3-cubes (3D).

            Computation time:
            - For dim=0 in 1D, 2D, 3D volumes / dim=1 in 2D volumes / dim=2 in 3D volumes:
                Typically on the order of computing the barcode in dimension `dim` of
                the input volume, plus the overhead for saving the cycles.
            - For dim=1 in 3D volumes:
                Utilizes the cache of the previous Betti matching computation, hence the
                computation time only comes from saving the cycles to arrays.

            Examples
            --------
            Get the first matched and unmatched dimension 0 representative cycles in the
            prediction volume:
            ```python
            bm = betti_matching.BettiMatching(prediction, target)
            bm.compute_matching()
            cycles = bm.compute_representative_cycles(input=0, dim=0, matched=0, unmatched=0)
            print(cycles.matched_cycles[0])
            print(cycles.unmatched_cycles[0])
            ```

            Get the all matched and all unmatched dimension 0 cycles:
            ```python
            ...
            cycles = bm.compute_representative_cycles(input=0, dim=0,
                matched=np.arange(3), unmatched="all")
            ```

            Create a cycle image showing matched connected components:
            ```python
            bm = betti_matching.BettiMatching(prediction, target)
            bm.compute_matching()
            cycles = bm.compute_representative_cycles(input=0, dim=0, matched="all")
            cycle_image = np.zeros(prediction.shape)
            for i, cycle in reversed(list(enumerate(cycles.matched_cycles))):
                cycle_image[*cycle.T] = i+1
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
    
    py::class_<RepresentativeCycleResult>(resultTypesModule, "RepresentativeCycleResult",
        R"(
        RepresentativeCycleResult

        Holds the result of the computation of representative cycles for matched and
        unmatched persistence pairs in the Betti matching.

        Each cycle is represented as a NumPy array of shape (len_cycle+1, d). The entries
        are voxel coordinates which outline the cycle boundary.

        Attributes
        ----------
        matched_cycles : list of numpy.ndarray, optional
            The representative cycles for the matched persistence pairs. The i-th entry
            is the representative cycle for the i-th matched persistence pair.
        unmatched_cycles : list of numpy.ndarray, optional
            The representative cycles for the unmatched persistence pairs. The i-th entry
            is the representative cycle for the i-th unmatched persistence pair.
        )")
        .def_readonly("matched_cycles", &RepresentativeCycleResult::matchedCycles)
        .def_readonly("unmatched_cycles", &RepresentativeCycleResult::unmatchedCycles)
        .def("__repr__", [](RepresentativeCycleResult &self) {
            return "RepresentativeCycleResult(" + (self.matchedCycles.has_value() ? ("matched_cycles=[...] (" + std::to_string(self.matchedCycles->size()) + " cycles)") : "")
                + (self.matchedCycles.has_value() && self.unmatchedCycles.has_value() ? ", " : "")
                + (self.unmatchedCycles.has_value() ? ("unmatched_cycles=[...] (" + std::to_string(self.unmatchedCycles->size()) + " cycles)") : "") + ")";
        });

    m.doc() = R"(
        Betti-matching-3D (betti_matching)
        ==================================

        Provides
        1. A fast algorithm for computing the Betti matching on 1D, 2D, 3D and
           n-D images, implemented in C++.
        2. A fast algorithm for computing barcodes on cubical grid complexes
           for 1D, 2D, 3D and n-D images - a subprocedure necessary for the
           Betti matching computation, but also exposed separately.
        3. Algorithms for computing representative cycles for matched and unmatched
           persistence pairs in the Betti matching.

        This package uses the V-construction (vertex construction), uses persistence
        modules that come from sublevel sets, and uses a x-y-z-type tiebreaking order
        (or x-y-type/x-type in 2D/1D) in the barcode computation algorithms.


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
        `BettiMatching` class. This allows to cache the computed matching and compute
        representative cycles.

        ```python
        bm = betti_matching.BettiMatching(a, b)
        bm.compute_matching()
        result = bm.get_matching()
        ```
    )";
}
