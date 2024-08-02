#include "BettiMatching.h"
#include "config.h"
#include "data_structures.h"
#include "src_1D/data_structures.h"
#include "utils.h"

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <future>
#include <iostream>
#include <iterator>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>

using namespace std;
namespace py = pybind11;

typedef py::array InputVolume;
typedef py::array_t<value_t, py::array::c_style> TypedInputVolume;
typedef std::tuple<vector<vector<VoxelMatch>>, vector<vector<VoxelPair>>,
                   vector<vector<VoxelPair>>>
    BettiMatchingPairsResult;

struct BettiMatchingResult {
    std::vector<py::array_t<int64_t>> input1MatchedBirthCoordinates;
    std::vector<py::array_t<int64_t>> input1MatchedDeathCoordinates;
    std::vector<py::array_t<int64_t>> input2MatchedBirthCoordinates;
    std::vector<py::array_t<int64_t>> input2MatchedDeathCoordinates;
    py::array_t<int64_t> numMatched;
    // Optional return values if include_input1_unmatched_pairs flag is set
    std::optional<std::vector<py::array_t<int64_t>>>
        input1UnmatchedBirthCoordinates;
    std::optional<std::vector<py::array_t<int64_t>>>
        input1UnmatchedDeathCoordinates;
    std::optional<py::array_t<int64_t>> numUnmatchedInput1;
    // Optional return values if include_input2_unmatched_pairs flag is set (not
    // needed for Betti matching training loss)
    std::optional<std::vector<py::array_t<int64_t>>>
        input2UnmatchedBirthCoordinates;
    std::optional<std::vector<py::array_t<int64_t>>>
        input2UnmatchedDeathCoordinates;
    std::optional<py::array_t<int64_t>> numUnmatchedInput2;

    // Default move constructor, deleted copy constructor
    BettiMatchingResult(BettiMatchingResult &&) noexcept = default;
    BettiMatchingResult(const BettiMatchingResult &) = delete;
    BettiMatchingResult &operator=(const BettiMatchingResult &) = delete;
};
struct BarcodeResult {
    std::vector<py::array_t<int64_t>> birthCoordinates;
    std::vector<py::array_t<int64_t>> deathCoordinates;
    py::array_t<int64_t> numPairs;

    // Default move constructor, deleted copy constructor
    BarcodeResult(BarcodeResult &&) noexcept = default;
    BarcodeResult(const BarcodeResult &) = delete;
    BarcodeResult &operator=(const BarcodeResult &) = delete;
};

struct RepresentativeCycleResult {
    optional<vector<py::array_t<int64_t>>> matchedCycles;
    optional<vector<py::array_t<int64_t>>> unmatchedCycles;
};

template <typename T>
string repr_vector(const vector<T> shape,
                   std::tuple<string, string> parentheses = make_tuple("(",
                                                                       ")"),
                   string separator = ", ") {
    stringstream out_stream;
    out_stream << std::get<0>(parentheses);
    for_each(shape.begin(), std::prev(shape.end()),
             [&out_stream, &separator](auto size) {
                 out_stream << std::to_string(size) << separator;
             });
    out_stream << shape[shape.size() - 1] << std::get<1>(parentheses);
    return out_stream.str();
};

BettiMatchingPairsResult
computeMatchingFromNumpyArrays(TypedInputVolume &input1,
                               TypedInputVolume &input2) {
    vector<index_t> shape1(input1.shape(), input1.shape() + input1.ndim());
    vector<index_t> shape2(input2.shape(), input2.shape() + input2.ndim());
    if (shape1 != shape2) {
        throw invalid_argument(
            "The shapes of the input volumes must agree. Got " +
            repr_vector(shape1) + " and " + repr_vector(shape2));
    }
    vector<value_t> input1Vector(input1.mutable_data(),
                                 input1.mutable_data() + input1.size());
    vector<value_t> input2Vector(input2.mutable_data(),
                                 input2.mutable_data() + input2.size());

    Config config;
    BettiMatching bettiMatching(std::move(input1Vector),
                                std::move(input2Vector), std::move(shape1),
                                std::move(config));
    bettiMatching.computeMatching();
    return bettiMatching.getMatching();
};

vector<vector<VoxelPair>>
computeBarcodeFromNumpyArrays(TypedInputVolume &input1) {
    vector<index_t> shape1(input1.shape(), input1.shape() + input1.ndim());
    vector<value_t> input1Vector(input1.mutable_data(),
                                 input1.mutable_data() + input1.size());
    vector<value_t> input1VectorCopy(input1.mutable_data(),
                                     input1.mutable_data() + input1.size());

    Config config;
    BettiMatching bettiMatching(std::move(input1Vector),
                                std::move(input1VectorCopy), std::move(shape1),
                                std::move(config));
    return bettiMatching.computePairsInput0();
};

BettiMatchingResult convertPairResultToArrayResult(
    BettiMatchingPairsResult &result, size_t numDimensions,
    bool includeInput1UnmatchedPairs, bool includeInput2UnmatchedPairs);
BarcodeResult
convertPairResultToBarcodeResult(vector<vector<VoxelPair>> &pairsByDimension);

vector<BettiMatchingResult> computeMatchingFromInputs(
    vector<InputVolume> &untypedInputs1, vector<InputVolume> &untypedInputs2,
    bool includeInput1UnmatchedPairs, bool includeInput2UnmatchedPairs) {
    // Validate and convert inputs
    vector<TypedInputVolume> inputs1;
    vector<TypedInputVolume> inputs2;
    if (untypedInputs1.size() != untypedInputs2.size()) {
        throw invalid_argument(
            "Different numbers of inputs where provided: " +
            std::to_string(untypedInputs1.size()) + " (inputs1) and " +
            std::to_string(untypedInputs2.size()) + " (inputs2)");
    }
    for (size_t i = 0; i < untypedInputs1.size(); i++) {
        vector<index_t> shape1(untypedInputs1[i].shape(),
                               untypedInputs1[i].shape() +
                                   untypedInputs1[i].ndim());
        vector<index_t> shape2(untypedInputs2[i].shape(),
                               untypedInputs2[i].shape() +
                                   untypedInputs2[i].ndim());
        if (shape1 != shape2) {
            throw invalid_argument(
                "The shapes of the input volumes must agree. Got " +
                repr_vector(shape1) + " and " + repr_vector(shape2));
        }
        inputs1.push_back(untypedInputs1[i].cast<TypedInputVolume>());
        inputs2.push_back(untypedInputs2[i].cast<TypedInputVolume>());
    }
    size_t numDimensions = inputs1[0].ndim();

    // Create one asynchronous task for Betti matching computation per input in
    // the batch
    std::vector<std::future<BettiMatchingPairsResult>> resultFutures;
    for (int i = 0; i < inputs1.size(); i++) {
        resultFutures.push_back(std::async(computeMatchingFromNumpyArrays,
                                           std::ref(inputs1[i]),
                                           std::ref(inputs2[i])));
    }

    // Now block on all of them one at a time.
    std::vector<BettiMatchingPairsResult> pairResults;
    for (auto &future : resultFutures) {
        pairResults.push_back(future.get());
    }

    vector<BettiMatchingResult> arrayResults;

    // Write the results into numpy arrays
    for (auto &result : pairResults) {
        auto arrayResult = convertPairResultToArrayResult(
            result, numDimensions, includeInput1UnmatchedPairs,
            includeInput2UnmatchedPairs);
        arrayResults.emplace_back(std::move(arrayResult));
    }
    return arrayResults;
}

vector<BarcodeResult>
computeBarcodeFromInputs(vector<InputVolume> &untypedInputs) {
    vector<TypedInputVolume> inputs;
    for (auto &untypedInput : untypedInputs) {
        inputs.push_back(untypedInput.cast<TypedInputVolume>());
    }

    std::vector<vector<vector<VoxelPair>>> pairResults;
    // Create one asynchronous task for barcode computation per input in the
    // batch
    std::vector<std::future<vector<vector<VoxelPair>>>> resultFutures;
    for (int i = 0; i < inputs.size(); i++) {
        resultFutures.push_back(
            std::async(computeBarcodeFromNumpyArrays, std::ref(inputs[i])));
    }

    // Now block on all of them one at a time.
    for (auto &future : resultFutures) {
        pairResults.push_back(future.get());
    }

    // Write the results into numpy arrays
    vector<BarcodeResult> arrayResults;
    for (auto &result : pairResults) {
        arrayResults.emplace_back(convertPairResultToBarcodeResult(result));
    }
    return arrayResults;
}

BettiMatchingResult convertPairResultToArrayResult(
    BettiMatchingPairsResult &result, size_t numDimensions,
    bool includeInput1UnmatchedPairs, bool includeInput2UnmatchedPairs) {
    auto &matchesByDimension = std::get<0>(result);
    auto &input1UnmatchedByDimension = std::get<1>(result);
    auto &input2UnmatchedByDimension = std::get<2>(result);

    vector<py::array_t<int64_t>> input1MatchedBirthCoordinates(numDimensions),
        input1MatchedDeathCoordinates(numDimensions),
        input2MatchedBirthCoordinates(numDimensions),
        input2MatchedDeathCoordinates(numDimensions),
        input1UnmatchedBirthCoordinates(numDimensions),
        input1UnmatchedDeathCoordinates(numDimensions),
        input2UnmatchedBirthCoordinates(numDimensions),
        input2UnmatchedDeathCoordinates(numDimensions);
    // Shape: (num dimensions,)
    py::array_t<int64_t> numMatchedByDim({numDimensions}, {sizeof(int64_t)});
    auto numMatchedByDimView = numMatchedByDim.mutable_unchecked();
    py::array_t<int64_t> numUnmatchedInput1ByDim;
    if (includeInput1UnmatchedPairs) {
        numUnmatchedInput1ByDim =
            py::array_t<int64_t>({numDimensions}, {sizeof(int64_t)});
    }
    py::array_t<int64_t> numUnmatchedInput2ByDim;
    if (includeInput2UnmatchedPairs) {
        numUnmatchedInput2ByDim =
            py::array_t<int64_t>({numDimensions}, {sizeof(int64_t)});
    }
    for (size_t dimension = 0; dimension < numDimensions; dimension++) {
        // Shape: (num matched in dimension, num dimensions)
        input1MatchedBirthCoordinates[dimension] = py::array_t<int64_t>(
            {matchesByDimension[dimension].size(), numDimensions},
            {numDimensions * sizeof(int64_t), sizeof(int64_t)});
        input1MatchedDeathCoordinates[dimension] = py::array_t<int64_t>(
            {matchesByDimension[dimension].size(), numDimensions},
            {numDimensions * sizeof(int64_t), sizeof(int64_t)});
        input2MatchedBirthCoordinates[dimension] = py::array_t<int64_t>(
            {matchesByDimension[dimension].size(), numDimensions},
            {numDimensions * sizeof(int64_t), sizeof(int64_t)});
        input2MatchedDeathCoordinates[dimension] = py::array_t<int64_t>(
            {matchesByDimension[dimension].size(), numDimensions},
            {numDimensions * sizeof(int64_t), sizeof(int64_t)});
        if (includeInput1UnmatchedPairs) {
            // Shape: (num unmatched (input1) in dimension, num dimensions)
            input1UnmatchedBirthCoordinates[dimension] = py::array_t<int64_t>(
                {input1UnmatchedByDimension[dimension].size(), numDimensions},
                {numDimensions * sizeof(int64_t), sizeof(int64_t)});
            input1UnmatchedDeathCoordinates[dimension] = py::array_t<int64_t>(
                {input1UnmatchedByDimension[dimension].size(), numDimensions},
                {numDimensions * sizeof(int64_t), sizeof(int64_t)});
        }
        if (includeInput2UnmatchedPairs) {
            // Shape: (num unmatched (inputd) in dimension, num dimensions)
            input2UnmatchedBirthCoordinates[dimension] = py::array_t<int64_t>(
                {input2UnmatchedByDimension[dimension].size(), numDimensions},
                {numDimensions * sizeof(int64_t), sizeof(int64_t)});
            input2UnmatchedDeathCoordinates[dimension] = py::array_t<int64_t>(
                {input2UnmatchedByDimension[dimension].size(), numDimensions},
                {numDimensions * sizeof(int64_t), sizeof(int64_t)});
        }

        auto input1MatchedBirthCoordinatesView =
            input1MatchedBirthCoordinates[dimension].mutable_unchecked();
        auto input1MatchedDeathCoordinatesView =
            input1MatchedDeathCoordinates[dimension].mutable_unchecked();
        auto input2MatchedBirthCoordinatesView =
            input2MatchedBirthCoordinates[dimension].mutable_unchecked();
        auto input2MatchedDeathCoordinatesView =
            input2MatchedDeathCoordinates[dimension].mutable_unchecked();

        int i = 0;
        for (auto &match : matchesByDimension[dimension]) {
            for (size_t d = 0; d < numDimensions; d++) {
                input1MatchedBirthCoordinatesView(i, d) = match.pair0.birth[d];
                input1MatchedDeathCoordinatesView(i, d) = match.pair0.death[d];
                input2MatchedBirthCoordinatesView(i, d) = match.pair1.birth[d];
                input2MatchedDeathCoordinatesView(i, d) = match.pair1.death[d];
            }
            i++;
        }
        numMatchedByDimView(dimension) = matchesByDimension[dimension].size();

        if (includeInput1UnmatchedPairs) {
            i = 0;
            auto input1UnmatchedBirthCoordinatesView =
                input1UnmatchedBirthCoordinates[dimension].mutable_unchecked();
            auto input1UnmatchedDeathCoordinatesView =
                input1UnmatchedDeathCoordinates[dimension].mutable_unchecked();
            auto numUnmatchedInput1ByDimView =
                numUnmatchedInput1ByDim.mutable_unchecked();
            for (auto &unmatched : input1UnmatchedByDimension[dimension]) {
                for (size_t d = 0; d < numDimensions; d++) {
                    input1UnmatchedBirthCoordinatesView(i, d) =
                        unmatched.birth[d];
                    input1UnmatchedDeathCoordinatesView(i, d) =
                        unmatched.death[d];
                }
                i++;
            }
            numUnmatchedInput1ByDimView(dimension) =
                input1UnmatchedByDimension[dimension].size();
        }
        if (includeInput2UnmatchedPairs) {
            i = 0;
            auto input2UnmatchedBirthCoordinatesView =
                input2UnmatchedBirthCoordinates[dimension].mutable_unchecked();
            auto input2UnmatchedDeathCoordinatesView =
                input2UnmatchedDeathCoordinates[dimension].mutable_unchecked();
            auto numUnmatchedInput2ByDimView =
                numUnmatchedInput2ByDim.mutable_unchecked();
            for (auto &unmatched : input2UnmatchedByDimension[dimension]) {
                for (int d = 0; d < numDimensions; d++) {
                    input2UnmatchedBirthCoordinatesView(i, d) =
                        unmatched.birth[d];
                    input2UnmatchedDeathCoordinatesView(i, d) =
                        unmatched.death[d];
                }
                i++;
            }
            numUnmatchedInput2ByDimView(dimension) =
                input2UnmatchedByDimension[dimension].size();
        }
    }
    auto arrayResult = BettiMatchingResult{
        std::move(input1MatchedBirthCoordinates),
        std::move(input1MatchedDeathCoordinates),
        std::move(input2MatchedBirthCoordinates),
        std::move(input2MatchedDeathCoordinates), std::move(numMatchedByDim)};
    if (includeInput1UnmatchedPairs) {
        arrayResult.input1UnmatchedBirthCoordinates =
            std::move(input1UnmatchedBirthCoordinates);
        arrayResult.input1UnmatchedDeathCoordinates =
            std::move(input1UnmatchedDeathCoordinates);
        arrayResult.numUnmatchedInput1 = std::move(numUnmatchedInput1ByDim);
    }
    if (includeInput2UnmatchedPairs) {
        arrayResult.input2UnmatchedBirthCoordinates =
            std::move(input2UnmatchedBirthCoordinates);
        arrayResult.input2UnmatchedDeathCoordinates =
            std::move(input2UnmatchedDeathCoordinates);
        arrayResult.numUnmatchedInput2 = std::move(numUnmatchedInput2ByDim);
    }
    return arrayResult;
}

BarcodeResult
convertPairResultToBarcodeResult(vector<vector<VoxelPair>> &pairsByDimension) {
    auto numDimensions = pairsByDimension.size();

    vector<py::array_t<int64_t>> birthCoordinates(numDimensions),
        deathCoordinates(numDimensions);
    // Shape: (num dimensions,)
    py::array_t<int64_t> numPairsByDim({numDimensions}, {sizeof(int64_t)});
    auto numPairsByDimView = numPairsByDim.mutable_unchecked();

    for (size_t dimension = 0; dimension < numDimensions; dimension++) {
        // Shape: (num pairs in dimension, num dimensions)
        birthCoordinates[dimension] = py::array_t<int64_t>(
            {pairsByDimension[dimension].size(), numDimensions},
            {numDimensions * sizeof(int64_t), sizeof(int64_t)});
        deathCoordinates[dimension] = py::array_t<int64_t>(
            {pairsByDimension[dimension].size(), numDimensions},
            {numDimensions * sizeof(int64_t), sizeof(int64_t)});

        auto birthCoordinatesView =
            birthCoordinates[dimension].mutable_unchecked();
        auto deathCoordinatesView =
            deathCoordinates[dimension].mutable_unchecked();

        int i = 0;
        for (auto &pair : pairsByDimension[dimension]) {
            for (int d = 0; d < numDimensions; d++) {
                birthCoordinatesView(i, d) = pair.birth[d];
                deathCoordinatesView(i, d) = pair.death[d];
            }
            i++;
        }
        numPairsByDimView(dimension) = pairsByDimension[dimension].size();
    }

    auto arrayResult = BarcodeResult{birthCoordinates = birthCoordinates,
                                     deathCoordinates = deathCoordinates,
                                     numPairsByDim = numPairsByDim};
    return arrayResult;
}

template <typename... Types>
py::array_t<int64_t> coordinateListToArray(
    typename vector<std::tuple<Types...>>::iterator coordinates1,
    typename vector<std::tuple<Types...>>::iterator coordinates2) {
    size_t size = coordinates2 - coordinates1;
    auto dim = (int)std::tuple_size<std::tuple<Types...>>::value;
    py::array_t<int64_t> coordinateArray(
        {(int)size, dim}, {dim * sizeof(int64_t), sizeof(int64_t)});
    auto coordinateArrayView = coordinateArray.mutable_unchecked();
    for (int i = 0; coordinates1 < coordinates2; i++, coordinates1++) {
        auto &coordinates = *coordinates1;
        std::size_t j = 0;
        std::apply(
            [&](const Types &...args) mutable {
                // Using a fold expression to compute the sum of elements
                // and their indices
                ((coordinateArrayView(i, j++) = args), ...);
            },
            coordinates);
    }
    return coordinateArray;
}

PYBIND11_MODULE(betti_matching, m) {
    m.def(
        "compute_matching",
        [](InputVolume &untypedInput1, InputVolume &untypedInput2,
           bool includeInput1UnmatchedPairs, bool includeInput2UnmatchedPairs) {
            vector<InputVolume> untypedInputs1 = {untypedInput1};
            vector<InputVolume> untypedInputs2 = {untypedInput2};
            return std::move(computeMatchingFromInputs(
                untypedInputs1, untypedInputs2, includeInput1UnmatchedPairs,
                includeInput2UnmatchedPairs)[0]);
        },
        py::arg("input1").noconvert(), py::arg("input2").noconvert(),
        py::arg("include_input1_unmatched_pairs") = true,
        py::arg("include_input2_unmatched_pairs") = true,
        R"(
        compute_matching(
            input1,
            input2,
            include_input1_unmatched_pairs=True,
            include_input2_unmatched_pairs=True
        )

        Compute the Betti matching between two input volumes.

        Parameters
        ----------
        input1 : numpy.ndarray
            The first input volume (the "prediction" in the machine
            learning context).
        input2 : numpy.ndarray
            The second input volume (the "target" in the machine
            learning context).
        include_input1_unmatched_pairs : bool, optional
            Whether to include the unmatched pairs in the input1 volume
            in the result. Default is True.
        include_input2_unmatched_pairs : bool, optional
            Whether to include the unmatched pairs in the input2 volume
            in the result. Default is True. Deactivate when the target
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
            result.input1_matched_birth_coordinates[0])
        a_dim0_matches_death_coordinates = (
            result.input1_matched_death_coordinates[0])
        a_dim0_matched_bars_lengths = (a[*a_dim0_matches_death_coordinates.T]
            - a[*a_dim0_matches_birth_coordinates.T])
        ```
        )");

    m.def("compute_matching", &computeMatchingFromInputs,
          py::arg("inputs1").noconvert(), py::arg("inputs2").noconvert(),
          py::arg("include_input1_unmatched_pairs") = true,
          py::arg("include_input2_unmatched_pairs") = true,
          R"(
        compute_matching(
            inputs1,
            inputs2,
            include_input1_unmatched_pairs=True,
            include_input2_unmatched_pairs=True
        )

        Compute the Betti matching between two batches of input volumes in
        parallel. The Betti matching computation are parallelized using
        std::async.

        Parameters
        ----------
        inputs1 : list of numpy.ndarray
            The batch of first input volumes (the "predictions" in the machine
            learning context).
        inputs2 : list of numpy.ndarray
            The batch of second input volumes (the "targets" in the machine
            learning context).
        include_input1_unmatched_pairs : bool, optional
            Whether to include the unmatched pairs in the input1 volumes
            in the result. Default is True.
        include_input2_unmatched_pairs : bool, optional
            Whether to include the unmatched pairs in the input2 volumes
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
        num_matched_a1_b1 = results[0].num_matched.sum()
        num_matched_a2_b2 = results[1].num_matched.sum()
        ```
        )");

    m.def(
        "compute_barcode",
        [](InputVolume &untypedInput) {
            vector<InputVolume> untypedInputs = {untypedInput};
            return std::move(computeBarcodeFromInputs(untypedInputs)[0]);
        },
        py::arg("input").noconvert(),
        R"(
        compute_barcode(input)

        Compute the barcode for a single input.

        Parameters
        ----------
        input: numpy.ndarray
            The input volume.

        Returns
        -------
        result: BarcodeResult
            The result of the barcode computation. Contains the birth and death
            coordinate for each persistence pair, from which the birth and death
            intensity values can be recovered.
            See the documentation of `betti_matching.return_types.BarcodeResult`
            for details about the contained data.

        Example
        -------
        ```
        a = np.random.rand(10, 10, 10)
        result = betti_matching.compute_barcode(a)
        births_deaths_dim1 = np.stack([
            a[*result.birth_coordinates[1].T],
            a[*result.death_coordinates[1].T]
        ], axis=0)
        ```
        )");

    m.def("compute_barcode", &computeBarcodeFromInputs,
          py::arg("inputs").noconvert(),
          R"(
            compute_barcode(inputs)

            Compute the barcode for a batch of inputs. The independent computations
            are parallelized using std::async.

            Parameters
            ----------
            inputs: list of numpy.ndarray
                The input volumes.

            Returns
            -------
            results: list of BarcodeResult
                The results of the barcode computations. Each result contains the birth
                and death coordinate for each persistence pair, from which the birth
                and death intensity values can be recovered.
                See the documentation of `betti_matching.return_types.BarcodeResult`
                for details about the contained data.

            Example
            -------
            ```
            a1, a2 = np.random.rand(10, 10, 10), np.random.rand(8, 8, 8)
            results = betti_matching.compute_barcode([a1, a2])
            births_deaths_dim1 = np.stack([
                a[*results[0].birth_coordinates[1].T],
                a[*results[0].death_coordinates[1].T]
            ], axis=0)
            ```
            )");

    py::class_<BettiMatching>(m, "BettiMatching",
                              R"(
        BettiMatching(input1, input2)

        Class for computing the Betti matching between two input volumes, and computing
        representative cycles. In comparison to `betti_matching.compute_matching()`,
        this class allows to cache the computed matching and use it to compute
        representative cycles for the matched and unmatched features.

        Parameters
        ----------
        input1 : numpy.ndarray
            The first input volume (the "prediction" in the machine learning context).
        input2 : numpy.ndarray
            The second input volume (the "target" in the machine learning context).

        Example
        -------
        ```python
        a, b = np.random.rand(10, 10, 10), np.random.rand(10, 10, 10)
        bm = betti_matching.BettiMatching(a, b)
        bm.compute_matching()
        result = bm.get_matching()
        ```
        )")
        .def(py::init([](InputVolume &untypedInput1,
                         InputVolume &untypedInput2) {
                 auto input1 = untypedInput1.cast<TypedInputVolume>();
                 auto input2 = untypedInput2.cast<TypedInputVolume>();
                 vector<index_t> shape1(input1.shape(),
                                        input1.shape() + input1.ndim());
                 vector<index_t> shape2(input2.shape(),
                                        input2.shape() + input2.ndim());
                 if (shape1 != shape2) {
                     throw invalid_argument(
                         "The shapes of the input volumes must agree. Got " +
                         repr_vector(shape1) + " and " + repr_vector(shape2));
                 }
                 vector<value_t> input1Vector(input1.mutable_data(),
                                              input1.mutable_data() +
                                                  input1.size());
                 vector<value_t> input2Vector(input2.mutable_data(),
                                              input2.mutable_data() +
                                                  input2.size());
                 Config config;
                 return BettiMatching(std::move(input1Vector),
                                      std::move(input2Vector),
                                      std::move(shape1), std::move(config));
             }),
             py::arg("input1").noconvert(), py::arg("input2").noconvert())

        .def_readwrite("shape", &BettiMatching::shape)

        .def("compute_matching", &BettiMatching::computeMatching,
             R"(
            compute_matching()

            Compute the Betti matching between the two input volumes.
            Skips the computation if it has already been performed previously on
            the same `BettiMatching` instance.
            )")

        .def(
            "get_matching",
            [](BettiMatching &self, bool includeInput1UnmatchedPairs,
               bool includeInput2UnmatchedPairs) {
                BettiMatchingPairsResult result = self.getMatching();
                return convertPairResultToArrayResult(
                    result, self.shape.size(), includeInput1UnmatchedPairs,
                    includeInput2UnmatchedPairs);
            },
            py::arg("include_input1_unmatched_pairs") = true,
            py::arg("include_input2_unmatched_pairs") = true,
            R"(
            get_matching()

            Retrieve the computed matching. `compute_matching()` must be called
            before calling this method.

            Parameters
            ----------
            include_input1_unmatched_pairs : bool, optional
                Whether to include the unmatched pairs in the input1 volume
                in the result. Default is True.
            include_input2_unmatched_pairs : bool, optional
                Whether to include the unmatched pairs in the input2 volume
                in the result. Default is True. Can be deactivated when the target
                unmatched pairs are not needed, such as in training with the Betti
                matching loss, where they do not contribute to the gradient.

            Returns
            -------
            result : BettiMatchingResult
                The result of the Betti matching between the two input volumes.
                See the documentation of `betti_matching.return_types.BettiMatchingResult`
                for details about the contained data.
            )")

        .def("print", &BettiMatching::printResult,
             R"(
            print()

            Print a summary of the computed matching. `compute_matching()` must be called
            before calling this method.
            )")

        .def(
            "compute_representative_cycles",
            [](BettiMatching &self, variant<size_t, string> input, size_t dim,
               variant<size_t, vector<size_t>, string> &matchedPairsIndices,
               variant<size_t, vector<size_t>, string> &unmatchedPairsIndices,
               bool includeDeathVoxel, bool deduplicateVoxels) {
                int inputNumber;
                if (input.index() == 0) {
                    inputNumber = std::get<size_t>(input) - 1;
                    if (inputNumber != 0 && inputNumber != 1) {
                        throw invalid_argument(
                            "Invalid value for input: " +
                            std::to_string(std::get<size_t>(input)));
                    }
                } else {
                    string inputKeyword = std::get<string>(input);
                    if (inputKeyword == "input1") {
                        inputNumber = 0;
                    } else if (inputKeyword == "input2") {
                        inputNumber = 1;
                    } else if (inputKeyword == "comparison") {
                        inputNumber = 2;
                    } else {
                        throw invalid_argument(
                            R"(Invalid value for input ")" + inputKeyword +
                            R"(": Valid keywords are "input1", "input2", and "comparison".)");
                    }
                }

                // Convert the input arguments to lists of indices (or the
                // empty optional for "all")
                bool matchedNoneRequested = false,
                     unmatchedNoneRequested = false;
                optional<vector<size_t>> matchedPairsIndicesVector,
                    unmatchedPairsIndicesVector;
                if (matchedPairsIndices.index() == 0) {
                    matchedPairsIndicesVector = {
                        std::get<size_t>(matchedPairsIndices)};
                } else if (matchedPairsIndices.index() == 1) {
                    matchedPairsIndicesVector =
                        std::get<vector<size_t>>(matchedPairsIndices);
                } else {
                    auto keyword = std::get<string>(matchedPairsIndices);
                    if (keyword == "all") {
                    } // Nothing to do: matchedPairsIndicesVector is already
                      // the empty optional
                    else if (keyword == "none") {
                        matchedNoneRequested = true;
                        matchedPairsIndicesVector = {vector<size_t>()};
                    } else {
                        throw invalid_argument(
                            R"(Invalid value for matched ")" + keyword +
                            R"(": Valid keywords are "all" and "none".)");
                    }
                }
                if (unmatchedPairsIndices.index() == 0) {
                    unmatchedPairsIndicesVector = {
                        std::get<size_t>(unmatchedPairsIndices)};
                } else if (unmatchedPairsIndices.index() == 1) {
                    unmatchedPairsIndicesVector =
                        std::get<vector<size_t>>(unmatchedPairsIndices);
                } else {
                    auto keyword = std::get<string>(unmatchedPairsIndices);
                    if (keyword == "all") {
                    } // Nothing to do: unmatchedPairsIndicesVector is
                      // already the empty optional
                    else if (keyword == "none") {
                        unmatchedNoneRequested = true;
                        unmatchedPairsIndicesVector = {vector<size_t>()};
                    } else {
                        throw invalid_argument(
                            R"(Invalid value for unmatched ")" + keyword +
                            R"(": Valid keywords are "all" and "none".)");
                    }
                }
                if (matchedNoneRequested && unmatchedNoneRequested) {
                    throw invalid_argument(
                        "No matched or unmatched cycles were requested.");
                }

                if (inputNumber == 2 && !unmatchedNoneRequested) {
                    throw invalid_argument(
                        "Computing representative cycles for unmatched "
                        "pairs "
                        "in the comparison image is not supported.");
                }

                // Compute representative cycles
                auto result = self.computeRepresentativeCycles(
                    inputNumber, dim, matchedPairsIndicesVector,
                    unmatchedPairsIndicesVector);

                // Convert the computed cycles to NumPy arrays
                vector<py::array_t<int64_t>> matchedCyclesArrays;
                vector<py::array_t<int64_t>> unmatchedCyclesArrays;
                switch (self.shape.size()) {
                case 1: {
                    auto &resultDim1 =
                        std::get<tuple<vector<dim1::RepresentativeCycle>,
                                       vector<dim1::RepresentativeCycle>>>(
                            result);
                    for (auto &matchedCycle : std::get<0>(resultDim1)) {
                        matchedCyclesArrays.emplace_back(
                            coordinateListToArray<index_t>(
                                matchedCycle.begin(),
                                matchedCycle.end() - (!includeDeathVoxel)));
                    }
                    for (auto &unmatchedCycle : std::get<1>(resultDim1)) {
                        unmatchedCyclesArrays.emplace_back(
                            coordinateListToArray<index_t>(
                                unmatchedCycle.begin(),
                                unmatchedCycle.end() - (!includeDeathVoxel)));
                    }
                    break;
                }
                case 2: {
                    auto &resultDim2 =
                        std::get<tuple<vector<dim2::RepresentativeCycle>,
                                       vector<dim2::RepresentativeCycle>>>(
                            result);
                    for (auto &matchedCycle : std::get<0>(resultDim2)) {
                        matchedCyclesArrays.emplace_back(
                            coordinateListToArray<index_t, index_t>(
                                matchedCycle.begin(),
                                matchedCycle.end() - (!includeDeathVoxel)));
                    }
                    for (auto &unmatchedCycle : std::get<1>(resultDim2)) {
                        unmatchedCyclesArrays.emplace_back(
                            coordinateListToArray<index_t, index_t>(
                                unmatchedCycle.begin(),
                                unmatchedCycle.end() - (!includeDeathVoxel)));
                    }
                    break;
                }
                case 3: {
                    auto &resultDim3 =
                        std::get<tuple<vector<dim3::RepresentativeCycle>,
                                       vector<dim3::RepresentativeCycle>>>(
                            result);
                    for (auto &matchedCycle : std::get<0>(resultDim3)) {
                        if (dim == 1 && deduplicateVoxels) {
                            // Sort everything (but the death voxel in the
                            // last position!), then remove duplicates
                            std::sort(matchedCycle.begin(),
                                      matchedCycle.end() - 1);
                            matchedCycle.erase(
                                std::unique(matchedCycle.begin(),
                                            matchedCycle.end() - 1),
                                matchedCycle.end() - 1);
                        }
                        matchedCyclesArrays.emplace_back(
                            coordinateListToArray<index_t, index_t, index_t>(
                                matchedCycle.begin(),
                                matchedCycle.end() - (!includeDeathVoxel)));
                    }
                    for (auto &unmatchedCycle : std::get<1>(resultDim3)) {
                        if (dim == 1 && deduplicateVoxels) {
                            // Sort everything (but the death voxel in the
                            // last position!), then remove duplicates
                            std::sort(unmatchedCycle.begin(),
                                      unmatchedCycle.end() - 1);
                            unmatchedCycle.erase(
                                std::unique(unmatchedCycle.begin(),
                                            unmatchedCycle.end() - 1),
                                unmatchedCycle.end() - 1);
                        }
                        unmatchedCyclesArrays.emplace_back(
                            coordinateListToArray<index_t, index_t, index_t>(
                                unmatchedCycle.begin(),
                                unmatchedCycle.end() - (!includeDeathVoxel)));
                    }
                    break;
                }
                default: {
                    throw runtime_error("Invalid value for dim: " +
                                        std::to_string(self.shape.size()));
                }
                }
                return RepresentativeCycleResult{
                    matchedNoneRequested
                        ? std::optional<vector<py::array_t<int64_t>>>{}
                        : matchedCyclesArrays,
                    unmatchedNoneRequested
                        ? std::optional<vector<py::array_t<int64_t>>>{}
                        : unmatchedCyclesArrays};
            },
            py::arg("input"), py::arg("dim"), py::arg("matched") = "none",
            py::arg("unmatched") = "none",
            py::arg("include_death_voxels") = false,
            py::arg("deduplicate_voxels") = false,
            R"(
            compute_representative_cycles(
                input,
                dim,
                matched="none",
                unmatched="none",
                include_death_voxels=False,
                deduplicate_voxels=False
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
            input : one of {0, 1}, or one of {"input1", "input2", "comparison"}
                The index of the input volume for which to compute the cycles.
                - First input volume: 0 or "input1".
                - Second input volume: 1 or "input2".
                - Comparison image: "comparison" (supports only matched cycles).
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
                Indices must be strictly ascending.
            include_death_voxels : bool, optional
                Whether to include the death voxel of the persistence pair at the last
                index of the respective representative cycle output array.
                Default is False.
            deduplicate_voxels : bool, optional
                Whether to deduplicate voxel coordinates in the representative cycles.
                Only relevant for dim=1 in 3D volumes, where cycles may have duplicate
                voxel coordinates.
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
            )");

    py::class_<dim3::Cube>(m, "Cube")
        .def("x", &dim3::Cube::x)
        .def("y", &dim3::Cube::y)
        .def("z", &dim3::Cube::z)
        .def("type", &dim3::Cube::type)
        .def("__repr__", [](dim3::Cube &self) {
            return "dim3::Cube(x=" + std::to_string(self.x()) +
                   ", y=" + std::to_string(self.y()) +
                   ", z=" + std::to_string(self.z()) + ")";
        });

    auto resultTypesModule =
        m.def_submodule("return_types", "Return types for betti_matching");

    py::class_<BettiMatchingResult>(resultTypesModule, "BettiMatchingResult",
                                    R"(
        BettiMatchingResult

        Holds the result of the Betti matching between two arrays.
        The result contains the coordinates of the birth and death voxels for each
        matched and unmatched persistence pair, represented as NumPy arrays.

        Each coordinate array contains the birth or death coordinates of persistence
        pairs in all dimensions: first of the 0-dim. pairs, then of the 1-dim. pairs,
        and so on. The starting index of the persistence pairs birth/death coordinates
        in dimension i can be recovered using `np.cumsum()` on `num_matched`,
        `num_unmatched_input1`, and `num_unmatched_input2`,
        respectively.

        Attributes
        ----------
        input1_matched_birth_coordinates : numpy.ndarray
        input1_matched_death_coordinates : numpy.ndarray
        input2_matched_birth_coordinates : numpy.ndarray
        input2_matched_death_coordinates : numpy.ndarray
            Arrays of shape (n_matched, d). The birth and death coordinates of matched
            persistence pairs in the first and second input, respectively.
        input1_unmatched_birth_coordinates : numpy.ndarray
        input1_unmatched_death_coordinates : numpy.ndarray
            Arrays of shape (n_unmatched_input1, d). The birth and death coordinates
            of unmatched persistence pairs in the first input.
        input2_unmatched_birth_coordinates : numpy.ndarray, optional
        input2_unmatched_death_coordinates : numpy.ndarray, optional
            Arrays of shape (n_unmatched_input2, d). The birth and death
            coordinates of unmatched persistence pairs in the second input. Only present
            if the `include_input2_unmatched_pairs` flag was set to `True` in the
            Betti matching computation (can be turned off since the target unmatched
            pairs do not contribute to the gradient when training with the Betti
            matching loss).
        num_matched : numpy.ndarray
            Array of shape (d,). The number of matched persistence pairs by dimension.
        num_unmatched_input1 : numpy.ndarray
        num_unmatched_input2 : numpy.ndarray
            Arrays of shape (d,). The number of unmatched persistence pairs in the
            first and second input, respectively, by dimension.
        )")
        .def_readonly("input1_matched_birth_coordinates",
                      &BettiMatchingResult::input1MatchedBirthCoordinates)
        .def_readonly("input1_matched_death_coordinates",
                      &BettiMatchingResult::input1MatchedDeathCoordinates)
        .def_readonly("input2_matched_birth_coordinates",
                      &BettiMatchingResult::input2MatchedBirthCoordinates)
        .def_readonly("input2_matched_death_coordinates",
                      &BettiMatchingResult::input2MatchedDeathCoordinates)
        .def_readonly("input1_unmatched_birth_coordinates",
                      &BettiMatchingResult::input1UnmatchedBirthCoordinates)
        .def_readonly("input1_unmatched_death_coordinates",
                      &BettiMatchingResult::input1UnmatchedDeathCoordinates)
        .def_readonly("input2_unmatched_birth_coordinates",
                      &BettiMatchingResult::input2UnmatchedBirthCoordinates)
        .def_readonly("input2_unmatched_death_coordinates",
                      &BettiMatchingResult::input2UnmatchedDeathCoordinates)
        .def_readonly("num_matched", &BettiMatchingResult::numMatched)
        .def_readonly("num_unmatched_input1",
                      &BettiMatchingResult::numUnmatchedInput1)
        .def_readonly("num_unmatched_input2",
                      &BettiMatchingResult::numUnmatchedInput2)
        .def("__repr__", [](BettiMatchingResult &self) {
            auto reprSmall1DimMemberArray = [](string name,
                                               py::array_t<int64_t> &array) {
                vector<int64_t> values(array.size());
                for (int i = 0; i < array.size(); i++) {
                    values[i] = array.at(i);
                }
                return name + "=" + repr_vector(values, make_tuple("[", "]"));
            };
            auto reprMemberVector = [](string name,
                                       vector<py::array_t<int64_t>> &vector) {
                return name + "=[...]";
            };
            return "BettiMatchingResult(" +
                   (reprMemberVector("input1_matched_birth_coordinates",
                                     self.input1MatchedBirthCoordinates) +
                    ", " +
                    reprMemberVector("input1_matched_death_coordinates",
                                     self.input1MatchedDeathCoordinates) +
                    ", " +
                    reprMemberVector("input2_matched_birth_coordinates",
                                     self.input2MatchedBirthCoordinates) +
                    ", " +
                    reprMemberVector("input2_matched_death_coordinates",
                                     self.input2MatchedDeathCoordinates) +
                    ", " +
                    (self.input1UnmatchedBirthCoordinates.has_value()
                         ? (reprMemberVector(
                                "input1_unmatched_birth_coordinates",
                                *self.input1UnmatchedBirthCoordinates) +
                            ", ")
                         : "") +
                    (self.input1UnmatchedDeathCoordinates.has_value()
                         ? (reprMemberVector(
                                "input1_unmatched_death_coordinates",
                                *self.input1UnmatchedDeathCoordinates) +
                            ", ")
                         : "") +
                    (self.input2UnmatchedBirthCoordinates.has_value()
                         ? (reprMemberVector(
                                "input2_unmatched_birth_coordinates",
                                *self.input2UnmatchedBirthCoordinates) +
                            ", ")
                         : "") +
                    (self.input2UnmatchedDeathCoordinates.has_value()
                         ? (reprMemberVector(
                                "input2_unmatched_death_coordinates",
                                *self.input2UnmatchedDeathCoordinates) +
                            ", ")
                         : "") +
                    reprSmall1DimMemberArray("num_matched", self.numMatched) +
                    (self.numUnmatchedInput1.has_value()
                         ? (", " +
                            reprSmall1DimMemberArray("num_unmatched_input1",
                                                     *self.numUnmatchedInput1))
                         : "") +
                    (self.numUnmatchedInput2.has_value()
                         ? (", " +
                            reprSmall1DimMemberArray("num_unmatched_input2",
                                                     *self.numUnmatchedInput2))
                         : "")) +
                   ")";
        });

    py::class_<BarcodeResult>(resultTypesModule, "BarcodeResult",
                              R"(
        BarcodeResult

        Holds the result of the barcode computation for a single input volume.

        Each coordinate array contains the birth or death coordinates of persistence
        pairs in all dimensions: first of the 0-dim. pairs, then of the 1-dim. pairs,
        and so on. The starting index of the persistence pairs birth/death coordinates
        in dimension i can be recovered using `np.cumsum()` on `num_matched`,
        `num_unmatched_input1`, and `num_unmatched_input2`,
        respectively.

        Attributes
        ----------
        birth_coordinates : numpy.ndarray
        death_coordinates : numpy.ndarray
            Arrays of shape (n_pairs, d). The birth and death coordinates of
            all persistence pairs in the input volume.
        num_pairs : numpy.ndarray
            Array of shape (d,). The number of persistence pairs by dimension.
        )")
        .def_readonly("birth_coordinates", &BarcodeResult::birthCoordinates)
        .def_readonly("death_coordinates", &BarcodeResult::deathCoordinates)
        .def_readonly("num_pairs", &BarcodeResult::numPairs)
        .def("__repr__", [](BarcodeResult &self) {
            auto reprSmall1DimMemberArray = [](string name,
                                               py::array_t<int64_t> &array) {
                vector<int64_t> values(array.size());
                for (int i = 0; i < array.size(); i++) {
                    values[i] = array.at(i);
                }
                return name + "=" + repr_vector(values, make_tuple("[", "]"));
            };
            auto reprMemberVector = [](string name,
                                       vector<py::array_t<int64_t>> &vector) {
                return name + "=[...]";
            };
            return "BarcodeResult(" +
                   (reprMemberVector("birth_coordinates",
                                     self.birthCoordinates) +
                    ", " +
                    reprMemberVector("death_coordinates",
                                     self.deathCoordinates) +
                    ", " +
                    reprSmall1DimMemberArray("num_pairs", self.numPairs)) +
                   ")";
        });

    py::class_<RepresentativeCycleResult>(resultTypesModule,
                                          "RepresentativeCycleResult",
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
        .def_readonly("matched_cycles",
                      &RepresentativeCycleResult::matchedCycles)
        .def_readonly("unmatched_cycles",
                      &RepresentativeCycleResult::unmatchedCycles)
        .def("__repr__", [](RepresentativeCycleResult &self) {
            return "RepresentativeCycleResult(" +
                   (self.matchedCycles.has_value()
                        ? ("matched_cycles=[...] (" +
                           std::to_string(self.matchedCycles->size()) +
                           " cycles)")
                        : "") +
                   (self.matchedCycles.has_value() &&
                            self.unmatchedCycles.has_value()
                        ? ", "
                        : "") +
                   (self.unmatchedCycles.has_value()
                        ? ("unmatched_cycles=[...] (" +
                           std::to_string(self.unmatchedCycles->size()) +
                           " cycles)")
                        : "") +
                   ")";
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
        modules that come from sublevel sets (i.e. the filtration starts at low values),
        and uses a x-y-z-type tiebreaking order (or x-y-type/x-type in 2D/1D)
        in the barcode computation algorithms.
        The essential interval [min(volume), +oo) that corresponds to the oldest
        connected component is not included in the results.

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
