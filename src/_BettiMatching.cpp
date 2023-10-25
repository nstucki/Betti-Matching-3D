#include <algorithm>
#include <cstddef>
#include <iterator>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <cassert>
#include <future>
#include <utility>
#include <string>
#include "utils.h"
#include "config.h"
#include "data_structures.h"
#include "BettiMatching.h"

#include <iostream>

namespace py = pybind11;

std::string repr_vector(const std::vector<index_t> shape)
{
    std::stringstream out_stream;
    out_stream << "(";
    std::for_each(shape.begin(), std::prev(shape.end()), [&out_stream, &shape](auto size)
                  { out_stream << std::to_string(size) << ", "; });
    out_stream << shape[shape.size() - 1] << ")";
    return out_stream.str();
};

std::tuple<vector<vector<VoxelMatch>>, vector<vector<VoxelPair>>, vector<vector<VoxelPair>>> computeMatchingWithVoxelsFromPyArray(
    py::array_t<value_t> &input0, py::array_t<value_t> &input1, py::array_t<value_t> &comparison, Config &config)
{
    const vector<index_t> shape0(input0.shape(), input0.shape() + input0.ndim());
    const vector<index_t> shape1(input1.shape(), input1.shape() + input1.ndim());
    const vector<index_t> shapeComparison(comparison.shape(), comparison.shape() + comparison.ndim());
    if (shape0 != shape1 || shape0 != shapeComparison) {
        throw invalid_argument("The shapes of the tree input volumes must agree. Got " + repr_vector(shape0) + ", " + repr_vector(shape1) + " and " + repr_vector(shapeComparison));
    }
    const vector<value_t> input0Vector(input0.mutable_data(), input0.mutable_data() + input0.size());
    const vector<value_t> input1Vector(input1.mutable_data(), input1.mutable_data() + input1.size());
    const vector<value_t> comparisonVector(comparison.mutable_data(), comparison.mutable_data() + comparison.size());

    return BettiMatching(input0Vector, input1Vector, comparisonVector, shape0, config).computeMatchingWithVoxels();
};

PYBIND11_MODULE(betti_matching, m)
{
    py::class_<Config>(m, "Config")
        .def(py::init<>())
        .def_readwrite("filename_0", &Config::filename0)
        .def_readwrite("filename_1", &Config::filename1)
        .def_readwrite("matched_filename", &Config::matchedFilename)
        .def_readwrite("unmatched_0_filename", &Config::unmatched0Filename)
        .def_readwrite("unmatched_1_filename", &Config::unmatched1Filename)
        .def_readwrite("format_0", &Config::format0)
        .def_readwrite("format_1", &Config::format1)
        .def_readwrite("threshold", &Config::threshold)
        .def_readwrite("min_recursion_to_cache", &Config::minRecursionToCache)
        .def_readwrite("cache_size", &Config::cacheSize);

    py::class_<BettiMatching>(m, "BettiMatching")
        .def(py::init([](py::array_t<value_t> &input0, py::array_t<value_t> &input1, py::array_t<value_t> &comparison,
                         Config &config)
                      {
                const vector<index_t> shape0(input0.shape(), input0.shape() + input0.ndim());
                const vector<index_t> shape1(input1.shape(), input1.shape() + input1.ndim());
                const vector<index_t> shapeComparison(comparison.shape(), comparison.shape() + comparison.ndim());
                if (shape0 != shape1 || shape0 != shapeComparison) {
                    throw invalid_argument("The shapes of the tree input volumes must agree. Got " + repr_vector(shape0) + ", " + repr_vector(shape1) + " and " + repr_vector(shapeComparison));
                }
                const vector<value_t> input0Vector(input0.mutable_data(), input0.mutable_data() + input0.size());
                const vector<value_t> input1Vector(input1.mutable_data(), input1.mutable_data() + input1.size());
                const vector<value_t> comparisonVector(comparison.mutable_data(), comparison.mutable_data() + comparison.size());

                return BettiMatching(input0Vector, input1Vector, comparisonVector, shape0, config); }))
        .def(py::init([](std::string input0_path, std::string input1_path, std::string comparison_path,
                         Config &config)
                      {

                vector<value_t> input0Vector;
                vector<value_t> input1Vector;
                vector<value_t> comparisonVector;

                vector<index_t> shape0;
                vector<index_t> shape1;
                vector<index_t> shapeComparison;

                readImage(input0_path, NUMPY, input0Vector, shape0);
                readImage(input1_path, NUMPY, input1Vector, shape1);
                readImage(comparison_path, NUMPY, comparisonVector, shapeComparison);

                if (shape0 != shape1 || shape0 != shapeComparison) {
                    throw invalid_argument("The shapes of the tree input volumes must agree. Got " + repr_vector(shape0) + ", " + repr_vector(shape1) + " and " + repr_vector(shapeComparison));
                }

                return BettiMatching(input0Vector, input1Vector, comparisonVector, shape0, config); }))
        .def("compute_matching_with_voxels",
             &BettiMatching::computeMatchingWithVoxels);

    m.def("compute_matching_with_voxels", computeMatchingWithVoxelsFromPyArray);

    m.def("compute_matching_with_voxels", [](std::string input0_path, std::string input1_path, std::string comparison_path,
                         Config &config)
    {

        vector<value_t> input0Vector;
        vector<value_t> input1Vector;
        vector<value_t> comparisonVector;

        vector<index_t> shape0;
        vector<index_t> shape1;
        vector<index_t> shapeComparison;

        readImage(input0_path, NUMPY, input0Vector, shape0);
        readImage(input1_path, NUMPY, input1Vector, shape1);
        readImage(comparison_path, NUMPY, comparisonVector, shapeComparison);

        if (shape0 != shape1 || shape0 != shapeComparison) {
            throw invalid_argument("The shapes of the tree input volumes must agree. Got " + repr_vector(shape0) + ", " + repr_vector(shape1) + " and " + repr_vector(shapeComparison));
        }

        return BettiMatching(input0Vector, input1Vector, comparisonVector, shape0, config).computeMatchingWithVoxels();
    });

    m.def("compute_matching_with_voxels", [](
        vector<py::array_t<value_t>> inputs0,
        vector<py::array_t<value_t>> inputs1,
        vector<py::array_t<value_t>> comparisons,
        Config &config)
    {
        if (inputs0.size() != inputs1.size() || inputs0.size() != comparisons.size()) {
            throw invalid_argument("Different numbers of inputs where provided: " + std::to_string(inputs0.size()) + " (inputs0), " + std::to_string(inputs1.size()) + " (inputs1) and " + std::to_string(comparisons.size()) + " (comparisons)");
        }
        size_t batch_size = inputs0.size();
        size_t num_dimensions = inputs0[0].ndim();

        std::vector<std::future<std::tuple<vector<vector<VoxelMatch>>, vector<vector<VoxelPair>>, vector<vector<VoxelPair>>>>> result_futures;
        for (int i = 0; i < inputs0.size(); i++) {
            result_futures.push_back(std::async(computeMatchingWithVoxelsFromPyArray, std::ref(inputs0[i]), std::ref(inputs1[i]), std::ref(comparisons[i]), std::ref(config)));
        }

        // Now block on all of them one at a time.
        std::vector<std::tuple<vector<vector<VoxelMatch>>, vector<vector<VoxelPair>>, vector<vector<VoxelPair>>>> results;
        for (auto& future : result_futures) {
            results.push_back(future.get());
        }

        // vector<py::array_t<double>> all_prediction_matches_birth_coordinates; // vector of size batch_size, elements of shape (N_k, 3)
        // vector<py::array_t<double>> all_prediction_matches_death_coordinates; // vector of size batch_size, elements of shape (N_k, 3)
        // vector<py::array_t<double>> all_target_matches_birth_coordinates; // vector of size batch_size, elements of shape (M_k, 3)
        // vector<py::array_t<double>> all_target_matches_death_coordinates; // vector of size batch_size, elements of shape (M_k, 3)
        // vector<py::array_t<double>> all_prediction_unmatched_birth_coordinates; // vector of size batch_size, elements of shape (L_k, 3)
        // vector<py::array_t<double>> all_prediction_unmatched_death_coordinates; // vector of size batch_size, elements of shape (L_k, 3)
        vector<std::tuple<py::array_t<int64_t>, py::array_t<int64_t>, py::array_t<int64_t>, py::array_t<int64_t>, py::array_t<int64_t>, py::array_t<int64_t>>> resultArrays;

        for (auto &result : results) {
            auto matches_by_dimension = std::get<0>(result);
            auto prediction_unmatched_by_dimension = std::get<1>(result);

            size_t num_matches = 0;
            size_t num_prediction_unmatched = 0;
            // cout << prediction_unmatched_by_dimension[0].size() << " " << prediction_unmatched_by_dimension[1].size() << " " << prediction_unmatched_by_dimension[2].size() << endl;
            for (size_t d = 0; d < num_dimensions; d++) {
                num_matches += matches_by_dimension[d].size();
                num_prediction_unmatched += prediction_unmatched_by_dimension[d].size();
            }
            size_t matches_coordinates_shape[2]{num_matches, num_dimensions};
            size_t matches_coordinates_stride[2]{num_dimensions * sizeof(int64_t), sizeof(int64_t)};
            size_t prediction_unmatched_coordinates_shape[2]{num_prediction_unmatched, num_dimensions};
            size_t prediction_unmatched_stride[2]{num_dimensions * sizeof(int64_t), sizeof(int64_t)};
            
            py::array_t<int64_t> prediction_matches_birth_coordinates(matches_coordinates_shape, matches_coordinates_stride); // shape (N_k, 3)
            py::array_t<int64_t> prediction_matches_death_coordinates(matches_coordinates_shape, matches_coordinates_stride); // shape (N_k, 3)
            py::array_t<int64_t> target_matches_birth_coordinates(matches_coordinates_shape, matches_coordinates_stride); // shape (M_k, 3)
            py::array_t<int64_t> target_matches_death_coordinates(matches_coordinates_shape, matches_coordinates_stride); // shape (M_k, 3)
            py::array_t<int64_t> prediction_unmatched_birth_coordinates(prediction_unmatched_coordinates_shape, prediction_unmatched_stride); // shape (L_k, 3)
            py::array_t<int64_t> prediction_unmatched_death_coordinates(prediction_unmatched_coordinates_shape, prediction_unmatched_stride); // shape (L_k, 3)

            auto prediction_matches_birth_coordinates_view = prediction_matches_birth_coordinates.mutable_unchecked();
            auto prediction_matches_death_coordinates_view = prediction_matches_death_coordinates.mutable_unchecked();
            auto target_matches_birth_coordinates_view = target_matches_birth_coordinates.mutable_unchecked();
            auto target_matches_death_coordinates_view = target_matches_death_coordinates.mutable_unchecked();
            auto prediction_unmatched_birth_coordinates_view = prediction_unmatched_birth_coordinates.mutable_unchecked();
            auto prediction_unmatched_death_coordinates_view = prediction_unmatched_death_coordinates.mutable_unchecked();

            int i = 0;
            for (auto &matches_in_dimension : matches_by_dimension) {
                for (auto &match : matches_in_dimension) {
                    for (int d = 0; d < num_dimensions; d++) {
                        prediction_matches_birth_coordinates_view(i, d) = match.pair0.birth[d];
                        prediction_matches_death_coordinates_view(i, d) = match.pair0.death[d];
                        target_matches_birth_coordinates_view(i, d) = match.pair1.birth[d];
                        target_matches_death_coordinates_view(i, d) = match.pair1.death[d];
                    }
                    i += 1;
                }
            }
            i = 0;
            for (auto &unmatched_in_dimension : prediction_unmatched_by_dimension) {
                for (auto &unmatched : unmatched_in_dimension) {
                    for (int d = 0; d < num_dimensions; d++) {
                        prediction_unmatched_birth_coordinates_view(i, d) = unmatched.birth[d];
                        prediction_unmatched_death_coordinates_view(i, d) = unmatched.death[d];
                    }
                    i += 1;
                }
            }

            std::tuple<py::array_t<int64_t>, py::array_t<int64_t>, py::array_t<int64_t>, py::array_t<int64_t>, py::array_t<int64_t>, py::array_t<int64_t>> resultsTuple = {
                std::move(prediction_matches_birth_coordinates),
                std::move(prediction_matches_death_coordinates),
                std::move(target_matches_birth_coordinates),
                std::move(target_matches_death_coordinates),
                std::move(prediction_unmatched_birth_coordinates),
                std::move(prediction_unmatched_death_coordinates)};
            resultArrays.emplace_back(resultsTuple);
        }
        return resultArrays;
    });

    py::class_<VoxelMatch>(m, "VoxelMatch")
        .def_readonly("pair0", &VoxelMatch::pair0)
        .def_readonly("pair1", &VoxelMatch::pair1)
        .def("__repr__", [](VoxelMatch &self)
             { return "VoxelMatch(pair0=(birth=" + repr_vector(self.pair0.birth) + ", death=" + repr_vector(self.pair0.death) + "), pair1=(birth=" + repr_vector(self.pair1.birth) + ", death=" + repr_vector(self.pair1.death) + ")"; });

    py::class_<VoxelPair>(m, "VoxelPair")
        .def_readonly("birth", &VoxelPair::birth)
        .def_readonly("death", &VoxelPair::death)
        .def("__repr__", [](VoxelPair &self) {
          return "VoxelPair(birth=" + repr_vector(self.birth) +
                 ", death=" + repr_vector(self.death) + ")";
        });
}
