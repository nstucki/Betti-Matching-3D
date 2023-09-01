#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <cassert>
#include <utility>
#include <string>
#include "utils.h"
#include "config.h"
#include "data_structures.h"
#include "BettiMatching.h"

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

    m.def("compute_matching_with_voxels", [](py::array_t<value_t> &input0, py::array_t<value_t> &input1, py::array_t<value_t> &comparison,
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

        return BettiMatching(input0Vector, input1Vector, comparisonVector, shape0, config).computeMatchingWithVoxels();
    });

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
