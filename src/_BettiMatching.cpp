#include "utils.h"
#include "config.h"
#include "data_structures.h"
#include "BettiMatching.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <cassert>
#include <utility>
#include <string>

using namespace std;
namespace py = pybind11;



string repr_vector(const vector<index_t> shape) {
    stringstream out_stream;
    out_stream << "(";
    for_each(shape.begin(), std::prev(shape.end()), [&out_stream, &shape](auto size)
                { out_stream << std::to_string(size) << ", "; });
    out_stream << shape[shape.size() - 1] << ")";
    return out_stream.str();
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


    m.def("compute_matching", [](py::array_t<value_t>& input0, py::array_t<value_t>& input1) {
        vector<index_t> shape0(input0.shape(), input0.shape() + input0.ndim());
        vector<index_t> shape1(input1.shape(), input1.shape() + input1.ndim());
        if (shape0 != shape1) {
            throw invalid_argument("The shapes of the tree input volumes must agree. Got " + repr_vector(shape0)
                                    + " and " + repr_vector(shape1));
        }
        vector<value_t> input0Vector(input0.mutable_data(), input0.mutable_data() + input0.size());
        vector<value_t> input1Vector(input1.mutable_data(), input1.mutable_data() + input1.size());
        Config config;
        BettiMatching BM(std::move(input0Vector), std::move(input1Vector), std::move(shape0), std::move(config));
        BM.computeMatching();
        return BM.getMatching();
    });


    m.def("compute_matching", [](string input0_path, string input1_path) {
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
        BettiMatching BM(std::move(input0Vector), std::move(input1Vector), std::move(shape0), std::move(config));
        BM.computeMatching();
        return BM.getMatching();
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
}
