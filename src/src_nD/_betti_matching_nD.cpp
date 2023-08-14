#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <cassert>
#include <utility>
#include <string>
#include "utils.h"
#include "config.h"
#include "data_structures.h"
#include "betti_matching_nD.h"

namespace py = pybind11;

std::string repr_coordinates(std::vector<index_t> coordinates)
{
    std::stringstream out_stream;
    out_stream << "(";
    std::for_each(coordinates.begin(), std::prev(coordinates.end()), [&out_stream, &coordinates](auto coordinate)
                  { out_stream << std::to_string(coordinate) << ", "; });
    out_stream << coordinates[coordinates.size() - 1] << ")";
    return out_stream.str();
};

std::string repr_cube(Cube &self)
{
    return "Cube(birth=" + std::to_string(self.birth) + ", coordinates=" + repr_coordinates(self.coordinates) + ")";
};

std::string repr_pair(Pair &self)
{
    return "Pair(birth=" + repr_cube(self.birth) + ", death=" + repr_cube(self.death) + ")";
};

std::string repr_match(Match &self)
{
    return "Match((" + std::to_string(self.pair0.birth.birth) + ", " + std::to_string(self.pair0.death.birth) + ") -> (" + std::to_string(self.pair1.birth.birth) + ", " + std::to_string(self.pair1.death.birth) + "))";
};

PYBIND11_MODULE(betti_matching, m)
{
    m.def(
        "read_image",
        [](string const &filename, fileFormat const &format)
        {
            vector<double> image;
            vector<index_t> shape;
            readImage(filename, format, image, shape);
            auto image_array = py::array_t<double>(
                shape,
                &image[0]);
            return image_array;
        });
    m.attr("file_format") =
        &py::enum_<fileFormat>(m, "file_format")
             .value("DIPHA", fileFormat::DIPHA)
             .value("PERSEUS", fileFormat::PERSEUS)
             .value("NUMPY", fileFormat::NUMPY);

    py::class_<CubicalGridComplex>(m, "CubicalGridComplex")
        .def(py::init<const vector<value_t>, const vector<index_t> &>())
        .def(py::init([](py::array_t<double> &image)
                      {
                double *image_pointer = image.mutable_data();
                // TODO this is probably not memory-safe, right? We'd need to std::move it?
                vector<double> _image(image_pointer, image_pointer + image.size());
                const vector<index_t> _shape(image.shape(), image.shape() + image.ndim());

                return std::move(CubicalGridComplex(_image, _shape)); }))
        .def("get_parent_voxel", &CubicalGridComplex::getParentVoxel)
        .def("print_cubes", &CubicalGridComplex::printCubes);

    py::class_<Config>(m, "Config")
        .def(py::init<>())
        .def_readwrite("filename_0", &Config::filename_0)
        .def_readwrite("filename_1", &Config::filename_1)
        .def_readwrite("matched_filename", &Config::matched_filename)
        .def_readwrite("unmatched_0_filename", &Config::unmatched_0_filename)
        .def_readwrite("unmatched_1_filename", &Config::unmatched_1_filename)
        .def_readwrite("format_0", &Config::format_0)
        .def_readwrite("format_1", &Config::format_1)
        .def_readwrite("threshold", &Config::threshold)
        .def_readwrite("min_recursion_to_cache", &Config::minRecursionToCache)
        .def_readwrite("cache_size", &Config::cacheSize)
        .def_readwrite("print", &Config::print)
        .def_readwrite("verbose", &Config::verbose);

    py::class_<BettiMatchingND>(m, "BettiMatchingND")
        .def(py::init<
             const CubicalGridComplex &,
             const CubicalGridComplex &,
             const CubicalGridComplex &,
             const Config &>())
        .def("compute_matching", &BettiMatchingND::computeMatching)
        .def_readonly("dim", &BettiMatchingND::dim)
        .def_property_readonly("cgc_0", [](BettiMatchingND &self)
                               { return self.cgc0; })
        .def_property_readonly("cgc_1", [](BettiMatchingND &self)
                               { return self.cgc1; })
        .def_property_readonly("cgc_comp", [](BettiMatchingND &self)
                               { return self.cgcComp; })
        .def_property_readonly("pairs_0", [](BettiMatchingND &self)
                               { return self.pairs0; })
        .def_property_readonly("pairs_1", [](BettiMatchingND &self)
                               { return self.pairs1; })
        .def_property_readonly("pairs_comp", [](BettiMatchingND &self)
                               { return self.pairsComp; })
        .def_property_readonly("matches", [](BettiMatchingND &self)
                               { return self.matches; })
        .def_property_readonly("is_matched_0", [](BettiMatchingND &self)
                               { return self.isMatched0; })
        .def_property_readonly("is_matched_1", [](BettiMatchingND &self)
                               { return self.isMatched1; })
        .def_property_readonly("unmatched_pairs_0", [](BettiMatchingND &self)
                               { return self.unmatchedPairs0; })
        .def_property_readonly("unmatched_pairs_1", [](BettiMatchingND &self)
                               { return self.unmatchedPairs1; });

    py::class_<Pair>(m, "Pair")
        .def_readonly("birth", &Pair::birth)
        .def_readonly("death", &Pair::death)
        .def("__repr__", &repr_pair);

    py::class_<Match>(m, "Match")
        .def_readonly("pair0", &Match::pair0)
        .def_readonly("pair1", &Match::pair1)
        .def("__repr__", &repr_match);

    py::class_<Cube>(m, "Cube")
        .def_readonly("birth", &Cube::birth)
        .def_readonly("coordinates", &Cube::coordinates)
        .def("__repr__", &repr_cube);

    m.def("print_result", &printResult);
}
