#pragma once

#include "config.h"
#include "src_3D/data_structures.h"

#include <vector>
#include <unordered_map>

using namespace std;

void readImage(const string& filename, const fileFormat& format, vector<double>& image, vector<index_t>& shape);

vector<vector<bool>> getSubsets(index_t n, index_t k);

// See https://stackoverflow.com/a/28411055
template <typename T, std::size_t... Indices>
auto vectorToTupleHelper(const std::vector<T> &v,
                         std::index_sequence<Indices...>) {
        return std::make_tuple(v[Indices]...);
}

template <std::size_t N, typename T>
auto vectorToTuple(const std::vector<T> &v) {
  return vectorToTupleHelper(v, std::make_index_sequence<N>());
}

// See https://stackoverflow.com/a/42495119
template <class Tuple,
   class T = std::decay_t<std::tuple_element_t<0, std::decay_t<Tuple>>>>
std::vector<T> tupleToVector(Tuple&& tuple)
{
    return std::apply([](auto&&... elems) {
        std::vector<T> result;
        result.reserve(sizeof...(elems));
        (result.push_back(std::forward<decltype(elems)>(elems)), ...);
        return result;
    }, std::forward<Tuple>(tuple));
}
