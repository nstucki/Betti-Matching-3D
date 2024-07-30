#pragma once

#include "config.h"
#include "src_3D/data_structures.h"

#include <vector>
#include <unordered_map>
#include <optional>
#include <functional>

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

// Function to convert a two lists of pair indices (relative to the matched and unmatched pairs) to a list of pair references. The function is templated to work with 1D, 2D and 3D BettiMatching classes,
// which have their own Pair and Match types, respectively.
template<typename Pair, typename Match>
vector<std::reference_wrapper<Pair>> assembleRequestedPairs(const optional<vector<size_t>> &matchedPairsIndices, const optional<vector<size_t>> &unmatchedPairsIndices, vector<Pair> &pairs, unordered_map<uint64_t, bool> &isMatched,
        vector<Match> &matches, const int input) {
    // Assemble the list of requested pairs: First the matched pairs (all if empty optional was passed, then the unmatched pairs (all if empty optional was passed)
    vector<std::reference_wrapper<Pair>> requestedPairs;

    if (matchedPairsIndices.has_value()) {
        for (size_t i : *matchedPairsIndices) {
            if (i >= matches.size()) {
                throw invalid_argument("Index " + to_string(i) + " is out of bounds for matched pairs");
            }
            requestedPairs.push_back(input == 0 ? matches[i].pair0 : matches[i].pair1);
        }
    } else {
        for (auto &match : matches) {
            requestedPairs.push_back(input == 0 ? match.pair0 : match.pair1);
        }
    }

    if (unmatchedPairsIndices.has_value()) {
        // As unmatched pairs indices are relative to the unmatched pairs only, which are not available as a list,
        // we need to iterate over all pairs to find the requested ones
        auto nextUnmatchedPairIndex = unmatchedPairsIndices->begin();
        int unmatchedPairsSeen = 0;
        for (auto pair = pairs.begin(); pair != pairs.end() && nextUnmatchedPairIndex != unmatchedPairsIndices->end(); pair++) {
            if (!isMatched[pair->birth.index]) {
                if (unmatchedPairsSeen == *nextUnmatchedPairIndex) {
                    requestedPairs.push_back(*pair);
                    nextUnmatchedPairIndex++;
                }
                unmatchedPairsSeen++;
            }
        }
        if (nextUnmatchedPairIndex != unmatchedPairsIndices->end()) {
            throw invalid_argument("Index " + to_string(*nextUnmatchedPairIndex) + " is out of bounds for unmatched pairs");
        }
    } else {
        for (auto &pair : pairs) {
            if (!isMatched[pair.birth.index]) {
                requestedPairs.push_back(pair);
            }
        }
    }

    return requestedPairs;
}
