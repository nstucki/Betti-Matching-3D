#pragma once

#include "config.h"
#include "data_structures.h"

#include <vector>
#include <unordered_map>

using namespace std;


vector<vector<bool>> getSubsets(index_t n, index_t k);

void readImage(string const &filename, fileFormat const &format, vector<double> &image, vector<index_t> &shape);

void printResult(index_t dim, const CubicalGridComplex& cgc0, const CubicalGridComplex& cgc1, const CubicalGridComplex& cgcComp, 
    const vector<vector<Pair>>& pairs0, const vector<vector<Pair>>& pairs1, const vector<vector<Pair>>& pairsComp, 
    const vector<vector<Match>>& matches, unordered_map<index_t, bool>& isMatched0, unordered_map<index_t, bool>& isMatched1);