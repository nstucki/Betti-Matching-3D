#pragma once

#include "config.h"
#include "data_structures.h"

#include <vector>
#include <unordered_map>

using namespace std;

void readImage(const string& filename, const fileFormat& format, vector<double>& image, vector<index_t>& shape);
void tokenize(const string& str, const char delim, vector<string> &out);
void printResult(const CubicalGridComplex* const cgc0, const CubicalGridComplex* const cgc1, const CubicalGridComplex* const cgcComp, 
                    const vector<vector<Pair>>& pairs0, const vector<vector<Pair>>& pairs1, const vector<vector<Pair>>& pairsComp,
                    const vector<vector<Match>>& matched, 
                    unordered_map<uint64_t, bool>& isMatched0, unordered_map<uint64_t, bool>& isMatched1);