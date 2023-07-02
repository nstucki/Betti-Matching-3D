#pragma once

#include "config.h"

#include <vector>

using namespace std;


vector<vector<bool>> getSubsets(uint64_t n, uint64_t k);

void readImage(string const &filename, fileFormat const &format, vector<double> &image, vector<uint64_t> &shape);