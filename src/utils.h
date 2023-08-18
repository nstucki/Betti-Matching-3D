#pragma once

#include "config.h"
#include "src_3D/data_structures.h"

#include <vector>
#include <unordered_map>

using namespace std;

void readImage(const string& filename, const fileFormat& format, vector<double>& image, vector<index_t>& shape);

vector<vector<bool>> getSubsets(index_t n, index_t k);