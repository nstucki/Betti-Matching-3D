#pragma once

#include "data_structures.h"

#include <string>

using namespace std;

enum fileFormat {DIPHA, PERSEUS, NUMPY};

struct Config {
	string filename_0 = "";
	string filename_1 = "";
	string matched_filename = "matched.csv";
	string unmatched_0_filename = "unmatched_0.csv";
	string unmatched_1_filename = "unmatched_1.csv";
	fileFormat format_0;
	fileFormat format_1;

	float threshold = numeric_limits<float>::infinity();
	index_t minRecursionToCache = 0;
	index_t cacheSize = numeric_limits<index_t>::max();
	
	bool print = false;
	bool verbose = false;
};
