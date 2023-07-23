#pragma once

#include "data_structures.h"

#include <string>

using namespace std;

enum file_format {DIPHA, PERSEUS, NUMPY};

struct Config {
	string filename_0 = "";
	string filename_1 = "";
	string matched_filename = "../matched.csv";
	string unmatched_0_filename = "../unmatched_0.csv";
	string unmatched_1_filename = "../unmatched_1.csv";
	file_format format_0;
	file_format format_1;

	value_t threshold = numeric_limits<value_t>::infinity();

	index_t min_recursion_to_cache = 0;
	index_t cache_size = 1 << 31;

	bool print = false;
	bool verbose = false;
};
