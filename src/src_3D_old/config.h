#pragma once

#include <string>

using namespace std;

enum output_location {LOC_NONE, LOC_YES};
enum file_format {DIPHA, PERSEUS, NUMPY};

struct Config {
	string filename_0 = "";
	string filename_1 = "";
	string matched_filename = "matched.csv";
	string unmatched_0_filename = "unmatched_0.csv";
	string unmatched_1_filename = "unmatched_1.csv";
	file_format format_0;
	file_format format_1;
	int min_recursion_to_cache = 0;
	uint64_t cache_size = 1 << 31;
	bool print = false;
	bool verbose = false;
};
