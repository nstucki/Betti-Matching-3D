#pragma once

//#include <cfloat>
#include <string>

using namespace std;

enum output_location { LOC_NONE, LOC_YES};
enum file_format { DIPHA, PERSEUS, NUMPY, CSV };

struct Config {
	string filename_0 = "";
	string filename_1 = "";
	string output_filename = "output.csv";
	file_format format_0;
	file_format format_1;
	bool print = false;
	bool verbose = false;
	int min_recursion_to_cache = 0;
	uint32_t cache_size = 1 << 31;
};
