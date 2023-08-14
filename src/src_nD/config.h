#pragma once

#include <string>
#include <limits>
#include <cstdint>

using namespace std;

#define USE_CLEARING_DIM_0
#define USE_DOUBLE
typedef uint32_t index_t;

#define INFTY numeric_limits<value_t>::infinity()
#define NONE numeric_limits<index_t>::max()
#ifdef USE_DOUBLE
typedef double value_t;
#else
typedef float value_t;
#endif

enum fileFormat {DIPHA, PERSEUS, NUMPY};


struct Config {
	string filename_0 = "";
	string filename_1 = "";
	string matched_filename = "../matched.csv";
	string unmatched_0_filename = "../unmatched_0.csv";
	string unmatched_1_filename = "../unmatched_1.csv";
	fileFormat format_0;
	fileFormat format_1;

	value_t threshold = INFTY;
	
	index_t minRecursionToCache = 0;
	index_t cacheSize = numeric_limits<index_t>::max();
	
	bool print = false;
	bool verbose = false;
};
