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
	string filename0 = "";
	string filename1 = "";
	string matchedFilename = "../matched.csv";
	string unmatched0Filename = "../unmatched_0.csv";
	string unmatched1Filename = "../unmatched_1.csv";
	fileFormat format0;
	fileFormat format1;

	value_t threshold = numeric_limits<value_t>::infinity();

	index_t minRecursionToCache = 0;
	index_t cacheSize = 1 << 31;

	bool print = false;
	bool verbose = false;
};
