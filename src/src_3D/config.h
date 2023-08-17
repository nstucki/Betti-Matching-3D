#pragma once

#include <string>
#include <limits>
#include <cstdint>

using namespace std;

#define RUNTIME
//#define COMPUTE_COMPARISON

#define USE_EMERGENT_PAIRS
#define USE_CLEARING_DIM_0

#define INFTY numeric_limits<value_t>::infinity()
#define NONE numeric_limits<index_t>::max()

typedef uint32_t index_t;
typedef double value_t;

enum fileFormat { DIPHA, PERSEUS, NUMPY };


struct Config {
	string filename0 = "";
	string filename1 = "";
	string matchedFilename = "matched.csv";
	string unmatched0Filename = "unmatched_0.csv";
	string unmatched1Filename = "unmatched_1.csv";
	fileFormat format0;
	fileFormat format1;

	value_t threshold = numeric_limits<value_t>::infinity();

	size_t minRecursionToCache = 1;
	size_t cacheSize = 1 << 31;

	bool verbose = false;
	bool saveResult = false;
};
