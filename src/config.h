#pragma once

#include <string>
#include <limits>
#include <cstdint>

#define RUNTIME
#define COMPUTE_COMPARISON
//#define USE_REDUCTION_MATRIX
#define USE_CACHE
#define USE_APPARENT_PAIRS
#define USE_APPARENT_PAIRS_COMP
#define USE_EMERGENT_PAIRS
#define USE_CLEARING_IMAGE
#define USE_CLEARING_DIM0

typedef uint32_t index_t;
typedef double value_t;

#define INFTY numeric_limits<value_t>::infinity()
#define NONE numeric_limits<index_t>::max()

enum fileFormat { DIPHA, PERSEUS, NUMPY };

using namespace std;


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
	size_t cacheSize = numeric_limits<size_t>::max();

	bool print = false;
	bool saveResult = false;
};
