#pragma once

#include <string>
#include <limits>
#include <cstdint>

#define RUNTIME
//#define COMPUTE_COMPARISON
//#define USE_REDUCTION_MATRIX
#define USE_APPARENT_PAIRS_COMP
#define USE_APPARENT_PAIRS
#define USE_EMERGENT_PAIRS
#define USE_CACHE
#define USE_CLEARING_IMAGE
#define USE_CLEARING_DIM0

typedef uint32_t index_t;
typedef double value_t;

#define INFTY numeric_limits<value_t>::infinity()
#define NONE numeric_limits<index_t>::max()

enum fileFormat { DIPHA, PERSEUS, NUMPY };

using namespace std;



struct Config {
	value_t threshold = INFTY;

	size_t minRecursionToCache = 1;
	size_t cacheSize = numeric_limits<size_t>::max();
};
