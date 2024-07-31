#pragma once

#include <cstdint>
#include <limits>
#include <string>

// #define RUNTIME
#define COMPUTE_COMPARISON
// #define USE_ISPAIRED
// #define USE_REDUCTION_MATRIX
// #define USE_APPARENT_PAIRS_COMP
// #define USE_APPARENT_PAIRS
#define USE_EMERGENT_PAIRS
#define USE_CACHE
// #define USE_CLEARING_IMAGE
#define USE_CLEARING_DIM0
#define USE_STABLE_SORT_OR_STABLE_PARTITION // enables stable sort/binary input
                                            // sorting optimization in dim3 edge
                                            // enumeration methods
#define PARALLELIZE_INDEPENDENT_BARCODES_DIM1 // enables parallelization of
                                              // independent barcode
                                              // computations in
                                              // dim3::Dimension1

typedef uint32_t index_t;
typedef double value_t;

#define INFTY numeric_limits<value_t>::infinity()
#define NONE numeric_limits<index_t>::max()
#define NONE_INDEX numeric_limits<uint64_t>::max()

enum fileFormat { DIPHA, PERSEUS, NUMPY };

using namespace std;

#if defined(PARALLELIZE_INDEPENDENT_BARCODES_DIM1) and                         \
    defined(USE_CLEARING_IMAGE)
static_assert(false,
              "PARALLELIZE_INDEPENDENT_BARCODES_DIM1 is active alongside "
              "incompatible options (race conditions could occur!)");
#endif

struct Config {
    value_t threshold = INFTY;

    size_t minRecursionToCache = 1;
    size_t cacheSize = numeric_limits<size_t>::max();
};
