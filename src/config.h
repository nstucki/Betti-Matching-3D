#ifndef CONFIG_H
#define CONFIG_H

#include <cfloat>
#include <string>

enum calculation_method { LINKFIND, COMPUTEPAIRS, ALEXANDER};
enum output_location { LOC_NONE, LOC_YES};
enum file_format { DIPHA, PERSEUS, NUMPY, CSV };


struct Config {
	std::string filename = "";
	std::string output_filename = "output.csv"; //default output filename
	file_format format;
//    calculation_method method = ALEXANDER;
	calculation_method method = LINKFIND;
	double threshold = DBL_MAX;
	int maxdim=2;;  // compute PH for these dimensions
	bool print = false; // flag for printing parsistence pairs to stdout
	bool verbose = false;
	bool tconstruction = false; // T-construction or V-construction
	bool embedded = false; // embed image in the sphere (for alexander duality)
	output_location location = LOC_YES; // flag for saving location
	int min_recursion_to_cache = 0; // num of minimum recursions for a reduced column to be cached
	uint32_t cache_size = 1 << 31; // the maximum number of reduced columns to be cached
};

#endif
