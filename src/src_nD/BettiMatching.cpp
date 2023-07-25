#include "template_functions.h"

#include "utils.h"
#include "top_dimension.h"
#include "inter_dimensions.h"
#include "dimension_0.h"

#include <iostream>
#include <fstream>
#include <chrono>

using namespace std;
using namespace std::chrono;


void print_usage_and_exit(int exit_code) {
    cout << endl;
    cerr << "Usage: "
         << "Betti Matching "
         << "[options] [input_0_filename] [input_1_filename]" << endl
         << endl
         << "Options:" << endl
         << endl
         << "  --help, -h                      print this screen" << endl
         << "  --threshold, -t                 compute persistent homology of up to threshold" << endl
         << "  --verbose, -v                   display detailed processing information" << endl
         << "  --min_recursion_to_cache, -mc   minimum number of recursion for a reduced column to be cached" << endl
         << "  --cache_size, -c	          maximum number of reduced columns to be cached" << endl
         << "  --matched, -m                   name of the file containing matched pairs" << endl
         << "  --unmatchedmatched_0, -u0       name of the file containing unmatched pairs of input 0" << endl
         << "  --unmatchedmatched_1, -u1       name of the file containing unmatched pairs of input 1" << endl
         << "  --print, -p                     print result in console" << endl
         << endl;
	exit(exit_code);
}


int main(int argc, char** argv) {
    Config config;
    
    for (int i = 1; i < argc; ++i) {
		const string arg(argv[i]);
		if (arg == "--help" || arg == "-h") {
			print_usage_and_exit(0);
        } else if (arg == "--threshold" || arg == "-t") {
            string parameter = string(argv[++i]);
			size_t next_pos;
			config.threshold = stod(parameter, &next_pos);
			if (next_pos != parameter.size()) print_usage_and_exit(-1);
		} else if (arg == "--verbose" || arg == "-v") {
			config.verbose = true;
		} else if (arg == "--matched" || arg == "-m") {
			config.matched_filename = string(argv[++i]);
        } else if (arg == "--unmatched_0" || arg == "-u0") {
			config.unmatched_0_filename = string(argv[++i]);
        } else if (arg == "--unmatched_1" || arg == "-u1") {
			config.unmatched_1_filename = string(argv[++i]);
		} else if (arg == "--min_recursion_to_cache" || arg == "-mc") {
            config.minRecursionToCache = stoi(argv[++i]);
		} else if (arg == "--cache_size" || arg == "-c") {
            config.cacheSize = stoi(argv[++i]);
		} else if (arg == "--print" || arg == "-p") {
			config.print = true;
		} else {
            if (config.filename_0.empty()) {
                config.filename_0 = argv[i];
            } else {
                config.filename_1 = argv[i];
            }
		}
	}
    if (config.filename_0.empty()) { print_usage_and_exit(-1); } 
    if (config.filename_0.find(".txt") != string::npos) {
		config.format_0 = PERSEUS;
	} else if (config.filename_0.find(".npy")!= string::npos) {
		config.format_0 = NUMPY;
	} else if (config.filename_0.find(".complex")!= string::npos) {
		config.format_0 = DIPHA;
	} else {
		cerr << "unknown input file format! (the filename extension should be .txt/.npy/.complex): " << config.filename_0 << endl;
		exit(-1);
	}
    ifstream file_stream_0(config.filename_0);
	if (!config.filename_0.empty() && file_stream_0.fail()) {
		cerr << "couldn't open file " << config.filename_0 << endl;
		exit(-1);
	}
    if (config.filename_1.empty()) { print_usage_and_exit(-1); }
    if(config.filename_1.find(".txt")!= string::npos) {
		config.format_1 = PERSEUS;
	} else if (config.filename_1.find(".npy")!= string::npos) {
		config.format_1 = NUMPY;
	} else if (config.filename_1.find(".complex")!= string::npos) {
		config.format_1 = DIPHA;
	} else {
		cerr << "unknown input file format! (the filename extension should be .txt/.npy/.complex): " << config.filename_1 << endl;
		exit(-1);
	}
    ifstream file_stream_1(config.filename_1);
	if (!config.filename_1.empty() && file_stream_1.fail()) {
		cerr << "couldn't open file " << config.filename_1 << endl;
		exit(-1);
	}

    if (config.verbose) { cout << "reading images ... "; }
    auto start = high_resolution_clock::now();

    vector<double> readImage0;
    vector<double> readImage1;
    vector<index_t> shape0;
    vector<index_t> shape1;
    readImage(config.filename_0, config.format_0, readImage0, shape0);
    readImage(config.filename_1, config.format_1, readImage1, shape1);
    #ifdef USE_FLOAT
    vector<value_t> image0(readImage0.begin(), readImage0.end());
    vector<value_t> image1(readImage1.begin(), readImage1.end());
    #endif
    #ifdef USE_DOUBLE
    vector<value_t>& image0 = readImage0;
    vector<value_t>& image1 = readImage1;
    #endif

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);
    if (config.verbose) { cout << "took " << duration.count() << " ms" << endl; }
    
    assert (shape0 == shape1);
    index_t dim = shape0.size();
    
    vector<value_t> imageComp;
    transform(image0.begin(), image0.end(), image1.begin(), back_inserter(imageComp), [](value_t a, value_t b) { return min(a,b); });

    CubicalGridComplex cgc0(move(image0), shape0);
    CubicalGridComplex cgc1(move(image1), shape1);
    CubicalGridComplex cgcComp(move(imageComp), shape0);

    vector<vector<Pair>> pairs0(dim);
	vector<vector<Pair>> pairs1(dim);
	vector<vector<Pair>> pairsComp(dim);
	vector<vector<Match>> matches(dim);
	unordered_map<index_t, bool> isMatched0;
	unordered_map<index_t, bool> isMatched1;

 	vector<Cube> ctr0;
 	vector<Cube> ctr1;
    vector<Cube> ctrComp;

    if (dim > 1) {
        TopDimension topDim(cgc0, cgc1, cgcComp,  config, pairs0[dim-1], pairs1[dim-1], pairsComp[dim-1], matches[dim-1],
                            isMatched0, isMatched1);       
        topDim.computePairsAndMatch(ctr0, ctr1, ctrComp);
    }

    if (dim > 2) {
        InterDimensions interDim(cgc0, cgc1, cgcComp, config, pairs0, pairs1, pairsComp, matches, isMatched0, isMatched1);
        interDim.computePairsAndMatch(ctr0, ctr1, ctrComp);
    }

    {
        Dimension0 dim0(cgc0, cgc1, cgcComp,  config, pairs0[0], pairs1[0], pairsComp[0], matches[0], isMatched0, isMatched1);       
        dim0.computePairsAndMatch(ctr0, ctr1, ctrComp);
    }

    if (config.print) {
        printResult(dim, cgc0, cgc1, cgcComp, pairs0, pairs1, pairsComp, matches, isMatched0, isMatched1);
    }
}
    