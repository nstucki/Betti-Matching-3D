#include "utils.h"
#include "dimension_0.h"
#include "dimension_1.h"
#include "dimension_2.h"

#include <iostream>
#include <sstream> 
#include <fstream>
#include <cfloat>
#include <chrono>
#include <algorithm>
#include <cassert>

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
         << "  --verbose, -v                   print result in console" << endl
         << "  --min_recursion_to_cache, -mc   minimum number of recursion for a reduced column to be cached" << endl
         << "  --cache_size, -c	          maximum number of reduced columns to be cached" << endl
         << "  --matched, -m                   name of the file containing matched pairs" << endl
         << "  --unmatchedmatched_0, -u0       name of the file containing unmatched pairs of input 0" << endl
         << "  --unmatchedmatched_1, -u1       name of the file containing unmatched pairs of input 1" << endl                
         << endl;
	exit(exit_code);
}


int main(int argc, char** argv) {
    #ifdef RUNTIME
    cout << endl <<  "reading config & images ... ";
    auto startTotal = high_resolution_clock::now();
    auto start = high_resolution_clock::now();
    #endif 
    Config config;
    for (int i = 1; i < argc; ++i) {
		const string arg(argv[i]);
		if (arg == "--help" || arg == "-h") { print_usage_and_exit(0); } 
        else if (arg == "--verbose" || arg == "-v") { config.verbose = true; } 
        else if (arg == "--matched" || arg == "-m") { config.matchedFilename = string(argv[++i]); } 
        else if (arg == "--unmatched_0" || arg == "-u0") { config.unmatched0Filename = string(argv[++i]); } 
        else if (arg == "--unmatched_1" || arg == "-u1") { config.unmatched1Filename = string(argv[++i]); } 
        else if (arg == "--min_recursion_to_cache" || arg == "-mc") { config.minRecursionToCache = stoi(argv[++i]); } 
        else if (arg == "--cache_size" || arg == "-c") { config.cacheSize = stoi(argv[++i]); } 
        else {
            if (config.filename0.empty()) { config.filename0 = argv[i]; } 
            else { config.filename1 = argv[i]; }
		}
	}
    if (config.filename0.empty()) { print_usage_and_exit(-1); } 
    if (config.filename0.find(".txt")!= std::string::npos) { config.format0 = PERSEUS; } 
    else if (config.filename0.find(".npy")!= std::string::npos) { config.format0 = NUMPY; } 
    else if(config.filename0.find(".complex")!= std::string::npos) { config.format0 = DIPHA; } 
    else {
		cerr << "unknown input file format! (the filename extension should be .txt/.npy/.complex): " << config.filename0 << endl;
		exit(-1);
	}
    ifstream fileStream0(config.filename0);
	if (!config.filename0.empty() && fileStream0.fail()) {
		cerr << "couldn't open file " << config.filename0 << endl;
		exit(-1);
	}
    if (config.filename1.empty()) { print_usage_and_exit(-1); }
    if (config.filename1.find(".txt")!= std::string::npos) { config.format1 = PERSEUS; } 
    else if (config.filename1.find(".npy")!= std::string::npos) { config.format1 = NUMPY; } 
    else if (config.filename1.find(".complex")!= std::string::npos) { config.format1 = DIPHA; } 
    else {
		cerr << "unknown input file format! (the filename extension should be .txt/.npy/.complex): " << config.filename1 << endl;
		exit(-1);
	}
    ifstream fileStream1(config.filename1);
	if (!config.filename1.empty() && fileStream1.fail()) {
		cerr << "couldn't open file " << config.filename1 << endl;
		exit(-1);
	}

    vector<value_t> image0;
    vector<value_t> image1;
    vector<value_t> imageComp;
    vector<index_t> shape;
    {   
        vector<index_t> shape0;
        vector<index_t> shape1;
        readImage(config.filename0, config.format0, image0, shape0);
        readImage(config.filename1, config.format1, image1, shape1);
        assert (shape0 == shape1);
        shape = shape0;
        transform(image0.begin(), image0.end(), image1.begin(), back_inserter(imageComp), [](value_t a, value_t b) 
        { return min(a, b); });
    }
    #ifdef RUNTIME
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);
    auto duration_total = duration;
    cout << "of shape (" << shape[0] << "," << shape[1] << "," << shape[2] << ") ... ";
    cout << duration.count() << " ms" << endl << endl;
    #endif

    #ifdef RUNTIME
    cout << "initializing CubicalGridComplex and results ... ";
    start = start = high_resolution_clock::now();
    #endif
    CubicalGridComplex* cgc0 = new CubicalGridComplex(std::move(image0), shape);
    CubicalGridComplex* cgc1 = new CubicalGridComplex(std::move(image1), shape);
    CubicalGridComplex* cgcComp = new CubicalGridComplex(std::move(imageComp), shape);

    vector<vector<Pair>> pairs0(3);
    vector<vector<Pair>> pairs1(3);
    vector<vector<Pair>> pairsComp(3);
    vector<vector<Match>> matches(3);
    vector<unordered_map<uint64_t, bool>> isMatched0(3);
	vector<unordered_map<uint64_t, bool>> isMatched1(3);

    vector<Cube> ctr0;
    vector<Cube> ctr1;
    vector<Cube> ctrComp;
    #ifdef RUNTIME
    stop = high_resolution_clock::now();
    duration = duration_cast<milliseconds>(stop - start);
    cout << duration.count() << " ms" << endl << endl;
    #endif

    {   
        #ifdef RUNTIME
        cout << "computing dimension 2 ... " << endl;
        start = high_resolution_clock::now();
        #endif
        Dimension2 dim2(cgc0, cgc1, cgcComp,  config, pairs0[2], pairs1[2], pairsComp[2], matches[2], isMatched0[2], isMatched1[2]);       
        dim2.computePairsAndMatch(ctr0, ctr1, ctrComp);
        #ifdef RUNTIME
        stop = high_resolution_clock::now();
        duration = duration_cast<milliseconds>(stop - start);
        cout << "total: " << duration.count() << " ms" << endl << endl;
        #endif
    }

    {   
        #ifdef RUNTIME
        start = high_resolution_clock::now();
        cout << "computing dimension 1 ... " << endl; 
        #endif
        Dimension1 dim1(cgc0, cgc1, cgcComp,  config, pairs0[1], pairs1[1], pairsComp[1], matches[1], isMatched0[1], isMatched1[1]);       
        dim1.computePairsAndMatch(ctr0, ctr1, ctrComp);
        #ifdef RUNTIME
        stop = high_resolution_clock::now();
        duration = duration_cast<milliseconds>(stop - start);
        cout << "total: " << duration.count() << " ms" << endl << endl;
        #endif
    }
    
    {   
        #ifdef RUNTIME
        start = high_resolution_clock::now();
        cout << "computing dimension 0 ... " << endl; 
        #endif
        Dimension0 dim0(cgc0, cgc1, cgcComp,  config, pairs0[0], pairs1[0], pairsComp[0], matches[0], isMatched0[0], isMatched1[0]);       
        dim0.computePairsAndMatch(ctr0, ctr1, ctrComp);
        #ifdef RUNTIME
        stop = high_resolution_clock::now();
        duration = duration_cast<milliseconds>(stop - start);
        cout << "total: " << duration.count() << " ms" << endl << endl;
        #endif
    }

    #ifdef RUNTIME
    stop = high_resolution_clock::now();
    duration = duration_cast<milliseconds>(stop - startTotal);
    cout << "Betti Matching runtime: " << duration.count() << " ms" << endl << endl;
    #endif

    if (config.verbose) { printResult(cgc0, cgc1, cgcComp, pairs0, pairs1, pairsComp, matches, isMatched0, isMatched1); }
}