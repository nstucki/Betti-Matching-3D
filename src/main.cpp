#include "utils.h"
#include "src_3D/BettiMatching.h"
#include "src_2D/BettiMatching.h"

#include <iostream>
#include <fstream>

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
         << "  --print, -p                     print result to console" << endl
         << "  --save, -s                      save result in .txt files" << endl
         << "  --min_recursion_to_cache, -mc   minimum number of recursion for a reduced column to be cached" << endl
         << "  --cache_size, -c	          maximum number of reduced columns to be cached" << endl
         << "  --matched, -m                   name of the .txt file containing matched pairs" << endl
         << "  --unmatchedmatched_0, -u0       name of the .txt file containing unmatched pairs of input 0" << endl
         << "  --unmatchedmatched_1, -u1       name of the .txt file containing unmatched pairs of input 1" << endl                
         << endl;
	exit(exit_code);
}


int main(int argc, char** argv) {
#ifdef RUNTIME
    cout << endl << "reading config & images ... ";
    auto startTotal = high_resolution_clock::now();
#endif 
    Config config;
    for (int i = 1; i < argc; ++i) {
		const string arg(argv[i]);
		if (arg == "--help" || arg == "-h") { print_usage_and_exit(0); }
        else if (arg == "--threshold" || arg == "-t") { config.threshold = stoi(argv[++i]); }
        else if (arg == "--save" || arg == "-s") { config.minRecursionToCache = stoi(argv[++i]); } 
        else if (arg == "--print" || arg == "-p") { config.print = true; } 
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
    vector<value_t> input0;
    vector<value_t> input1;
    vector<value_t> comparison;
    vector<index_t> shape;
    size_t dim;
    {   
        vector<index_t> shape0;
        vector<index_t> shape1;
        readImage(config.filename0, config.format0, input0, shape0);
        readImage(config.filename1, config.format1, input1, shape1);
        assert (shape0 == shape1);
        shape = shape0;
        dim = shape.size();
        transform(input0.begin(), input0.end(), input1.begin(), back_inserter(comparison), 
                    [](value_t a, value_t b) { return min(a, b); });
    }
#ifdef RUNTIME
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - startTotal);
    cout << "of shape (" << shape[0];
    for (uint8_t i = 1; i < shape.size(); i++) { cout << "," << shape[i]; }
    cout << ") ... " << duration.count() << " ms" << endl << endl;
#endif

    if (dim == 2) {
#ifdef RUNTIME
        cout << "initializing BettiMatching ... ";
        auto start = high_resolution_clock::now();
#endif
        dim2::BettiMatching BM{std::move(input0), std::move(input1), std::move(comparison), std::move(shape), config};
#ifdef RUNTIME
        stop = high_resolution_clock::now();
        duration = duration_cast<milliseconds>(stop - startTotal);
        cout << duration.count() << " ms" << endl << endl;
#endif

#ifdef RUNTIME
        cout << "computing Betti Matching ..." << endl;
#endif
        BM.computeMatching();
        BM.computeVoxels();
        
#ifdef RUNTIME
        stop = high_resolution_clock::now();
        duration = duration_cast<milliseconds>(stop - startTotal);
        cout << "total runtime: " << duration.count() << " ms" << endl << endl;
#endif

        if (config.print) { BM.printResult(); }
        if (config.saveResult) {}
    } 

    else if (dim == 3) {
#ifdef RUNTIME
        cout << "initializing BettiMatching ... ";
        auto start = high_resolution_clock::now();
#endif
        dim3::BettiMatching BM{std::move(input0), std::move(input1), std::move(comparison), std::move(shape), config};
#ifdef RUNTIME
        stop = high_resolution_clock::now();
        duration = duration_cast<milliseconds>(stop - startTotal);
        cout << duration.count() << " ms" << endl << endl;
#endif

#ifdef RUNTIME
        cout << "computing Betti Matching ..." << endl;
#endif
        BM.computeMatching();
        BM.computeVoxels();

#ifdef RUNTIME
        stop = high_resolution_clock::now();
        duration = duration_cast<milliseconds>(stop - startTotal);
        cout << "Betti Matching runtime: " << duration.count() << " ms" << endl << endl;
#endif

        if (config.print) { BM.printResult(); }
        if (config.saveResult) {}
    }
}