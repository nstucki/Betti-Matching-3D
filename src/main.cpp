#include "BettiMatching.h"
#include "utils.h"

#include <algorithm>
#include <cassert>
#include <chrono>
#include <fstream>
#include <iostream>

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
         << "  --threshold, -t                 compute persistent homology of "
            "up to threshold"
         << endl
         << "  --print, -p                     print result to console" << endl
         << "  --save, -s                      save result in .txt files"
         << endl
         << "  --min_recursion_to_cache, -mc   minimum number of recursion for "
            "a reduced column to be cached"
         << endl
         << "  --cache_size, -c	          maximum number of reduced columns to "
            "be cached"
         << endl
         << "  --matched, -m                   name of the .txt file "
            "containing matched pairs"
         << endl
         << "  --unmatchedmatched_0, -u0       name of the .txt file "
            "containing unmatched pairs of input 0"
         << endl
         << "  --unmatchedmatched_1, -u1       name of the .txt file "
            "containing unmatched pairs of input 1"
         << endl
         << endl;
    exit(exit_code);
}

int main(int argc, char **argv) {
#ifdef RUNTIME
    cout << endl << "reading config & images ... ";
    auto start = high_resolution_clock::now();
#endif

    Config config;
    string filename0 = "";
    string filename1 = "";
    string matchedFilename = "matched.csv";
    string unmatched0Filename = "unmatched_0.csv";
    string unmatched1Filename = "unmatched_1.csv";
    fileFormat format0;
    fileFormat format1;
    bool print = false;
    bool saveResult = false;

    for (int i = 1; i < argc; ++i) {
        const string arg(argv[i]);
        if (arg == "--help" || arg == "-h") {
            print_usage_and_exit(0);
        } else if (arg == "--min_recursion_to_cache" || arg == "-mc") {
            config.minRecursionToCache = stoi(argv[++i]);
        } else if (arg == "--cache_size" || arg == "-c") {
            config.cacheSize = stoi(argv[++i]);
        } else if (arg == "--print" || arg == "-p") {
            print = true;
        } else if (arg == "--matched" || arg == "-m") {
            matchedFilename = string(argv[++i]);
        } else if (arg == "--unmatched_0" || arg == "-u0") {
            unmatched0Filename = string(argv[++i]);
        } else if (arg == "--unmatched_1" || arg == "-u1") {
            unmatched1Filename = string(argv[++i]);
        } else {
            if (filename0.empty()) {
                filename0 = argv[i];
            } else {
                filename1 = argv[i];
            }
        }
    }

    if (filename0.empty()) {
        print_usage_and_exit(-1);
    }
    if (filename0.find(".txt") != string::npos) {
        format0 = PERSEUS;
    } else if (filename0.find(".npy") != string::npos) {
        format0 = NUMPY;
    } else if (filename0.find(".complex") != string::npos) {
        format0 = DIPHA;
    } else {
        cerr << "unknown input file format! (the filename extension should be "
                ".txt/.npy/.complex): "
             << filename0 << endl;
        exit(-1);
    }
    ifstream fileStream0(filename0);
    if (!filename0.empty() && fileStream0.fail()) {
        cerr << "couldn't open file " << filename0 << endl;
        exit(-1);
    }
    if (filename1.empty()) {
        print_usage_and_exit(-1);
    }
    if (filename1.find(".txt") != string::npos) {
        format1 = PERSEUS;
    } else if (filename1.find(".npy") != string::npos) {
        format1 = NUMPY;
    } else if (filename1.find(".complex") != string::npos) {
        format1 = DIPHA;
    } else {
        cerr << "unknown input file format! (the filename extension should be "
                ".txt/.npy/.complex): "
             << filename1 << endl;
        exit(-1);
    }
    ifstream fileStream1(filename1);
    if (!filename1.empty() && fileStream1.fail()) {
        cerr << "couldn't open file " << filename1 << endl;
        exit(-1);
    }

    vector<value_t> input0;
    vector<value_t> input1;
    vector<value_t> comparison;
    vector<index_t> shape;

    {
        vector<index_t> shape0;
        vector<index_t> shape1;
        readImage(filename0, format0, input0, shape0);
        readImage(filename1, format1, input1, shape1);
        assert(shape0 == shape1);
        shape = shape0;
    }

#ifdef RUNTIME
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);
    cout << "of shape (" << shape[0];
    for (uint8_t i = 1; i < shape.size(); i++) {
        cout << "," << shape[i];
    }
    cout << ") ... " << duration.count() << " ms" << endl << endl;
#endif

    BettiMatching BM(std::move(input0), std::move(input1), std::move(shape),
                     std::move(config));
    BM.computeMatching();

    if (print) {
        BM.printResult();
    }
    if (saveResult) {
    }
}
