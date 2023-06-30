#include "template_functions.h"

#include "top_dimension.h"
#include "inter_dimensions.h"

#include "npy.hpp"

#include <iostream>
#include <fstream>
#include <cfloat>
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


void read_image(string const &filename, file_format const &format, vector<float> &image, vector<uint64_t> &shape) {
    switch(format){
        case DIPHA: {
            ifstream fin(filename, ios::in | ios::binary );
            int64_t d;
            fin.read((char *) &d, sizeof(int64_t));
            assert(d == 8067171840);
            fin.read((char *) &d, sizeof(int64_t));
            assert(d == 1);
            fin.read((char *) &d, sizeof(int64_t));
            fin.read((char *) &d, sizeof(int64_t));
            uint8_t dim = d;
            assert(dim < 4);
            uint64_t n;
            fin.read((char *) &d, sizeof(int64_t));
            shape.push_back(d);
            n = d;
            if (dim>1) {
                fin.read((char *) &d, sizeof(int64_t));
                shape.push_back(d);
                n *= d;
            }
            if (dim>2) {
                fin.read((char *)&d, sizeof(int64_t));
                shape.push_back(d);
                n *= d;
            }
            double value;
            image.reserve(n);
            while (!fin.eof()){
                fin.read((char *)&value, sizeof(double));
                image.push_back(value);
            }
            fin.close();
            return;
        }
        case PERSEUS: {
            ifstream reading_file; 
            reading_file.open(filename.c_str(), ios::in); 
            string reading_line_buffer; 
            getline(reading_file, reading_line_buffer); 
            uint8_t dim = atoi(reading_line_buffer.c_str());
            assert(dim < 4);
            uint64_t n;
            getline(reading_file, reading_line_buffer);
            shape.push_back(atoi(reading_line_buffer.c_str()));
            n = shape[0];
            if (dim>1) {
                getline(reading_file, reading_line_buffer); 
                shape.push_back(atoi(reading_line_buffer.c_str()));
                n *= shape[1];
            }
            if (dim>2) {
                getline(reading_file, reading_line_buffer);
                shape.push_back(atoi(reading_line_buffer.c_str()));
                n *= shape[2];
            }
            image.reserve(n);
            double value;
            while(!reading_file.eof()){
                getline(reading_file, reading_line_buffer);
                value = atof(reading_line_buffer.c_str());
                if (value != -1) {
                    image.push_back(value);
                }else{
                    image.push_back(DBL_MAX);
                }
            }
            reading_file.close();
            return;
		}
        case NUMPY: {
            vector<unsigned long> _shape;
            try{
                npy::LoadArrayFromNumpy(filename.c_str(), _shape, image);
            } catch (...) {
                cerr << "The data type of an numpy array should be numpy.float64." << endl;
                exit(-2);
            }
            uint8_t dim = shape.size();
            assert (dim < 4);
            for (uint32_t i : _shape) {
                shape.push_back(i);
            }
            return;
        }
    }
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
            config.min_recursion_to_cache = stoi(argv[++i]);
		} else if (arg == "--cache_size" || arg == "-c") {
            config.cache_size = stoi(argv[++i]);
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
    vector<float> image0;
    vector<float> image1;
    vector<uint64_t> shape0;
    vector<uint64_t> shape1;
    read_image(config.filename_0, config.format_0, image0, shape0);
    read_image(config.filename_1, config.format_1, image1, shape1);
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);
    auto duration_total = duration;
    if (config.verbose) { cout << "took " << duration.count() << " ms" << endl; }
    
    assert (shape0 == shape1);
    uint64_t dim = shape0.size();

    vector<float> imageComp;
    transform(image0.begin(), image0.end(), image1.begin(), back_inserter(imageComp), [](float a, float b){return min(a,b);});
    
    CubicalGridComplex cgc0(move(image0), shape0);
    CubicalGridComplex cgc1(move(image1), shape1);
    CubicalGridComplex cgcComp(move(imageComp), shape0);

    vector<Cube> ctrComp;
	vector<Cube> ctr0;
	vector<Cube> ctr1;

    if (dim > 0) {
        if (config.verbose) { cout << "comoputing top dimension ... "; }
        start = high_resolution_clock::now(); 
        TopDimension topDim(cgc0, cgc1, cgcComp, config);
        topDim.computePairsAndMatch(ctr0, ctr1, ctrComp);
        stop = high_resolution_clock::now();
        duration = duration_cast<milliseconds>(stop - start);
        duration_total += duration;
        if (config.verbose) { cout << "took " << duration.count() << " ms" << endl; }

        cout << "matches in topdim" << endl;
        for (auto& m : topDim.matches) {
            m.print(); cout << endl;
        }
    }

    if (dim > 1) {
        if (config.verbose) { cout << "comoputing top dimension ... "; }
        start = high_resolution_clock::now();
        InterDimensions interDim(cgc0, cgc1, cgcComp, config);
        interDim.computePairsAndMatch(ctr0, ctr1, ctrComp);
        stop = high_resolution_clock::now();
        duration = duration_cast<milliseconds>(stop - start);
        duration_total += duration;
        if (config.verbose) { cout << "took " << duration.count() << " ms" << endl; }

        cout << "matches in interdims" << endl;
        for (auto& m : interDim.matches[1]) {
            m.print();
        }
    
    }
}
    